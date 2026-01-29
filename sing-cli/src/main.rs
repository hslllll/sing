use mimalloc::MiMalloc;
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use anyhow::{bail, Result};
use clap::{Parser, Subcommand};
use crossbeam_channel::{bounded, Receiver, Sender};
use needletail::parse_fastx_file;
use sing_core::{
    align, build_index_from_sequences, cigar_ref_span, oriented_bases, write_cigar_string, AlignmentResult, Index,
    IndexLike, MemoryIndex, State, CONFIG,
};
use std::fs::File;
use std::io::{BufWriter, Write};
use std::path::PathBuf;
use std::sync::Arc;
use std::thread;

const PROG_NAME: &str = env!("CARGO_PKG_NAME");
const PROG_VERSION: &str = env!("CARGO_PKG_VERSION");
const MAX_SEQ_LEN: usize = 600;

#[derive(Parser)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    Build { reference: PathBuf, output: PathBuf },
    Map {
        index: PathBuf,
        #[arg(short = '1', long = "r1", num_args = 1.., value_name = "READ1", required = true)]
        r1: Vec<PathBuf>,
        #[arg(short = '2', long = "r2", num_args = 1.., value_name = "READ2")]
        r2: Option<Vec<PathBuf>>,
        #[arg(short, long)]
        output: Option<PathBuf>,
        #[arg(short, long)]
        threads: Option<usize>,
    },
}

struct ReadBatch {
    seqs1: Vec<Vec<u8>>,
    quals1: Vec<Vec<u8>>,
    seqs2: Vec<Vec<u8>>,
    quals2: Vec<Vec<u8>>,
    names: Vec<String>,
    count: usize,
    max_size: usize,
}

impl ReadBatch {
    fn new(max_size: usize) -> Self {
        Self {
            seqs1: Vec::with_capacity(max_size),
            quals1: Vec::with_capacity(max_size),
            seqs2: Vec::with_capacity(max_size),
            quals2: Vec::with_capacity(max_size),
            names: Vec::with_capacity(max_size),
            count: 0,
            max_size,
        }
    }

    fn clear(&mut self) {
        self.count = 0;
    }
}

fn load_index(path: &PathBuf) -> Result<MemoryIndex> {
    MemoryIndex::from_path(path)
}

fn save_index(index: &Index, path: &PathBuf) -> Result<()> {
    let file = File::create(path)?;
    let mut writer = BufWriter::with_capacity(16*1024*1024, file);
    index.to_writer(&mut writer)
}

fn build_from_reference(reference: PathBuf, output: PathBuf) -> Result<()> {
    eprintln!("Reading reference sequences from {:?}...", reference);
    let mut reader = parse_fastx_file(&reference)?;
    let mut records = Vec::new();
    while let Some(r) = reader.next() {
        let r = r?;
        let name = String::from_utf8_lossy(r.id())
            .split_whitespace()
            .next()
            .unwrap()
            .to_string();
        records.push((name, r.seq().to_vec()));
    }

    eprintln!("Building index...");
    let index = build_index_from_sequences(records)?;
    eprintln!("Saving index to {:?}...", output);
    save_index(&index, &output)?;
    eprintln!("Index build completed.");
    Ok(())
}

fn try_send_retry<T>(sender: &Sender<T>, value: T) {
    let _ = sender.send(value);
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Commands::Build { reference, output } => {
            build_from_reference(reference, output)?;
        }
        Commands::Map { index, r1, r2, output, threads } => {
            if r1.is_empty() {
                bail!("Provide at least one -1 input file");
            }

            let max_hw = std::thread::available_parallelism().map(|n| n.get()).unwrap_or(4);
            let worker_threads = threads.unwrap_or(max_hw).max(1);
            
            let reader_threads = (r1.len()*2).max(1);
            
            let batch_size = reader_threads * 512;

            eprintln!("Loading index from {:?}...", index);
            let idx = Arc::new(load_index(&index)?);
            eprintln!("Index loaded.");

            let channel_cap = worker_threads * 4;
            let (tx, rx): (Sender<ReadBatch>, Receiver<ReadBatch>) = bounded(channel_cap);
            let (recycle_tx, recycle_rx): (Sender<ReadBatch>, Receiver<ReadBatch>) = bounded(channel_cap);
            let paired_mode = r2.is_some();
            let mut reader_handles = Vec::new();

            eprintln!("Mapping with {} worker threads...", worker_threads);
            let total_reads = Arc::new(std::sync::atomic::AtomicUsize::new(0));
            let mapped_reads = Arc::new(std::sync::atomic::AtomicUsize::new(0));

            if paired_mode {
                let r2 = r2.unwrap();
                if r1.len() != r2.len() {
                    bail!("The number of -1 and -2 files must match");
                }

                let pairs: Vec<(PathBuf, PathBuf)> = r1.into_iter().zip(r2.into_iter()).collect();
                let pair_idx = Arc::new(std::sync::atomic::AtomicUsize::new(0));

                for _ in 0..reader_threads {
                    let pairs = pairs.clone();
                    let pair_idx = pair_idx.clone();
                    let tx = tx.clone();
                    let recycle_rx = recycle_rx.clone();

                    reader_handles.push(thread::spawn(move || {
                        let mut batch = recycle_rx
                            .try_recv()
                            .unwrap_or_else(|_| ReadBatch::new(batch_size));
                        batch.clear();
                        loop {
                            let i = pair_idx.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                            if i >= pairs.len() {
                                break;
                            }
                            let (r1_path, r2_path) = &pairs[i];

                            let mut r1 = parse_fastx_file(r1_path).unwrap();
                            let mut r2 = parse_fastx_file(r2_path).unwrap();

                            while let (Some(Ok(rec1)), Some(Ok(rec2))) = (r1.next(), r2.next()) {
                                if batch.count >= batch.seqs1.len() {
                                    batch.seqs1.push(Vec::with_capacity(MAX_SEQ_LEN));
                                    batch.quals1.push(Vec::with_capacity(MAX_SEQ_LEN));
                                    batch.seqs2.push(Vec::with_capacity(MAX_SEQ_LEN));
                                    batch.quals2.push(Vec::with_capacity(MAX_SEQ_LEN));
                                    batch.names.push(String::with_capacity(100));
                                }
                                let idx = batch.count;
                                let name_str = String::from_utf8_lossy(rec1.id());
                                let name_part = name_str.split_whitespace().next().unwrap();
                                batch.names[idx].clear();
                                batch.names[idx].push_str(name_part);
                                if batch.names[idx].ends_with("/1") || batch.names[idx].ends_with("/2") {
                                    let new_len = batch.names[idx].len() - 2;
                                    batch.names[idx].truncate(new_len);
                                }

                                batch.seqs1[idx].clear();
                                batch.seqs1[idx].extend_from_slice(&rec1.seq());
                                batch.quals1[idx].clear();
                                if let Some(q) = rec1.qual() {
                                    batch.quals1[idx].extend_from_slice(q);
                                } else {
                                    batch.quals1[idx].push(b'*');
                                }

                                batch.seqs2[idx].clear();
                                batch.seqs2[idx].extend_from_slice(&rec2.seq());
                                batch.quals2[idx].clear();
                                if let Some(q) = rec2.qual() {
                                    batch.quals2[idx].extend_from_slice(q);
                                } else {
                                    batch.quals2[idx].push(b'*');
                                }

                                batch.count += 1;

                                if batch.count >= batch.max_size {
                                    try_send_retry(&tx, batch);
                                    batch = recycle_rx
                                        .try_recv()
                                        .unwrap_or_else(|_| ReadBatch::new(batch_size));
                                    batch.clear();
                                }
                            }
                        }
                        if batch.count > 0 {
                            try_send_retry(&tx, batch);
                        }
                    }));
                }
            } else {
                let r1_paths = r1.clone();
                let r1_idx = Arc::new(std::sync::atomic::AtomicUsize::new(0));
                for _ in 0..reader_threads {
                    let r1_paths = r1_paths.clone();
                    let r1_idx = r1_idx.clone();
                    let tx = tx.clone();
                    let recycle_rx = recycle_rx.clone();

                    reader_handles.push(thread::spawn(move || {
                        let mut batch = recycle_rx
                            .try_recv()
                            .unwrap_or_else(|_| ReadBatch::new(batch_size));
                        batch.clear();
                        loop {
                            let i = r1_idx.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                            if i >= r1_paths.len() {
                                break;
                            }
                            let r1_path = &r1_paths[i];
                            let mut r1 = parse_fastx_file(r1_path).unwrap();

                            while let Some(Ok(rec1)) = r1.next() {
                                if batch.count >= batch.seqs1.len() {
                                    batch.seqs1.push(Vec::with_capacity(MAX_SEQ_LEN));
                                    batch.quals1.push(Vec::with_capacity(MAX_SEQ_LEN));
                                    batch.seqs2.push(Vec::new());
                                    batch.quals2.push(Vec::new());
                                    batch.names.push(String::with_capacity(64));
                                }
                                let idx = batch.count;
                                let name_str = String::from_utf8_lossy(rec1.id());
                                let name_part = name_str.split_whitespace().next().unwrap();
                                batch.names[idx].clear();
                                batch.names[idx].push_str(name_part);
                                if batch.names[idx].ends_with("/1") || batch.names[idx].ends_with("/2") {
                                    let new_len = batch.names[idx].len() - 2;
                                    batch.names[idx].truncate(new_len);
                                }

                                batch.seqs1[idx].clear();
                                batch.seqs1[idx].extend_from_slice(&rec1.seq());
                                batch.quals1[idx].clear();
                                if let Some(q) = rec1.qual() {
                                    batch.quals1[idx].extend_from_slice(q);
                                } else {
                                    batch.quals1[idx].push(b'*');
                                }

                                batch.seqs2[idx].clear();
                                batch.quals2[idx].clear();

                                batch.count += 1;

                                if batch.count >= batch.max_size {
                                    try_send_retry(&tx, batch);
                                    batch = recycle_rx
                                        .try_recv()
                                        .unwrap_or_else(|_| ReadBatch::new(batch_size));
                                    batch.clear();
                                }
                            }
                        }
                        if batch.count > 0 {
                            try_send_retry(&tx, batch);
                        }
                    }));
                }
            }

            let idx_writer = idx.clone();
            
            let mut writer: Box<dyn Write + Send> = match &output {
                Some(p) => Box::new(BufWriter::with_capacity(4 * 1024 * 1024, File::create(p)?)),
                None => Box::new(BufWriter::with_capacity(2 * 1024 * 1024, std::io::stdout())),
            };
            
            writeln!(writer, "@HD\tVN:1.6\tSO:unsorted")?;
            for i in 0..idx_writer.ref_count() {
                let name = idx_writer.ref_name(i);
                let seq = idx_writer.ref_seq(i);
                writeln!(writer, "@SQ\tSN:{}\tLN:{}", name, seq.len())?;
            }
            writeln!(writer, "@RG\tID:sing\tSM:song\tPL:ILLUMINA")?;
            writeln!(writer, "@PG\tID:{}\tPN:{}\tVN:{}", PROG_NAME, PROG_NAME, PROG_VERSION)?;

            let writer_queue_cap = worker_threads * 64;
            let (wtx, wrx): (Sender<Vec<u8>>, Receiver<Vec<u8>>) = bounded(writer_queue_cap);
            let (recycle_out_tx, recycle_out_rx): (Sender<Vec<u8>>, Receiver<Vec<u8>>) = bounded(writer_queue_cap);
            let writer_handle = thread::spawn(move || {
                let mut out = writer;
                while let Ok(mut chunk) = wrx.recv() {
                    out.write_all(&chunk).unwrap();
                    chunk.clear();
                    let _ = recycle_out_tx.try_send(chunk);
                }
                out.flush().unwrap();
            });

            let mut handles = Vec::new();
            for _ in 0..worker_threads {
                let rx = rx.clone();
                let recycle_tx = recycle_tx.clone();
                let idx = idx.clone();
                let total_reads = total_reads.clone();
                let mapped_reads = mapped_reads.clone();
                let wtx = wtx.clone();
                let recycle_out_rx = recycle_out_rx.clone();
                
                handles.push(thread::spawn(move || {
                    let mut state = State::new();
                    let mut rev_buf = Vec::new();
                    while let Ok(batch) = rx.recv() {
                        let mut output_chunk = match recycle_out_rx.try_recv() {
                            Ok(mut chunk) => {
                                chunk.clear();
                                chunk
                            }
                            Err(_) => Vec::with_capacity(batch.count * 800),
                        };
                        let mut batch_total = 0usize;
                        let mut batch_mapped = 0usize;
                        for i in 0..batch.count {
                            let s1 = &batch.seqs1[i];
                            let q1 = &batch.quals1[i];
                            let name = &batch.names[i];

                            let a1_candidates = align(s1, &idx, &mut state, &mut rev_buf);
                            let is_paired = !batch.seqs2[i].is_empty();
                            let a2_candidates = if is_paired {
                                align(&batch.seqs2[i], &idx, &mut state, &mut rev_buf)
                            } else {
                                Vec::new()
                            };

                            let (final_a1, final_a2) = if is_paired {
                                let mut valid_pairs = Vec::with_capacity(4);

                                for (i, res1) in a1_candidates.iter().enumerate() {
                                    for (j, res2) in a2_candidates.iter().enumerate() {
                                        if res1.ref_id == res2.ref_id && res1.is_rev != res2.is_rev {
                                            let span1 = cigar_ref_span(&res1.cigar);
                                            let span2 = cigar_ref_span(&res2.cigar);
                                            let end1 = res1.pos + span1;
                                            let end2 = res2.pos + span2;
                                            let left = res1.pos.min(res2.pos);
                                            let right = end1.max(end2);
                                            let dist = right - left;

                                            if dist < CONFIG.pair_max_dist {
                                                let score = res1.as_score + res2.as_score;
                                                valid_pairs.push((score, i, j));
                                            }
                                        }
                                    }
                                }

                                if !valid_pairs.is_empty() {
                                    valid_pairs.sort_by(|a, b| b.0.cmp(&a.0));
                                    let (best_score, i, j) = valid_pairs[0];
                                    let mut best_res1 = a1_candidates[i].clone();
                                    let mut best_res2 = a2_candidates[j].clone();

                                    // Set paired flag (0x1) and properly matched pair flag (0x2)
                                    let mapq = if valid_pairs.len() > 1 {
                                        let (second_score, _, _) = valid_pairs[1];
                                        let diff = best_score - second_score;
                                        if diff > 30 { 60 } else { (diff * 2) as u8 }
                                    } else {
                                        60
                                    };
                                    
                                    best_res1.mapq = mapq;
                                    best_res2.mapq = mapq;
                                    best_res1.paired = true;
                                    best_res2.paired = true;
                                    best_res1.proper_pair = true;
                                    best_res2.proper_pair = true;
                                    
                                    (Some(best_res1), Some(best_res2))
                                } else {
                                    // Paired-end mode requires concordant pair
                                    (None, None)
                                }
                            } else {
                                (a1_candidates.first().cloned(), None)
                            };

                            let a1 = final_a1;
                            let a2 = final_a2;

                            batch_total += 1;
                            if a1.is_some() {
                                batch_mapped += 1;
                            }
                            if is_paired {
                                batch_total += 1;
                                if a2.is_some() {
                                    batch_mapped += 1;
                                }
                            }

                            if is_paired {
                                let mut tlen = 0;
                                let mut proper_pair = false;

                                if let (Some(res1), Some(res2)) = (&a1, a2.as_ref()) {
                                    if res1.ref_id == res2.ref_id {
                                        let span1 = cigar_ref_span(&res1.cigar);
                                        let span2 = cigar_ref_span(&res2.cigar);
                                        let end1 = res1.pos + span1;
                                        let end2 = res2.pos + span2;
                                        let left = res1.pos.min(res2.pos);
                                        let right = end1.max(end2);
                                        let dist = right - left;

                                        if res1.pos <= res2.pos {
                                            tlen = dist;
                                        } else {
                                            tlen = -dist;
                                        }

                                        if res1.is_rev != res2.is_rev && dist.abs() < CONFIG.pair_max_dist {
                                            proper_pair = true;
                                        }
                                    }
                                }

                                let q2 = &batch.quals2[i];
                                let s2 = &batch.seqs2[i];

                                format_sam(&mut output_chunk, name, s1, q1, &a1, &a2, true, tlen, proper_pair, &idx);
                                format_sam(
                                    &mut output_chunk,
                                    name,
                                    s2,
                                    q2,
                                    &a2,
                                    &a1,
                                    false,
                                    -tlen,
                                    proper_pair,
                                    &idx,
                                );
                            } else {
                                format_sam_single(&mut output_chunk, name, s1, q1, &a1, &idx);
                            }
                        }
                        
                        
                        if batch_total > 0 {
                            total_reads.fetch_add(batch_total, std::sync::atomic::Ordering::Relaxed);
                        }
                        if batch_mapped > 0 {
                            mapped_reads.fetch_add(batch_mapped, std::sync::atomic::Ordering::Relaxed);
                        }

                        let _ = wtx.send(output_chunk);
                        
                        try_send_retry(&recycle_tx, batch);
                    }
                }));
            }

            for h in reader_handles {
                h.join().unwrap();
            }
            drop(tx);
            drop(recycle_rx);  

            for h in handles {
                h.join().unwrap();
            }
            drop(recycle_tx);  
            drop(wtx);
            writer_handle.join().unwrap();

            let total = total_reads.load(std::sync::atomic::Ordering::Relaxed);
            let mapped = mapped_reads.load(std::sync::atomic::Ordering::Relaxed);
            let rate = if total > 0 {
                (mapped as f64 / total as f64) * 100.0
            } else {
                0.0
            };
            eprintln!("Total reads: {}", total);
            eprintln!("Mapped reads: {} ({:.2}%)", mapped, rate);
        }
    }
    std::process::exit(0);
}

fn format_sam<I: IndexLike>(
    out: &mut Vec<u8>,
    name: &str,
    seq: &[u8],
    qual: &[u8],
    res: &Option<AlignmentResult>,
    mate_res: &Option<AlignmentResult>,
    first: bool,
    tlen: i32,
    proper_pair: bool,
    idx: &I,
) {
    out.extend_from_slice(name.as_bytes());
    out.push(b'\t');
    let mut flag = 0x1;
    if first {
        flag |= 0x40;
    } else {
        flag |= 0x80;
    }
    let valid_res = res.as_ref().filter(|r| (r.ref_id as usize) < idx.ref_count());
    let valid_mate = mate_res.as_ref().filter(|m| (m.ref_id as usize) < idx.ref_count());
    let is_mapped = valid_res.is_some();
    let is_mate_mapped = valid_mate.is_some();
    let (final_seq, final_qual) = oriented_bases(seq, qual, &valid_res.cloned());
    if let Some(r) = valid_res {
        if r.is_rev {
            flag |= 0x10;
        }
    } else {
        flag |= 0x4;
    }
    if let Some(m) = valid_mate {
        if m.is_rev {
            flag |= 0x20;
        }
    } else {
        flag |= 0x8;
    }
    if proper_pair && is_mapped && is_mate_mapped {
        flag |= 0x2;
    }
    // Also check the proper_pair flag from the result itself
    if let Some(r) = valid_res {
        if r.proper_pair {
            flag |= 0x2;
        }
    }
    let mut rname = "*";
    let mut pos = 0;
    let mut mapq = 0;
    let mut cigar_ref: Option<&Vec<u32>> = None;
    if let Some(r) = valid_res {
        rname = idx.ref_name(r.ref_id as usize);
        pos = r.pos + 1;
        mapq = r.mapq;
        cigar_ref = Some(&r.cigar);
    }
    let mut rnext = "*";
    let mut pnext = 0;
    if let Some(m) = valid_mate {
        pnext = m.pos + 1;
        if let Some(r) = valid_res {
            if m.ref_id == r.ref_id {
                rnext = "=";
            } else {
                rnext = idx.ref_name(m.ref_id as usize);
            }
        } else {
            rnext = idx.ref_name(m.ref_id as usize);
        }
    }
    write!(out, "{}\t{}\t{}\t{}\t", flag, rname, pos, mapq).unwrap();
    if let Some(c) = cigar_ref {
        write_cigar_string(c, out);
    } else {
        out.push(b'*');
    }
    write!(out, "\t{}\t{}\t{}\t", rnext, pnext, tlen).unwrap();
    out.extend_from_slice(&final_seq);
    out.push(b'\t');
    out.extend_from_slice(&final_qual);
    if let Some(r) = valid_res {
        write!(out, "\tNM:i:{}\tAS:i:{}", r.nm, r.as_score).unwrap();
        if !r.md.is_empty() {
            write!(out, "\tMD:Z:{}", r.md).unwrap();
        }
    }
    write!(out, "\tRG:Z:sing").unwrap();
    out.push(b'\n');
}

fn format_sam_single<I: IndexLike>(
    out: &mut Vec<u8>,
    name: &str,
    seq: &[u8],
    qual: &[u8],
    res: &Option<AlignmentResult>,
    idx: &I,
) {
    out.extend_from_slice(name.as_bytes());
    out.push(b'\t');
    let mut flag = 0u16;
    let valid_res = res.as_ref().filter(|r| (r.ref_id as usize) < idx.ref_count());
    if valid_res.is_none() {
        flag |= 0x4;
    }
    if let Some(r) = valid_res {
        if r.is_rev {
            flag |= 0x10;
        }
    }
    let (final_seq, final_qual) = oriented_bases(seq, qual, &valid_res.cloned());
    let mut rname = "*";
    let mut pos = 0;
    let mut mapq = 0;
    let mut cigar_ref: Option<&Vec<u32>> = None;
    if let Some(r) = valid_res {
        rname = idx.ref_name(r.ref_id as usize);
        pos = r.pos + 1;
        mapq = r.mapq;
        cigar_ref = Some(&r.cigar);
    }
    write!(out, "{}\t{}\t{}\t{}\t", flag, rname, pos, mapq).unwrap();
    if let Some(c) = cigar_ref {
        write_cigar_string(c, out);
    } else {
        out.push(b'*');
    }
    write!(out, "\t*\t0\t0\t").unwrap();
    out.extend_from_slice(&final_seq);
    out.push(b'\t');
    out.extend_from_slice(&final_qual);
    if let Some(r) = valid_res {
        write!(out, "\tNM:i:{}\tAS:i:{}", r.nm, r.as_score).unwrap();
        if !r.md.is_empty() {
            write!(out, "\tMD:Z:{}", r.md).unwrap();
        }
    }
    write!(out, "\tRG:Z:sing").unwrap();
    out.push(b'\n');
}
