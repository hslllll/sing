use mimalloc::MiMalloc;
#[global_allocator]
static GLOBAL: MiMalloc = MiMalloc;

use anyhow::{bail, Result};
use clap::{Parser, Subcommand};
use crossbeam_channel::{bounded, Receiver, Sender};
use needletail::parse_fastx_file;
use sing_core::{
    align, build_index_from_sequences, cigar_ref_span, oriented_bases, write_cigar_string, AlignmentResult, Index,
    IndexLike, MemoryIndex, State,
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
            let total_threads = threads.unwrap_or(max_hw).max(2);
            
            let worker_threads = ((total_threads - 1) * 3) / 4;  
            let reader_threads = total_threads - 1 - worker_threads;  
            let reader_threads = reader_threads.max(1);  
            
            let batch_size = 4096;
            let batch_cap = worker_threads * 16;

            eprintln!("Loading index from {:?}...", index);
            let idx = Arc::new(load_index(&index)?);
            eprintln!("Index loaded.");

            let (tx, rx): (Sender<ReadBatch>, Receiver<ReadBatch>) = bounded(batch_cap);
            let (recycle_tx, recycle_rx): (Sender<ReadBatch>, Receiver<ReadBatch>) = bounded(batch_cap);
            let (w_tx, w_rx): (Sender<Vec<u8>>, Receiver<Vec<u8>>) = bounded(batch_cap);
            let paired_mode = r2.is_some();
            let mut reader_handles = Vec::new();

            eprintln!("Mapping started with {} threads (readers: {}, workers: {}, writer: 1)...", 
                      total_threads, reader_threads, worker_threads);
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
            let max_write_merge = std::cmp::min(std::cmp::max(worker_threads * 1024 * 1024, 8 * 1024 * 1024), 64 * 1024 * 1024);
            let writer = thread::spawn(move || {
                let mut out: Box<dyn Write> = match output {
                    Some(p) => Box::new(BufWriter::with_capacity(16 * 1024 * 1024, File::create(p).unwrap())),
                    None => Box::new(BufWriter::with_capacity(16 * 1024 * 1024, std::io::stdout())),
                };

                writeln!(out, "@HD\tVN:1.6\tSO:unsorted").unwrap();
                for i in 0..idx_writer.ref_count() {
                    let name = idx_writer.ref_name(i);
                    let seq = idx_writer.ref_seq(i);
                    writeln!(out, "@SQ\tSN:{}\tLN:{}", name, seq.len()).unwrap();
                }
                writeln!(out, "@RG\tID:sing\tSM:song\tPL:ILLUMINA").unwrap();
                writeln!(out, "@PG\tID:{}\tPN:{}\tVN:{}", PROG_NAME, PROG_NAME, PROG_VERSION).unwrap();

                while let Ok(mut chunk) = w_rx.recv() {
                    while let Ok(next) = w_rx.try_recv() {
                        if chunk.len() + next.len() > max_write_merge {
                            break;
                        }
                        chunk.extend_from_slice(&next);
                    }
                    out.write_all(&chunk).unwrap();
                }
            });

            let mut handles = Vec::new();
            for _ in 0..worker_threads {
                let rx = rx.clone();
                let w_tx = w_tx.clone();
                let recycle_tx = recycle_tx.clone();
                let idx = idx.clone();
                let total_reads = total_reads.clone();
                let mapped_reads = mapped_reads.clone();
                handles.push(thread::spawn(move || {
                    let mut state = State::new();
                    let mut rev_buf = Vec::new();
                    while let Ok(batch) = rx.recv() {
                        let mut output_chunk = Vec::with_capacity(batch.count * 800);
                        for i in 0..batch.count {
                            let s1 = &batch.seqs1[i];
                            let q1 = &batch.quals1[i];
                            let name = &batch.names[i];

                            let a1 = align(s1, &idx, &mut state, &mut rev_buf);
                            let is_paired = !batch.seqs2[i].is_empty();
                            let a2 = if is_paired {
                                align(&batch.seqs2[i], &idx, &mut state, &mut rev_buf)
                            } else {
                                None
                            };

                            total_reads.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                            if a1.is_some() {
                                mapped_reads.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                            }
                            if is_paired {
                                total_reads.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
                                if a2.is_some() {
                                    mapped_reads.fetch_add(1, std::sync::atomic::Ordering::Relaxed);
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

                                        if res1.is_rev != res2.is_rev && dist.abs() < 2000 {
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
                        try_send_retry(&w_tx, output_chunk);
                        try_send_retry(&recycle_tx, batch);
                    }
                }));
            }

            for h in reader_handles {
                h.join().unwrap();
            }
            drop(tx);

            for h in handles {
                h.join().unwrap();
            }
            drop(w_tx);
            writer.join().unwrap();

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
