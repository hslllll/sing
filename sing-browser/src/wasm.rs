use wasm_bindgen::prelude::*;

use anyhow::{anyhow, Result};
use flate2::read::GzDecoder;
use needletail::parse_fastx_reader;
use noodles::{bam, sam, bgzf};
use noodles::sam::alignment::io::Write as _;
use sing_core::{
    align, build_index_from_sequences, oriented_bases, write_cigar_string, AlignmentResult, Index, State,
};
use std::io::{BufRead, BufReader, Cursor, Write};

#[wasm_bindgen(getter_with_clone)]
pub struct MappingResult {
    pub bam_data: Vec<u8>,
    pub mapped_reads: usize,
    pub mapped_bases: usize,
    pub total_reads: usize,
    pub mapped_mask: Vec<u8>,
}

#[wasm_bindgen]
pub struct SingWebEngine {
    index: Option<Index>,
    state: State,
    rev_buf: Vec<u8>,
    total_reads: usize,
    mapped_reads: usize,
    mapped_bases: usize,
    genome_size: usize,
}

#[wasm_bindgen]
impl SingWebEngine {
    #[wasm_bindgen(constructor)]
    pub fn new() -> SingWebEngine {
        SingWebEngine {
            index: None,
            state: State::new(),
            rev_buf: Vec::new(),
            total_reads: 0,
            mapped_reads: 0,
            mapped_bases: 0,
            genome_size: 0,
        }
    }

    #[wasm_bindgen]
    pub fn process_reference_chunk(&mut self, chunk_data: &[u8], _chunk_id: usize) -> Result<(), JsValue> {
        let records = parse_reference_chunk(chunk_data).map_err(to_js)?;
        let idx = build_index_from_sequences(records).map_err(to_js)?;
        self.genome_size = idx.ref_seqs.iter().map(|s| s.len()).sum();
        self.index = Some(idx);
        Ok(())
    }

    #[wasm_bindgen]
    pub fn export_index(&self) -> Result<Vec<u8>, JsValue> {
        if let Some(idx) = &self.index {
            let mut buf = Vec::new();
            idx.to_writer(&mut buf).map_err(to_js)?;
            Ok(buf)
        } else {
            Err(JsValue::from_str("No index loaded"))
        }
    }

    #[wasm_bindgen]
    pub fn import_index(&mut self, data: &[u8]) -> Result<(), JsValue> {
        let idx = Index::from_bytes(data).map_err(to_js)?;
        self.genome_size = idx.ref_seqs.iter().map(|s| s.len()).sum();
        self.index = Some(idx);
        Ok(())
    }

    #[wasm_bindgen]
    pub fn run_mapping_chunk(
        &mut self,
        ref_index_ptr: *const Index,
        read1_chunk: &[u8],
        read2_chunk: Option<Box<[u8]>>,
        output_format: String,
        sort_output: bool,
        write_header: bool,
    ) -> Result<MappingResult, JsValue> {
        let idx = resolve_index(ref_index_ptr, self.index.as_ref()).map_err(to_js)?;
        let header = header_from_index(idx, sort_output).map_err(to_js)?;
        let format = OutputFormat::parse(&output_format).map_err(to_js)?;

        let mut strategy: Box<dyn OutputStrategy> = if sort_output {
            Box::new(SortWriter::new(format, header.clone(), write_header))
        } else {
            Box::new(StreamWriter::new(format, header.clone(), write_header).map_err(to_js)?)
        };

        let reads1 = parse_reads_chunk(read1_chunk).map_err(to_js)?;
        
        let chunk_start_total = self.total_reads;
        let chunk_start_mapped = self.mapped_reads;
        let chunk_start_bases = self.mapped_bases;

        if let Some(read2) = read2_chunk {
            let reads2 = parse_reads_chunk(&read2).map_err(to_js)?;
            if reads1.len() != reads2.len() {
                return Err(JsValue::from_str("paired-end inputs have mismatched read counts"));
            }

            for ((name1, seq1, qual1), (name2, seq2, qual2)) in reads1.into_iter().zip(reads2.into_iter()) {
                let res1_candidates = align(&seq1, idx, &mut self.state, &mut self.rev_buf);
                let res2_candidates = align(&seq2, idx, &mut self.state, &mut self.rev_buf);
                let res1 = res1_candidates.first().cloned();
                let res2 = res2_candidates.first().cloned();

                self.total_reads += 2;
                if res1.is_some() { self.mapped_reads += 1; self.mapped_bases += seq1.len(); }
                if res2.is_some() { self.mapped_reads += 1; self.mapped_bases += seq2.len(); }

                let sam1 = sam_record(&name1, &seq1, &qual1, &res1, Some(&res2), true, idx).map_err(to_js)?;
                let sam2 = sam_record(&name2, &seq2, &qual2, &res2, Some(&res1), false, idx).map_err(to_js)?;

                strategy.write(&sam1).map_err(to_js)?;
                strategy.write(&sam2).map_err(to_js)?;
            }
        } else {
            for (name, seq, qual) in reads1 {
                let res_candidates = align(&seq, idx, &mut self.state, &mut self.rev_buf);
                let res = res_candidates.first().cloned();
                self.total_reads += 1;
                if res.is_some() { self.mapped_reads += 1; self.mapped_bases += seq.len(); }

                let sam_rec = sam_record(&name, &seq, &qual, &res, None, true, idx).map_err(to_js)?;
                strategy.write(&sam_rec).map_err(to_js)?;
            }
        }

        let MappingResultInternal { bam_data } = strategy.finalize().map_err(to_js)?;
        
        // Create a mapped_mask - 1 byte per read, indicating if it was mapped
        let chunk_total_reads = self.total_reads - chunk_start_total;
        let chunk_mapped_reads = self.mapped_reads - chunk_start_mapped;
        let chunk_mapped_bases = self.mapped_bases - chunk_start_bases;
        
        let mask_bytes = (chunk_total_reads + 7) / 8;
        let mut mapped_mask = vec![0u8; mask_bytes];
        
        // Mark the reads that were mapped in this chunk
        for i in 0..chunk_mapped_reads {
            let byte_idx = i / 8;
            let bit_idx = i % 8;
            if byte_idx < mapped_mask.len() {
                mapped_mask[byte_idx] |= 1 << bit_idx;
            }
        }
        
        Ok(MappingResult {
            bam_data,
            mapped_reads: chunk_mapped_reads,
            mapped_bases: chunk_mapped_bases,
            total_reads: chunk_total_reads,
            mapped_mask,
        })
    }

    #[wasm_bindgen]
    pub fn get_total_reads(&self) -> usize {
        self.total_reads
    }

    #[wasm_bindgen]
    pub fn get_mapped_reads(&self) -> usize {
        self.mapped_reads
    }

    #[wasm_bindgen]
    pub fn get_mapped_bases(&self) -> usize {
        self.mapped_bases
    }

    #[wasm_bindgen]
    pub fn get_genome_size(&self) -> usize {
        self.genome_size
    }

    #[wasm_bindgen]
    pub fn get_avg_coverage(&self) -> f64 {
        if self.genome_size == 0 {
            0.0
        } else {
            self.mapped_bases as f64 / self.genome_size as f64
        }
    }

    #[wasm_bindgen]
    pub fn map_reads_chunk(&mut self, reads_chunk: &[u8]) -> Result<MappingResult, JsValue> {
        self.run_mapping_chunk(std::ptr::null(), reads_chunk, None, "bam".to_string(), false, true)
    }
}

fn to_js<E: std::fmt::Display>(err: E) -> JsValue {
    JsValue::from_str(&err.to_string())
}

fn parse_reference_chunk(data: &[u8]) -> Result<Vec<(String, Vec<u8>)>> {
    let reader = get_smart_reader(data);
    let mut parser = parse_fastx_reader(reader)?;
    let mut out = Vec::new();
    while let Some(rec) = parser.next() {
        let rec = rec?;
        let name = String::from_utf8_lossy(rec.id())
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        out.push((name, rec.seq().to_vec()));
    }
    Ok(out)
}

fn parse_reads_chunk(data: &[u8]) -> Result<Vec<(String, Vec<u8>, Vec<u8>)>> {
    let reader = get_smart_reader(data);
    let mut parser = parse_fastx_reader(reader)?;
    let mut out = Vec::new();
    while let Some(rec) = parser.next() {
        let rec = rec?;
        let name = String::from_utf8_lossy(rec.id())
            .split_whitespace()
            .next()
            .unwrap_or("")
            .to_string();
        let seq = rec.seq().to_vec();
        let qual = rec
            .qual()
            .map(|q| q.to_vec())
            .unwrap_or_else(|| vec![b'*'; seq.len().max(1)]);
        out.push((name, seq, qual));
    }
    Ok(out)
}

fn header_from_index(idx: &Index, sorted: bool) -> Result<sam::Header> {
    let sort_tag = if sorted { "coordinate" } else { "unsorted" };
    let mut header_txt = format!("@HD\tVN:1.6\tSO:{}\n", sort_tag);
    for (name, seq) in idx.ref_names.iter().zip(idx.ref_seqs.iter()) {
        header_txt.push_str(&format!("@SQ\tSN:{}\tLN:{}\n", name, seq.len()));
    }
    header_txt.push_str("@RG\tID:sing-web\tSM:web\tPL:ILLUMINA\n");
    header_txt.push_str(&format!("@PG\tID:sing-web\tPN:sing-web\tVN:{}\n", env!("CARGO_PKG_VERSION")));
    Ok(header_txt.parse()?)
}

fn sam_record(
    name: &str,
    seq: &[u8],
    qual: &[u8],
    res: &Option<AlignmentResult>,
    mate: Option<&Option<AlignmentResult>>,
    is_first: bool,
    idx: &Index,
) -> Result<sam::Record> {
    let sam_line_str = sam_line(name, seq, qual, res, mate, is_first, idx);
    let mut reader = sam::io::Reader::new(Cursor::new(sam_line_str.as_bytes()));
    let mut sam_record = sam::Record::default();
    reader.read_record(&mut sam_record)?;
    Ok(sam_record)
}

fn sam_line(
    name: &str,
    seq: &[u8],
    qual: &[u8],
    res: &Option<AlignmentResult>,
    mate: Option<&Option<AlignmentResult>>,
    is_first: bool,
    idx: &Index,
) -> String {
    let mut out = Vec::new();
    format_single_sam(&mut out, name, seq, qual, res, mate, is_first, idx);
    String::from_utf8(out).unwrap_or_default()
}

fn format_single_sam(
    out: &mut Vec<u8>,
    name: &str,
    seq: &[u8],
    qual: &[u8],
    res: &Option<AlignmentResult>,
    mate: Option<&Option<AlignmentResult>>,
    is_first: bool,
    idx: &Index,
) {
    out.extend_from_slice(name.as_bytes());
    out.push(b'\t');
    let mut flag = 0u16;
    let (final_seq, final_qual) = oriented_bases(seq, qual, res);
    if res.is_none() {
        flag |= 0x4;
    }
    if let Some(r) = res {
        if r.is_rev {
            flag |= 0x10;
        }
    }

    let is_paired = mate.is_some();
    if is_paired {
        flag |= 0x1;
        if is_first {
            flag |= 0x40;
        } else {
            flag |= 0x80;
        }

        match mate {
            Some(Some(m)) => {
                if m.is_rev {
                    flag |= 0x20;
                }
            }
            _ => flag |= 0x8,
        }
    }

    let mut rname = "*";
    let mut pos = 0;
    let mut mapq = 0;
    let mut cigar_ref: Option<&Vec<u32>> = None;
    if let Some(r) = res {
        rname = &idx.ref_names[r.ref_id as usize];
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

    if let Some(Some(m)) = mate {
        let mate_rname = &idx.ref_names[m.ref_id as usize];
        write!(out, "\t{}\t{}\t0\t", mate_rname, m.pos + 1).unwrap();
    } else {
        write!(out, "\t*\t0\t0\t").unwrap();
    }

    out.extend_from_slice(&final_seq);
    out.push(b'\t');
    out.extend_from_slice(&final_qual);
    if let Some(r) = res {
        write!(out, "\tNM:i:{}\tAS:i:{}", r.nm, r.as_score).unwrap();
        if !r.md.is_empty() {
            write!(out, "\tMD:Z:{}", r.md).unwrap();
        }
    }
    out.push(b'\n');
}

fn get_smart_reader<'a>(data: &'a [u8]) -> Box<dyn BufRead + Send + 'a> {
    if data.len() >= 2 && data[0] == 0x1f && data[1] == 0x8b {
        let decoder = GzDecoder::new(Cursor::new(data));
        Box::new(BufReader::new(decoder))
    } else {
        Box::new(BufReader::new(Cursor::new(data)))
    }
}

fn resolve_index<'a>(ptr: *const Index, fallback: Option<&'a Index>) -> Result<&'a Index> {
    if !ptr.is_null() {
        return unsafe { ptr.as_ref().ok_or_else(|| anyhow!("reference pointer was null")) };
    }
    fallback.ok_or_else(|| anyhow!("reference chunk not loaded"))
}

#[derive(Clone, Copy)]
enum OutputFormat {
    Bam,
    Sam,
}

impl OutputFormat {
    fn parse(fmt: &str) -> Result<Self> {
        match fmt.to_ascii_lowercase().as_str() {
            "bam" => Ok(OutputFormat::Bam),
            "sam" => Ok(OutputFormat::Sam),
            other => Err(anyhow!("unsupported output format: {}", other)),
        }
    }
}

struct MappingResultInternal {
    bam_data: Vec<u8>,
}

trait OutputStrategy {
    fn write(&mut self, record: &sam::Record) -> Result<()>;
    fn finalize(self: Box<Self>) -> Result<MappingResultInternal>;
}

type BamWriter = bam::io::Writer<bgzf::io::writer::Writer<Cursor<Vec<u8>>>>;

struct StreamWriter {
    header: sam::Header,
    writer: StreamWriterKind,
}

enum StreamWriterKind {
    Bam(BamWriter),
    Sam(sam::io::Writer<Cursor<Vec<u8>>>),
}

impl StreamWriter {
    fn new(format: OutputFormat, header: sam::Header, write_header: bool) -> Result<Self> {
        let cursor = Cursor::new(Vec::new());
        let writer = match format {
            OutputFormat::Bam => {
                let mut w = bam::io::Writer::new(cursor);
                if write_header {
                    w.write_header(&header)?;
                }
                StreamWriterKind::Bam(w)
            }
            OutputFormat::Sam => {
                let mut w = sam::io::Writer::new(cursor);
                if write_header {
                    w.write_header(&header)?;
                }
                StreamWriterKind::Sam(w)
            }
        };
        Ok(StreamWriter { header, writer })
    }
}

impl OutputStrategy for StreamWriter {
    fn write(&mut self, record: &sam::Record) -> Result<()> {
        match &mut self.writer {
            StreamWriterKind::Bam(w) => w.write_alignment_record(&self.header, record)?,
            StreamWriterKind::Sam(w) => w.write_record(&self.header, record)?,
        }
        Ok(())
    }

    fn finalize(self: Box<Self>) -> Result<MappingResultInternal> {
        let bam_data = match self.writer {
            StreamWriterKind::Bam(w) => w.into_inner().finish()?.into_inner(),
            StreamWriterKind::Sam(w) => w.into_inner().into_inner(),
        };
        Ok(MappingResultInternal { bam_data })
    }
}

struct SortWriter {
    header: sam::Header,
    format: OutputFormat,
    records: Vec<sam::Record>,
    write_header: bool,
}

impl SortWriter {
    fn new(format: OutputFormat, header: sam::Header, write_header: bool) -> Self {
        Self { header, format, records: Vec::new(), write_header }
    }
}

impl OutputStrategy for SortWriter {
    fn write(&mut self, record: &sam::Record) -> Result<()> {
        self.records.push(record.clone());
        Ok(())
    }

    fn finalize(mut self: Box<Self>) -> Result<MappingResultInternal> {
        self.records.sort_by(|a, b| {
            let refs = self.header.reference_sequences();

            let a_ref = a
                .reference_sequence_name()
                .and_then(|name| refs.get_index_of(name.as_ref() as &[u8]).map(|i| i as u32));
            let b_ref = b
                .reference_sequence_name()
                .and_then(|name| refs.get_index_of(name.as_ref() as &[u8]).map(|i| i as u32));

            let a_key = a_ref.unwrap_or(u32::MAX);
            let b_key = b_ref.unwrap_or(u32::MAX);

            a_key.cmp(&b_key).then_with(|| {
                let a_pos = a
                    .alignment_start()
                    .and_then(|p| p.ok())
                    .map(|p| p.get() as u32)
                    .unwrap_or(u32::MAX);
                let b_pos = b
                    .alignment_start()
                    .and_then(|p| p.ok())
                    .map(|p| p.get() as u32)
                    .unwrap_or(u32::MAX);
                a_pos.cmp(&b_pos)
            })
        });

        match self.format {
             OutputFormat::Bam => {
                let cursor = Cursor::new(Vec::new());
                let mut writer = bam::io::Writer::new(cursor);
                if self.write_header {
                    writer.write_header(&self.header)?;
                }

                for record in &self.records {
                    writer.write_alignment_record(&self.header, record)?;
                }

                let bam_data = writer.into_inner().finish()?.into_inner();

                     Ok(MappingResultInternal { bam_data })
             },
             OutputFormat::Sam => {
                let cursor = Cursor::new(Vec::new());
                let mut writer = sam::io::Writer::new(cursor);
                if self.write_header {
                    writer.write_header(&self.header)?;
                }
                for record in &self.records {
                     writer.write_record(&self.header, record)?;
                }
                Ok(MappingResultInternal { bam_data: writer.into_inner().into_inner() })
             }
        }
    }
}
