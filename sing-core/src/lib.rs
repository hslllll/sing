use anyhow::{bail, Context, Result};
use bytemuck::try_cast_slice;
use memmap2::MmapOptions;
use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read, Write};
use std::ops::Range;
use std::path::Path;

// SIMD support detection - disabled for WASM targets
#[cfg(all(target_arch = "x86_64", not(target_family = "wasm")))]
use std::arch::x86_64::*;

#[cfg(all(target_arch = "x86_64", not(target_family = "wasm")))]
fn has_avx2() -> bool {
    is_x86_feature_detected!("avx2")
}

#[cfg(not(all(target_arch = "x86_64", not(target_family = "wasm"))))]
#[allow(dead_code)]
fn has_avx2() -> bool { false }

#[derive(Clone, Copy)]
pub struct Config {
    pub hash_window: usize,
    pub sync_s: usize,
    pub match_score: i32,
    pub mismatch_pen: i32,
    pub gap_open: i32,
    pub gap_ext: i32,
    pub x_drop: i32,
    pub max_seed_occ: usize,
    pub max_candidates: usize,
    pub maxindel: usize,
    pub pair_max_dist: i32,
    pub min_seeds: usize,
}

pub static CONFIG: Config = Config {
    hash_window: 21,
    sync_s: 6,
    match_score: 2,
    mismatch_pen: -2,
    gap_open: -2,
    gap_ext: -1,
    x_drop: 15,
    max_seed_occ: 100, // maximum occurrences of a seed to be used in genome
    max_candidates: 5000, // maximum number of candidate hits to consider per read
    pair_max_dist: 1000,
    maxindel: 15,
    min_seeds: 2,
};

const HASH_WINDOW: usize = CONFIG.hash_window;
const SYNC_S: usize = CONFIG.sync_s;

const BASES: [u64; 4] = [
    (std::f64::consts::E - 2f64).to_bits(),
    (std::f64::consts::PI - 3f64).to_bits(),
    (std::f64::consts::SQRT_2 - 1f64).to_bits(),
    (1.7320508075688772f64 - 1f64).to_bits(),
];

const BASE_LUT: [i8; 256] = {
    let mut lut = [(-1i8); 256];
    lut[b'A' as usize] = 0;
    lut[b'a' as usize] = 0;
    lut[b'C' as usize] = 1;
    lut[b'c' as usize] = 1;
    lut[b'G' as usize] = 2;
    lut[b'g' as usize] = 2;
    lut[b'T' as usize] = 3;
    lut[b't' as usize] = 3;
    lut
};

const ROT: u32 = 1;

const fn rot(val: u64, n: u32) -> u64 {
    (val << n) | (val >> (64 - n))
}

const REMOVE: [u64; 4] = [
    rot(BASES[0], (HASH_WINDOW as u32 * ROT) % 64),
    rot(BASES[1], (HASH_WINDOW as u32 * ROT) % 64),
    rot(BASES[2], (HASH_WINDOW as u32 * ROT) % 64),
    rot(BASES[3], (HASH_WINDOW as u32 * ROT) % 64),
];

const REMOVE_S: [u64; 4] = [
    rot(BASES[0], (SYNC_S as u32 * ROT) % 64),
    rot(BASES[1], (SYNC_S as u32 * ROT) % 64),
    rot(BASES[2], (SYNC_S as u32 * ROT) % 64),
    rot(BASES[3], (SYNC_S as u32 * ROT) % 64),
];

// Seed layout (u64): [hash_remainder: 24 bits] | [ref_id: 12 bits] | [pos: 28 bits]
pub const RADIX: usize = 20;
pub const SHIFT: usize = 64 - RADIX;
pub const HASH_BITS: usize = 24;
pub const RID_BITS: usize = 8;
pub const POS_BITS: usize = 32;
pub const RID_MASK: u64 = (1 << RID_BITS) - 1;
pub const POS_MASK: u64 = (1 << POS_BITS) - 1;
pub const HASH_MASK: u64 = (1 << HASH_BITS) - 1;

#[derive(Debug, Clone)]
pub struct Index {
    pub offsets: Vec<u32>,
    pub seeds: Vec<u64>,
    pub genome_size: u64,
    pub max_hits: u32,
    pub total_seed_occurrences: u64,
    pub seed_counts: Vec<u64>,
    pub ref_seqs: Vec<Vec<u8>>,
    pub ref_names: Vec<String>,
}

pub trait IndexLike {
    fn offsets(&self) -> &[u32];
    fn seeds(&self) -> &[u64];
    fn genome_size(&self) -> u64;
    fn max_hits(&self) -> usize;
    fn total_seed_occurrences(&self) -> u64;
    fn seed_count(&self, hash: u64) -> u32;
    fn ref_count(&self) -> usize;
    fn ref_seq(&self, id: usize) -> &[u8];
    fn ref_name(&self, id: usize) -> &str;
}

impl Index {
    pub fn from_bytes(data: &[u8]) -> Result<Self> {
        let mut pos = 0usize;

        fn read_u64(data: &[u8], pos: &mut usize, label: &str) -> Result<u64> {
            if *pos + 8 > data.len() {
                bail!("{}: unexpected EOF", label);
            }
            let val = u64::from_le_bytes(data[*pos..*pos + 8].try_into().unwrap());
            *pos += 8;
            Ok(val)
        }

        let offsets_len = read_u64(data, &mut pos, "read offsets len")? as usize;
        let offsets_bytes = offsets_len
            .checked_mul(4)
            .and_then(|n| pos.checked_add(n))
            .filter(|&end| end <= data.len())
            .ok_or_else(|| anyhow::anyhow!("read offsets: unexpected EOF"))?;
        let mut offsets = Vec::with_capacity(offsets_len);
        for chunk in data[pos..offsets_bytes].chunks_exact(4) {
            offsets.push(u32::from_le_bytes(chunk.try_into().unwrap()));
        }
        pos = offsets_bytes;
        let padding = (8 - (pos % 8)) % 8;
        let padded = pos.checked_add(padding).filter(|&end| end <= data.len())
            .ok_or_else(|| anyhow::anyhow!("read offsets padding: unexpected EOF"))?;
        pos = padded;

        let seeds_len = read_u64(data, &mut pos, "read seeds len")? as usize;
        let seeds_bytes = seeds_len
            .checked_mul(8)
            .and_then(|n| pos.checked_add(n))
            .filter(|&end| end <= data.len())
            .ok_or_else(|| anyhow::anyhow!("read seeds: unexpected EOF"))?;
        let mut seeds = Vec::with_capacity(seeds_len);
        for chunk in data[pos..seeds_bytes].chunks_exact(8) {
            seeds.push(u64::from_le_bytes(chunk.try_into().unwrap()));
        }
        pos = seeds_bytes;

        let genome_size = read_u64(data, &mut pos, "read genome_size")?;

        let ref_count = read_u64(data, &mut pos, "read ref seq count")? as usize;
        let mut ref_seqs = Vec::with_capacity(ref_count);
        for _ in 0..ref_count {
            let slen = read_u64(data, &mut pos, "read seq len")? as usize;
            let end = pos.checked_add(slen).filter(|&e| e <= data.len())
                .ok_or_else(|| anyhow::anyhow!("read seq bytes: unexpected EOF"))?;
            ref_seqs.push(data[pos..end].to_vec());
            pos = end;
        }

        let name_count = read_u64(data, &mut pos, "read name count")? as usize;
        let mut ref_names = Vec::with_capacity(name_count);
        for _ in 0..name_count {
            let slen = read_u64(data, &mut pos, "read name len")? as usize;
            let end = pos.checked_add(slen).filter(|&e| e <= data.len())
                .ok_or_else(|| anyhow::anyhow!("read name bytes: unexpected EOF"))?;
            let name_buf = &data[pos..end];
            ref_names.push(String::from_utf8(name_buf.to_vec()).context("utf8 ref name")?);
            pos = end;
        }

        let max_hits = if pos + 8 <= data.len() {
            read_u64(data, &mut pos, "read max_hits")? as u32
        } else {
            u32::MAX
        };

        let (total_seed_occurrences, seed_counts) = if pos + 16 <= data.len() {
            let total = read_u64(data, &mut pos, "read total_seed_occurrences")?;
            let counts_len = read_u64(data, &mut pos, "read seed_counts len")? as usize;
            let counts_bytes = counts_len
                .checked_mul(8)
                .and_then(|n| pos.checked_add(n))
                .filter(|&end| end <= data.len())
                .ok_or_else(|| anyhow::anyhow!("read seed_counts: unexpected EOF"))?;
            let mut counts = Vec::with_capacity(counts_len);
            for chunk in data[pos..counts_bytes].chunks_exact(8) {
                counts.push(u64::from_le_bytes(chunk.try_into().unwrap()));
            }
            (total, counts)
        } else {
            (0, Vec::new())
        };

        Ok(Index {
            offsets,
            seeds,
            genome_size,
            max_hits,
            total_seed_occurrences,
            seed_counts,
            ref_seqs,
            ref_names,
        })
    }

    pub fn from_reader<R: Read>(reader: R) -> Result<Self> {
        let mut r = BufReader::new(reader);
        let mut buf8 = [0u8; 8];
        let mut buf4 = [0u8; 4];

        r.read_exact(&mut buf8).context("read offsets len")?;
        let len = u64::from_le_bytes(buf8) as usize;
        let mut offsets = Vec::with_capacity(len);
        for _ in 0..len {
            r.read_exact(&mut buf4).context("read offset")?;
            offsets.push(u32::from_le_bytes(buf4));
        }
        let pos = 8usize + len.saturating_mul(4);
        let padding = (8 - (pos % 8)) % 8;
        if padding != 0 {
            let mut pad = [0u8; 8];
            r.read_exact(&mut pad[..padding]).context("read offsets padding")?;
        }

        r.read_exact(&mut buf8).context("read seeds len")?;
        let len = u64::from_le_bytes(buf8) as usize;
        let mut seeds = Vec::with_capacity(len);
        for _ in 0..len {
            r.read_exact(&mut buf8).context("read seed")?;
            seeds.push(u64::from_le_bytes(buf8));
        }

        r.read_exact(&mut buf8).context("read genome_size")?;
        let genome_size = u64::from_le_bytes(buf8);

        r.read_exact(&mut buf8).context("read ref seq count")?;
        let len = u64::from_le_bytes(buf8) as usize;
        let mut ref_seqs = Vec::with_capacity(len);
        for _ in 0..len {
            r.read_exact(&mut buf8).context("read seq len")?;
            let slen = u64::from_le_bytes(buf8) as usize;
            let mut seq = vec![0u8; slen];
            r.read_exact(&mut seq).context("read seq bytes")?;
            ref_seqs.push(seq);
        }

        r.read_exact(&mut buf8).context("read name count")?;
        let len = u64::from_le_bytes(buf8) as usize;
        let mut ref_names = Vec::with_capacity(len);
        for _ in 0..len {
            r.read_exact(&mut buf8).context("read name len")?;
            let slen = u64::from_le_bytes(buf8) as usize;
            let mut name_buf = vec![0u8; slen];
            r.read_exact(&mut name_buf).context("read name bytes")?;
            ref_names.push(String::from_utf8(name_buf).context("utf8 ref name")?);
        }

        let max_hits = {
            let mut buf = [0u8; 8];
            match r.read(&mut buf).context("read max_hits")? {
                0 => u32::MAX,
                n if n < 8 => {
                    r.read_exact(&mut buf[n..]).context("read max_hits")?;
                    u64::from_le_bytes(buf) as u32
                }
                _ => u64::from_le_bytes(buf) as u32,
            }
        };

        let (total_seed_occurrences, seed_counts) = {
            let mut buf = [0u8; 8];
            let total = match r.read(&mut buf).context("read total_seed_occurrences")? {
                0 => 0,
                n if n < 8 => {
                    r.read_exact(&mut buf[n..]).context("read total_seed_occurrences")?;
                    u64::from_le_bytes(buf)
                }
                _ => u64::from_le_bytes(buf),
            };
            if total == 0 {
                (0, Vec::new())
            } else {
                r.read_exact(&mut buf).context("read seed_counts len")?;
                let counts_len = u64::from_le_bytes(buf) as usize;
                let mut counts = Vec::with_capacity(counts_len);
                for _ in 0..counts_len {
                    r.read_exact(&mut buf).context("read seed_counts entry")?;
                    counts.push(u64::from_le_bytes(buf));
                }
                (total, counts)
            }
        };

        Ok(Index {
            offsets,
            seeds,
            genome_size,
            max_hits,
            total_seed_occurrences,
            seed_counts,
            ref_seqs,
            ref_names,
        })
    }

    pub fn to_writer<W: Write>(&self, mut w: W) -> Result<()> {
        w.write_all(&(self.offsets.len() as u64).to_le_bytes())?;
        for v in &self.offsets {
            w.write_all(&v.to_le_bytes())?;
        }
        let pos = 8usize + self.offsets.len().saturating_mul(4);
        let padding = (8 - (pos % 8)) % 8;
        if padding != 0 {
            let pad = [0u8; 8];
            w.write_all(&pad[..padding])?;
        }

        w.write_all(&(self.seeds.len() as u64).to_le_bytes())?;
        for s in &self.seeds {
            w.write_all(&s.to_le_bytes())?;
        }

        w.write_all(&self.genome_size.to_le_bytes())?;

        w.write_all(&(self.ref_seqs.len() as u64).to_le_bytes())?;
        for seq in &self.ref_seqs {
            w.write_all(&(seq.len() as u64).to_le_bytes())?;
            w.write_all(seq)?;
        }

        w.write_all(&(self.ref_names.len() as u64).to_le_bytes())?;
        for name in &self.ref_names {
            let bytes = name.as_bytes();
            w.write_all(&(bytes.len() as u64).to_le_bytes())?;
            w.write_all(bytes)?;
        }

        w.write_all(&(self.max_hits as u64).to_le_bytes())?;
        w.write_all(&self.total_seed_occurrences.to_le_bytes())?;
        w.write_all(&(self.seed_counts.len() as u64).to_le_bytes())?;
        for v in &self.seed_counts {
            w.write_all(&v.to_le_bytes())?;
        }

        Ok(())
    }
}

impl IndexLike for Index {
    fn offsets(&self) -> &[u32] { &self.offsets }
    fn seeds(&self) -> &[u64] { &self.seeds }
    fn genome_size(&self) -> u64 { self.genome_size }
    fn max_hits(&self) -> usize { self.max_hits as usize }
    fn total_seed_occurrences(&self) -> u64 { self.total_seed_occurrences }

    fn seed_count(&self, hash: u64) -> u32 {
        let key = hash >> 32;
        match self.seed_counts.binary_search_by_key(&key, |v| v >> 32) {
            Ok(idx) => (self.seed_counts[idx] & 0xFFFF_FFFF) as u32,
            Err(_) => 0,
        }
    }

    fn ref_count(&self) -> usize { self.ref_seqs.len() }
    fn ref_seq(&self, id: usize) -> &[u8] { &self.ref_seqs[id] }
    fn ref_name(&self, id: usize) -> &str { &self.ref_names[id] }
}

enum IndexBuffer {
    Mmap(memmap2::Mmap),
}

impl IndexBuffer {
    #[inline]
    fn as_slice(&self) -> &[u8] {
        match self { Self::Mmap(mmap) => mmap.as_ref() }
    }
}

pub struct MemoryIndex {
    buf: IndexBuffer,
    offsets_ptr: *const u32,
    offsets_len: usize,
    seeds_ptr: *const u64,
    seeds_len: usize,
    genome_size: u64,
    max_hits: u32,
    total_seed_occurrences: u64,
    seed_counts: Vec<u64>,
    ref_seq_ranges: Vec<Range<usize>>,
    ref_name_ranges: Vec<Range<usize>>,
}

unsafe impl Send for MemoryIndex {}
unsafe impl Sync for MemoryIndex {}

impl MemoryIndex {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe {
            MmapOptions::new().populate().map(&file).context("mmap index file")?
        };
        Self::from_data(IndexBuffer::Mmap(mmap))
    }

    fn from_data(buf: IndexBuffer) -> Result<Self> {
        let data = buf.as_slice();
        let mut pos = 0usize;

        fn read_u64(data: &[u8], pos: &mut usize, label: &str) -> Result<u64> {
            if *pos + 8 > data.len() { bail!("{}: unexpected EOF", label); }
            let val = u64::from_le_bytes(data[*pos..*pos + 8].try_into().unwrap());
            *pos += 8;
            Ok(val)
        }

        let offsets_len = read_u64(data, &mut pos, "read offsets len")? as usize;
        let offsets_bytes = offsets_len.checked_mul(4).and_then(|n| pos.checked_add(n))
            .filter(|&end| end <= data.len())
            .ok_or_else(|| anyhow::anyhow!("read offsets: unexpected EOF"))?;
        let offsets_slice = try_cast_slice::<u8, u32>(&data[pos..offsets_bytes])
            .map_err(|_| anyhow::anyhow!("read offsets: unaligned data"))?;
        let offsets_ptr = offsets_slice.as_ptr();
        let offsets_len = offsets_slice.len();
        pos = offsets_bytes;

        let seeds_len = read_u64(data, &mut pos, "read seeds len")? as usize;
        let seeds_bytes = seeds_len.checked_mul(8).and_then(|n| pos.checked_add(n))
            .filter(|&end| end <= data.len())
            .ok_or_else(|| anyhow::anyhow!("read seeds: unexpected EOF"))?;
        let seeds_slice = try_cast_slice::<u8, u64>(&data[pos..seeds_bytes])
            .map_err(|_| anyhow::anyhow!("read seeds: unaligned data"))?;
        let seeds_ptr = seeds_slice.as_ptr();
        let seeds_len = seeds_slice.len();
        pos = seeds_bytes;

        let genome_size = read_u64(data, &mut pos, "read genome_size")?;

        let ref_count = read_u64(data, &mut pos, "read ref seq count")? as usize;
        let mut ref_seq_ranges = Vec::with_capacity(ref_count);
        for _ in 0..ref_count {
            let slen = read_u64(data, &mut pos, "read seq len")? as usize;
            let end = pos.checked_add(slen).filter(|&e| e <= data.len())
                .ok_or_else(|| anyhow::anyhow!("read seq bytes: unexpected EOF"))?;
            ref_seq_ranges.push(pos..end);
            pos = end;
        }

        let name_count = read_u64(data, &mut pos, "read name count")? as usize;
        let mut ref_name_ranges = Vec::with_capacity(name_count);
        for _ in 0..name_count {
            let slen = read_u64(data, &mut pos, "read name len")? as usize;
            let end = pos.checked_add(slen).filter(|&e| e <= data.len())
                .ok_or_else(|| anyhow::anyhow!("read name bytes: unexpected EOF"))?;
            std::str::from_utf8(&data[pos..end]).context("utf8 ref name")?;
            ref_name_ranges.push(pos..end);
            pos = end;
        }

        let max_hits = if pos + 8 <= data.len() {
            read_u64(data, &mut pos, "read max_hits")? as u32
        } else { u32::MAX };

        let (total_seed_occurrences, seed_counts) = if pos + 16 <= data.len() {
            let total = read_u64(data, &mut pos, "read total_seed_occurrences")?;
            let counts_len = read_u64(data, &mut pos, "read seed_counts len")? as usize;
            let counts_bytes = counts_len.checked_mul(8).and_then(|n| pos.checked_add(n))
                .filter(|&end| end <= data.len())
                .ok_or_else(|| anyhow::anyhow!("read seed_counts: unexpected EOF"))?;
            let mut counts = Vec::with_capacity(counts_len);
            for chunk in data[pos..counts_bytes].chunks_exact(8) {
                counts.push(u64::from_le_bytes(chunk.try_into().unwrap()));
            }
            (total, counts)
        } else { (0, Vec::new()) };

        Ok(Self { buf, offsets_ptr, offsets_len, seeds_ptr, seeds_len, genome_size, max_hits,
            total_seed_occurrences, seed_counts, ref_seq_ranges, ref_name_ranges })
    }
}

impl IndexLike for MemoryIndex {
    fn offsets(&self) -> &[u32] { unsafe { std::slice::from_raw_parts(self.offsets_ptr, self.offsets_len) } }
    fn seeds(&self) -> &[u64] { unsafe { std::slice::from_raw_parts(self.seeds_ptr, self.seeds_len) } }
    fn genome_size(&self) -> u64 { self.genome_size }
    fn max_hits(&self) -> usize { self.max_hits as usize }
    fn total_seed_occurrences(&self) -> u64 { self.total_seed_occurrences }
    fn seed_count(&self, hash: u64) -> u32 {
        if self.seed_counts.is_empty() { return 0; }
        let key = hash >> 32;
        match self.seed_counts.binary_search_by_key(&key, |v| v >> 32) {
            Ok(idx) => (self.seed_counts[idx] & 0xFFFF_FFFF) as u32,
            Err(_) => 0,
        }
    }
    fn ref_count(&self) -> usize { self.ref_seq_ranges.len() }
    fn ref_seq(&self, id: usize) -> &[u8] { &self.buf.as_slice()[self.ref_seq_ranges[id].clone()] }
    fn ref_name(&self, id: usize) -> &str {
        unsafe { std::str::from_utf8_unchecked(&self.buf.as_slice()[self.ref_name_ranges[id].clone()]) }
    }
}

impl<T: IndexLike> IndexLike for std::sync::Arc<T> {
    fn offsets(&self) -> &[u32] { self.as_ref().offsets() }
    fn seeds(&self) -> &[u64] { self.as_ref().seeds() }
    fn genome_size(&self) -> u64 { self.as_ref().genome_size() }
    fn max_hits(&self) -> usize { self.as_ref().max_hits() }
    fn total_seed_occurrences(&self) -> u64 { self.as_ref().total_seed_occurrences() }
    fn seed_count(&self, hash: u64) -> u32 { self.as_ref().seed_count(hash) }
    fn ref_count(&self) -> usize { self.as_ref().ref_count() }
    fn ref_seq(&self, id: usize) -> &[u8] { self.as_ref().ref_seq(id) }
    fn ref_name(&self, id: usize) -> &str { self.as_ref().ref_name(id) }
}

#[derive(Clone, Default, Debug)]
pub struct AlignmentResult {
    pub ref_id: i32,
    pub pos: i32,
    pub is_rev: bool,
    pub mapq: u8,
    pub cigar: Vec<u32>,
    pub nm: i32,
    pub md: String,
    pub as_score: i32,
    pub paired: bool,
    pub proper_pair: bool,
}

pub fn build_index_from_sequences<I, S>(records: I) -> Result<Index>
where I: IntoIterator<Item = (S, Vec<u8>)>, S: Into<String>,
{
    let mut ref_seqs = Vec::new();
    let mut ref_names = Vec::new();
    let mut total_bases: u64 = 0;
    let mut total_seed_occurrences: u64 = 0;
    let mut mins = Vec::with_capacity(1024);
    let mut counts: HashMap<u64, u32> = HashMap::new();

    eprintln!("Building streaming index...");

    for (name, seq) in records.into_iter() {
        ref_names.push(name.into());
        total_bases += seq.len() as u64;
        get_syncmers(&seq, &mut mins);
        total_seed_occurrences = total_seed_occurrences.saturating_add(mins.len() as u64);
        for &(h, _p) in mins.iter() {
            *counts.entry(h).or_insert(0) = counts.get(&h).unwrap_or(&0).saturating_add(1);
        }
        ref_seqs.push(seq);
    }

    eprintln!("  Loaded {} bp across {} references", total_bases, ref_seqs.len());
    eprintln!("  Counted {} unique seeds", counts.len());

    let max_hits = CONFIG.max_seed_occ;

    eprintln!("  Total seed occurrences: {}, max_hits={}", total_seed_occurrences, max_hits);

    let mut bucket_counts = vec![0u32; 1 << RADIX];
    for seq in ref_seqs.iter() {
        get_syncmers(seq, &mut mins);
        for &(h, _p) in mins.iter() {
            let bucket = (h >> SHIFT) as usize;
            if bucket < bucket_counts.len() {
                bucket_counts[bucket] = bucket_counts[bucket].saturating_add(1);
            }
        }
    }

    let mut offsets = vec![0u32; 1 << RADIX];
    let mut running = 0u32;
    for i in 0..offsets.len() {
        offsets[i] = running;
        running = running.saturating_add(bucket_counts[i]);
    }

    let mut final_seeds: Vec<u64> = vec![0u64; running as usize];
    let mut write_offsets = offsets.clone();

    for (rid, seq) in ref_seqs.iter().enumerate() {
        get_syncmers(seq, &mut mins);
        for &(h, p) in mins.iter() {
            if counts.get(&h).map_or(false, |&c| c as usize > max_hits) { continue; }
            let bucket = (h >> SHIFT) as usize;
            if bucket >= write_offsets.len() { continue; }
            let idx = write_offsets[bucket] as usize;
            if idx >= final_seeds.len() { continue; }
            
            let hash_rem = (h >> (SHIFT - HASH_BITS)) & HASH_MASK;
            let rid_part = (rid as u64) & RID_MASK;
            let pos_part = (p as u64) & POS_MASK;
            final_seeds[idx] = (hash_rem << (RID_BITS + POS_BITS)) | (rid_part << POS_BITS) | pos_part;
            write_offsets[bucket] += 1;
        }
    }

    for bucket in 0..offsets.len() {
        let start = offsets[bucket] as usize;
        let end = if bucket + 1 < offsets.len() { offsets[bucket + 1] as usize } else { final_seeds.len() };
        if start < end { final_seeds[start..end].sort_unstable_by_key(|s| s >> (RID_BITS + POS_BITS)); }
    }

    eprintln!("  Index built: {} seeds", final_seeds.len());

    let mut seed_counts: Vec<u64> = counts.iter().map(|(&h, &c)| ((h >> 32) << 32) | (c as u64)).collect();
    seed_counts.sort_unstable_by_key(|v| v >> 32);

    Ok(Index { offsets, seeds: final_seeds, genome_size: total_bases, max_hits: max_hits as u32,
        total_seed_occurrences, seed_counts, ref_seqs, ref_names })
}

#[derive(Clone, Copy, Default, Debug)]
pub struct Hit { pub id_strand: u32, pub diag: i32, pub read_pos: u32, pub ref_pos: u32 }

pub struct State {
    pub mins: Vec<(u64, u32)>,
    pub candidates: Vec<Hit>,
    pub diag_counts: Vec<(i32, u16)>,
    pub group_buf: Vec<(u32, u32, u32)>, // (id_strand, start, len)
    pub cigar_buf: Vec<u32>,
}

impl State {
    pub fn new() -> Self {
        Self { 
            mins: Vec::with_capacity(128), 
            candidates: Vec::with_capacity(2048), 
            diag_counts: Vec::with_capacity(64),
            group_buf: Vec::with_capacity(32),
            cigar_buf: Vec::with_capacity(32),
        }
    }
}

#[inline(always)]
fn base_to_index(b: u8) -> Option<usize> {
    let v = BASE_LUT[b as usize];
    if v >= 0 { Some(v as usize) } else { None }
}

#[inline(always)]
pub fn get_syncmers(seq: &[u8], out: &mut Vec<(u64, u32)>) {
    out.clear();
    if seq.len() < HASH_WINDOW || HASH_WINDOW == 0 || SYNC_S > HASH_WINDOW { return; }

    let mut s_queue = [(0u32, 0u64); 32];
    let mut s_start = 0usize;
    let mut s_end = 0usize;
    let mut h_k = 0u64;
    let mut h_s = 0u64;
    let mut base_buf_k = [0usize; HASH_WINDOW];
    let mut base_buf_s = [0usize; SYNC_S];
    let mut ambig_k = 0u32;
    let mut ambig_s = 0u32;
    let mut prev_kpos = u32::MAX;

    for i in 0..seq.len() {
        let (base_idx, is_ambig) = match base_to_index(seq[i]) {
            Some(idx) => (idx, false),
            None => (0, true),
        };

        let k_slot = i % HASH_WINDOW;
        let prev_k = base_buf_k[k_slot];
        base_buf_k[k_slot] = base_idx;
        if is_ambig { ambig_k |= 1 << k_slot; } else { ambig_k &= !(1 << k_slot); }
        h_k = if i + 1 <= HASH_WINDOW {
            h_k.rotate_left(ROT) ^ BASES[base_idx]
        } else {
            h_k.rotate_left(ROT) ^ BASES[base_idx] ^ REMOVE[prev_k]
        };

        let s_slot = i % SYNC_S;
        let prev_s = base_buf_s[s_slot];
        base_buf_s[s_slot] = base_idx;
        if is_ambig { ambig_s |= 1 << s_slot; } else { ambig_s &= !(1 << s_slot); }
        h_s = if i + 1 <= SYNC_S {
            h_s.rotate_left(ROT) ^ BASES[base_idx]
        } else {
            h_s.rotate_left(ROT) ^ BASES[base_idx] ^ REMOVE_S[prev_s]
        };

        if i + 1 >= SYNC_S && ambig_s == 0 {
            let s_pos = (i + 1 - SYNC_S) as u32;
            let s_hash = h_s.wrapping_mul(0x517cc1b727220a95);
            while s_end > s_start && s_queue[(s_end - 1) & 31].1 >= s_hash { s_end -= 1; }
            s_queue[s_end & 31] = (s_pos, s_hash);
            s_end += 1;
        }

        if i + 1 >= HASH_WINDOW && ambig_k == 0 {
            let k_pos = (i + 1 - HASH_WINDOW) as u32;
            let max_s_pos = k_pos + (HASH_WINDOW - SYNC_S) as u32;
            while s_start < s_end && s_queue[s_start & 31].0 < k_pos { s_start += 1; }

            if s_start < s_end {
                let min_s_pos = s_queue[s_start & 31].0;
                if (min_s_pos == k_pos || min_s_pos == max_s_pos) && k_pos != prev_kpos {
                    out.push((h_k.wrapping_mul(0x517cc1b727220a95), k_pos));
                    prev_kpos = k_pos;
                }
            }
        }
    }
}

#[inline(always)]
fn collect_candidates<I: IndexLike>(idx: &I, mins: &[(u64, u32)], is_rev: bool, max_hits: usize, out: &mut Vec<Hit>) -> bool {
    let offsets = idx.offsets();
    let seeds = idx.seeds();
    let mut capped = false;
    
    for &(h, r_pos) in mins {
        if out.len() >= CONFIG.max_candidates { capped = true; break; }
        
        let bucket = (h >> SHIFT) as usize;
        if bucket >= offsets.len() { continue; }
        let start = offsets[bucket] as usize;
        let end = if bucket + 1 < offsets.len() { offsets[bucket + 1] as usize } else { seeds.len() };
        if start >= end { continue; }

        let target = (h >> (SHIFT - HASH_BITS)) & HASH_MASK;
        let search_key = target << (RID_BITS + POS_BITS);
        let bucket_seeds = &seeds[start..end];
        let pos = bucket_seeds.partition_point(|&s| s < search_key);
        
        let mut count = 0;
        let mut idx_in_bucket = pos;
        while idx_in_bucket < bucket_seeds.len() && count < max_hits {
            let seed = bucket_seeds[idx_in_bucket];
            if (seed >> (RID_BITS + POS_BITS)) != target { break; }
            let ref_id = ((seed >> POS_BITS) & RID_MASK) as u32;
            let ref_pos = (seed & POS_MASK) as u32;
            out.push(Hit {
                id_strand: (ref_id << 1) | (is_rev as u32),
                diag: (ref_pos as i32).wrapping_sub(r_pos as i32),
                read_pos: r_pos, ref_pos,
            });
            count += 1;
            idx_in_bucket += 1;
            if out.len() >= CONFIG.max_candidates { capped = true; break; }
        }
        if count >= max_hits { capped = true; }
    }
    capped
}

#[inline(always)]
fn push_cigar(cigar: &mut Vec<u32>, len: u32, op: u32) {
    if len == 0 { return; }
    if let Some(last) = cigar.last_mut() {
        if (*last & 0xF) == op { *last += len << 4; return; }
    }
    cigar.push((len << 4) | op);
}

#[inline(always)]
pub fn cigar_ref_span(cigar: &[u32]) -> i32 {
    cigar.iter().map(|&c| match c & 0xF { 0 | 2 | 3 | 7 | 8 => (c >> 4) as i32, _ => 0 }).sum()
}

#[inline(always)]
fn best_diag_for_group(hits: &[Hit], counts: &mut Vec<(i32, u16)>) -> Option<(i32, usize)> {
    if hits.is_empty() { return None; }
    counts.clear();
    
    for hit in hits {
        let d = hit.diag;
        if let Some(entry) = counts.iter_mut().find(|(diag, _)| *diag == d) {
            entry.1 += 1;
        } else {
            counts.push((d, 1));
        }
    }
    
    let (best_diag, best_count) = counts.iter().max_by_key(|&&(_, c)| c).copied()?;
    Some((best_diag, best_count as usize))
}

#[inline(always)]
fn anchor_for_diag(group: &[Hit], diag: i32) -> Option<(i32, u32)> {
    let mut first = None;
    let mut last = None;
    for hit in group { if hit.diag == diag { if first.is_none() { first = Some(hit.read_pos); } last = Some(hit.read_pos); } }
    Some((diag, (first? + last?) / 2))
}

const CASE_MASK: u64 = 0x2020202020202020;

#[inline(always)]
fn eq_base(a: u8, b: u8) -> bool {
    let a = if a == b'N' || a == b'n' { b'A' } else { a };
    let b = if b == b'N' || b == b'n' { b'A' } else { b };
    (a | 0x20) == (b | 0x20)
}

// SIMD-accelerated sequence comparison for x86_64
#[cfg(all(target_arch = "x86_64", not(target_family = "wasm")))]
#[target_feature(enable = "avx2")]
#[inline]
unsafe fn simd_cmp_32(a: &[u8], b: &[u8]) -> u32 {
    unsafe {
        let case_mask = _mm256_set1_epi8(0x20);
        let va = _mm256_or_si256(_mm256_loadu_si256(a.as_ptr() as *const __m256i), case_mask);
        let vb = _mm256_or_si256(_mm256_loadu_si256(b.as_ptr() as *const __m256i), case_mask);
        let cmp = _mm256_cmpeq_epi8(va, vb);
        _mm256_movemask_epi8(cmp) as u32
    }
}

#[cfg(all(target_arch = "x86_64", not(target_family = "wasm")))]
#[target_feature(enable = "sse2")]
#[inline]
unsafe fn simd_cmp_16(a: &[u8], b: &[u8]) -> u16 {
    unsafe {
        let case_mask = _mm_set1_epi8(0x20);
        let va = _mm_or_si128(_mm_loadu_si128(a.as_ptr() as *const __m128i), case_mask);
        let vb = _mm_or_si128(_mm_loadu_si128(b.as_ptr() as *const __m128i), case_mask);
        let cmp = _mm_cmpeq_epi8(va, vb);
        _mm_movemask_epi8(cmp) as u16
    }
}

// Scalar fallback for WASM and other platforms
#[inline(always)]
fn scalar_cmp_8(a: &[u8], b: &[u8]) -> bool {
    unsafe {
        let r_val = std::ptr::read_unaligned(a.as_ptr() as *const u64) | CASE_MASK;
        let g_val = std::ptr::read_unaligned(b.as_ptr() as *const u64) | CASE_MASK;
        r_val == g_val
    }
}

#[inline(always)]
fn extend_left(read: &[u8], rseq: &[u8], mut rp: usize, mut gp: usize, cfg: &Config) -> (i32, usize, usize, Vec<u32>, usize) {
    let (mut score, mut max_score) = (0i32, 0i32);
    let mut cigar = Vec::with_capacity(8);
    let mut match_run = 0u32;
    let (mut best_rp, mut best_gp, mut best_cigar_len, mut best_match_run) = (rp, gp, 0, 0);
    let (mut nm, mut best_nm) = (0usize, 0usize);

    // Use SIMD on x86_64 (non-WASM)
    #[cfg(all(target_arch = "x86_64", not(target_family = "wasm")))]
    let use_avx2 = has_avx2();
    #[cfg(not(all(target_arch = "x86_64", not(target_family = "wasm"))))]
    let _use_avx2 = false;

    while rp > 0 && gp > 0 {
        // Try 32-byte SIMD comparison (AVX2)
        #[cfg(all(target_arch = "x86_64", not(target_family = "wasm")))]
        if use_avx2 && rp >= 32 && gp >= 32 {
            let mask = unsafe { simd_cmp_32(&read[rp - 32..rp], &rseq[gp - 32..gp]) };
            if mask == 0xFFFFFFFF {
                score += cfg.match_score * 32; match_run += 32; rp -= 32; gp -= 32;
                if score > max_score { max_score = score; best_rp = rp; best_gp = gp; best_cigar_len = cigar.len(); best_match_run = match_run; best_nm = nm; }
                continue;
            }
        }

        // Try 16-byte SIMD comparison (SSE2)
        #[cfg(all(target_arch = "x86_64", not(target_family = "wasm")))]
        if rp >= 16 && gp >= 16 {
            let mask = unsafe { simd_cmp_16(&read[rp - 16..rp], &rseq[gp - 16..gp]) };
            if mask == 0xFFFF {
                score += cfg.match_score * 16; match_run += 16; rp -= 16; gp -= 16;
                if score > max_score { max_score = score; best_rp = rp; best_gp = gp; best_cigar_len = cigar.len(); best_match_run = match_run; best_nm = nm; }
                continue;
            }
        }

        // 8-byte scalar comparison (portable, works on WASM)
        if rp >= 8 && gp >= 8 && scalar_cmp_8(&read[rp - 8..rp], &rseq[gp - 8..gp]) {
            score += cfg.match_score * 8; match_run += 8; rp -= 8; gp -= 8;
            if score > max_score { max_score = score; best_rp = rp; best_gp = gp; best_cigar_len = cigar.len(); best_match_run = match_run; best_nm = nm; }
            continue;
        }

        let (rb, gb) = (read[rp - 1], rseq[gp - 1]);
        if eq_base(rb, gb) { score += cfg.match_score; match_run += 1; rp -= 1; gp -= 1; }
        else {
            if match_run > 0 { cigar.push((match_run << 4) | 0); match_run = 0; }
            let (mm, del, ins) = (score + cfg.mismatch_pen, score + cfg.gap_open, score + cfg.gap_open);
            if mm >= del && mm >= ins { score = mm; cigar.push((1 << 4) | 0); nm += 1; rp -= 1; gp -= 1; }
            else if del >= ins && gp > 0 {
                score = del; let mut glen = 1usize; nm += 1;
                while glen < cfg.maxindel && gp > glen && score + cfg.gap_ext > max_score - cfg.x_drop {
                    if rp > 0 && eq_base(read[rp - 1], rseq[gp - glen - 1]) { break; }
                    glen += 1; nm += 1; score += cfg.gap_ext;
                }
                cigar.push(((glen as u32) << 4) | 2); gp = gp.saturating_sub(glen);
            } else if rp > 0 { score = ins; cigar.push((1 << 4) | 1); nm += 1; rp -= 1; }
            else { break; }
        }
        if score > max_score { max_score = score; best_rp = rp; best_gp = gp; best_cigar_len = cigar.len(); best_match_run = match_run; best_nm = nm; }
        if score < max_score - cfg.x_drop { break; }
    }
    cigar.truncate(best_cigar_len);
    if best_match_run > 0 { cigar.push((best_match_run << 4) | 0); }
    cigar.reverse();
    (max_score, best_rp, best_gp, cigar, best_nm)
}

#[inline(always)]
fn extend_right(read: &[u8], rseq: &[u8], mut rp: usize, mut gp: usize, cfg: &Config) -> (i32, usize, usize, Vec<u32>, usize) {
    let (rlen, glen) = (read.len(), rseq.len());
    let (mut score, mut max_score) = (0i32, 0i32);
    let mut cigar = Vec::with_capacity(8);
    let mut match_run = 0u32;
    let (mut best_rp, mut best_gp, mut best_cigar_len, mut best_match_run) = (rp, gp, 0, 0);
    let (mut nm, mut best_nm) = (0usize, 0usize);

    // Use SIMD on x86_64 (non-WASM)
    #[cfg(all(target_arch = "x86_64", not(target_family = "wasm")))]
    let use_avx2 = has_avx2();
    #[cfg(not(all(target_arch = "x86_64", not(target_family = "wasm"))))]
    let _use_avx2 = false;

    while rp < rlen && gp < glen {
        // Try 32-byte SIMD comparison (AVX2)
        #[cfg(all(target_arch = "x86_64", not(target_family = "wasm")))]
        if use_avx2 && rlen - rp >= 32 && glen - gp >= 32 {
            let mask = unsafe { simd_cmp_32(&read[rp..rp + 32], &rseq[gp..gp + 32]) };
            if mask == 0xFFFFFFFF {
                rp += 32; gp += 32; score += cfg.match_score * 32; match_run += 32;
                if score > max_score { max_score = score; best_rp = rp; best_gp = gp; best_cigar_len = cigar.len(); best_match_run = match_run; best_nm = nm; }
                continue;
            }
        }

        // Try 16-byte SIMD comparison (SSE2)
        #[cfg(all(target_arch = "x86_64", not(target_family = "wasm")))]
        if rlen - rp >= 16 && glen - gp >= 16 {
            let mask = unsafe { simd_cmp_16(&read[rp..rp + 16], &rseq[gp..gp + 16]) };
            if mask == 0xFFFF {
                rp += 16; gp += 16; score += cfg.match_score * 16; match_run += 16;
                if score > max_score { max_score = score; best_rp = rp; best_gp = gp; best_cigar_len = cigar.len(); best_match_run = match_run; best_nm = nm; }
                continue;
            }
        }

        // 8-byte scalar comparison (portable, works on WASM)
        if rlen - rp >= 8 && glen - gp >= 8 && scalar_cmp_8(&read[rp..rp + 8], &rseq[gp..gp + 8]) {
            rp += 8; gp += 8; score += cfg.match_score * 8; match_run += 8;
            if score > max_score { max_score = score; best_rp = rp; best_gp = gp; best_cigar_len = cigar.len(); best_match_run = match_run; best_nm = nm; }
            continue;
        }

        let (rb, gb) = (read[rp], rseq[gp]);
        if eq_base(rb, gb) { score += cfg.match_score; match_run += 1; rp += 1; gp += 1; }
        else {
            if match_run > 0 { cigar.push((match_run << 4) | 0); match_run = 0; }
            let (mm, del, ins) = (score + cfg.mismatch_pen, score + cfg.gap_open, score + cfg.gap_open);
            if mm >= del && mm >= ins { score = mm; cigar.push((1 << 4) | 0); nm += 1; rp += 1; gp += 1; }
            else if del >= ins && gp < glen {
                score = del; let mut dlen = 1usize; nm += 1;
                while dlen < cfg.maxindel && gp + dlen < glen && score + cfg.gap_ext > max_score - cfg.x_drop {
                    if rp < rlen && eq_base(read[rp], rseq[gp + dlen]) { break; }
                    dlen += 1; nm += 1; score += cfg.gap_ext;
                }
                cigar.push(((dlen as u32) << 4) | 2); gp += dlen;
            } else if rp < rlen { score = ins; cigar.push((1 << 4) | 1); nm += 1; rp += 1; }
            else { break; }
        }
        if score > max_score { max_score = score; best_rp = rp; best_gp = gp; best_cigar_len = cigar.len(); best_match_run = match_run; best_nm = nm; }
        if score < max_score - cfg.x_drop { break; }
    }
    cigar.truncate(best_cigar_len);
    if best_match_run > 0 { cigar.push((best_match_run << 4) | 0); }
    (max_score, best_rp, best_gp, cigar, best_nm)
}

#[inline(always)]
fn greedy_extend(anchor_rp: usize, anchor_diag: i32, read: &[u8], rseq: &[u8], cfg: &Config) -> (Vec<u32>, i32, i32, usize, usize) {
    let a_gp = if anchor_diag >= 0 { (anchor_diag as usize).saturating_add(anchor_rp).min(rseq.len().saturating_sub(1)) }
               else { anchor_rp.saturating_sub((-anchor_diag) as usize) };
    let (l_score, l_rp, l_gp, l_cigar, l_nm) = extend_left(read, rseq, anchor_rp, a_gp, cfg);
    let (r_score, r_rp, _, r_cigar, r_nm) = extend_right(read, rseq, anchor_rp, a_gp, cfg);
    let mut cigar = Vec::with_capacity(l_cigar.len() + r_cigar.len() + 2);
    if l_rp > 0 { cigar.push(((l_rp as u32) << 4) | 4); }
    for &c in &l_cigar { push_cigar(&mut cigar, c >> 4, c & 0xF); }
    for &c in &r_cigar { push_cigar(&mut cigar, c >> 4, c & 0xF); }
    let trailing = read.len().saturating_sub(r_rp);
    if trailing > 0 { push_cigar(&mut cigar, trailing as u32, 4); }
    (cigar, l_gp as i32, l_score + r_score, r_rp.saturating_sub(l_rp), l_nm + r_nm)
}

pub fn align<I: IndexLike>(seq: &[u8], idx: &I, state: &mut State, rev: &mut Vec<u8>) -> Vec<AlignmentResult> {
    let State { mins, candidates, diag_counts, group_buf, cigar_buf: _ } = state;
    let cfg = &CONFIG;
    candidates.clear();
    let max_hits = idx.max_hits();
    let mut capped = false;

    get_syncmers(seq, mins);
    capped |= collect_candidates(idx, mins, false, max_hits, candidates);

    rev.clear();
    rev.extend(seq.iter().rev().map(|b| match b { b'A'|b'a' => b'T', b'C'|b'c' => b'G', b'G'|b'g' => b'C', b'T'|b't' => b'A', _ => b'N' }));
    get_syncmers(rev, mins);
    capped |= collect_candidates(idx, mins, true, max_hits, candidates);
    
    if candidates.is_empty() { return Vec::new(); }

    // Sort by id_strand to group hits
    candidates.sort_unstable_by_key(|h| h.id_strand);
    
    // Build group index without HashMap
    group_buf.clear();
    let mut prev_id = candidates[0].id_strand;
    let mut start = 0u32;
    for (i, hit) in candidates.iter().enumerate() {
        if hit.id_strand != prev_id {
            group_buf.push((prev_id, start, (i as u32) - start));
            prev_id = hit.id_strand;
            start = i as u32;
        }
    }
    group_buf.push((prev_id, start, (candidates.len() as u32) - start));

    // Find best diagonal for each group and sort by vote count
    let mut group_votes: Vec<(u32, i32, usize, u32, u32)> = Vec::with_capacity(group_buf.len());
    for &(id, start, len) in group_buf.iter() {
        let hits = &candidates[start as usize..(start + len) as usize];
        if let Some((diag, count)) = best_diag_for_group(hits, diag_counts) {
            group_votes.push((id, diag, count, start, len));
        }
    }
    if group_votes.is_empty() { return Vec::new(); }

    group_votes.sort_unstable_by(|a, b| b.2.cmp(&a.2));
    let (best_id, best_diag, best_count, best_start, best_len) = group_votes[0];
    if best_count < cfg.min_seeds { return Vec::new(); }
    let second_count = group_votes.get(1).map(|x| x.2).unwrap_or(0);
    let total_seeds = candidates.len();

    let attempt = |id: u32, diag: i32, start: u32, len: u32, rev_buf: &[u8]| -> Option<AlignmentResult> {
        let ref_id = (id >> 1) as i32;
        if (ref_id as usize) >= idx.ref_count() { return None; }
        let is_rev = (id & 1) == 1;
        let hits = &candidates[start as usize..(start + len) as usize];
        let (anchor_diag, anchor_rp) = anchor_for_diag(hits, diag)?;
        let target_seq = if is_rev { rev_buf } else { seq };
        let ref_seq = idx.ref_seq(ref_id as usize);
        let (cigar, pos, as_score, aligned_len, nm) = greedy_extend(anchor_rp as usize, anchor_diag, target_seq, ref_seq, cfg);
        if aligned_len == 0 || pos < 0 || (pos as usize) >= ref_seq.len() { return None; }
        Some(AlignmentResult { ref_id, pos, is_rev, mapq: 0, cigar, nm: nm as i32, md: String::new(), as_score, paired: false, proper_pair: false })
    };

    let mut results = Vec::with_capacity(2);
    if let Some(mut res) = attempt(best_id, best_diag, best_start, best_len, rev) {
        let diff = best_count.saturating_sub(second_count);
        let mapq = ((diff as f32 / total_seeds.max(1) as f32) * 60.0 + (diff.min(10) as f32) * 3.0).round().min(60.0) as u8;
        res.mapq = if capped { mapq.min(30) } else { mapq };
        results.push(res);
    }

    if group_votes.len() > 1 && second_count >= cfg.min_seeds {
        let (sec_id, sec_diag, _, sec_start, sec_len) = group_votes[1];
        if let Some(mut res) = attempt(sec_id, sec_diag, sec_start, sec_len, rev) { 
            res.mapq = 0; 
            results.push(res); 
        }
    }

    results.sort_by(|a, b| b.as_score.cmp(&a.as_score));
    results
}

pub fn write_cigar_string(cigar: &[u32], out: &mut Vec<u8>) {
    if cigar.is_empty() { out.push(b'*'); return; }
    for &c in cigar {
        let op = match c & 0xF { 0 => b'M', 1 => b'I', 2 => b'D', 3 => b'N', 4 => b'S', 5 => b'H', 6 => b'P', 7 => b'=', 8 => b'X', _ => b'?' };
        write!(out, "{}{}", c >> 4, op as char).unwrap();
    }
}

pub fn reverse_qual(qual: &[u8]) -> Vec<u8> { qual.iter().rev().cloned().collect() }

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|b| match b { b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A', b'a' => b't', b'c' => b'g', b'g' => b'c', b't' => b'a', _ => b'N' }).collect()
}

pub fn oriented_bases<'a>(seq: &'a [u8], qual: &'a [u8], res: &Option<AlignmentResult>) -> (Cow<'a, [u8]>, Cow<'a, [u8]>) {
    if let Some(r) = res { if r.is_rev { return (Cow::Owned(reverse_complement(seq)), Cow::Owned(reverse_qual(qual))); } }
    (Cow::Borrowed(seq), Cow::Borrowed(qual))
}
