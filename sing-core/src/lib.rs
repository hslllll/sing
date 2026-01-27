

use anyhow::{bail, Context, Result};
use bytemuck::try_cast_slice;
use memmap2::MmapOptions;
use std::borrow::Cow;
use std::fs::File;
use std::io::{BufReader, Read, Write};
use std::ops::Range;
use std::path::Path;

const WINDOW: usize = 16;
const SYNC_S: usize = 4;

const BASES: [u64; 4] = [
    (std::f64::consts::E - 2f64).to_bits(),
    (std::f64::consts::PI - 3.0).to_bits(),
    (std::f64::consts::SQRT_2 - 1.0).to_bits(),
    (1.7320508075688772f64 - 1.0).to_bits(),
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
    rot(BASES[0], (WINDOW as u32 * ROT) % 64),
    rot(BASES[1], (WINDOW as u32 * ROT) % 64),
    rot(BASES[2], (WINDOW as u32 * ROT) % 64),
    rot(BASES[3], (WINDOW as u32 * ROT) % 64),
];

const REMOVE_S: [u64; 4] = [
    rot(BASES[0], (SYNC_S as u32 * ROT) % 64),
    rot(BASES[1], (SYNC_S as u32 * ROT) % 64),
    rot(BASES[2], (SYNC_S as u32 * ROT) % 64),
    rot(BASES[3], (SYNC_S as u32 * ROT) % 64),
];

pub const RADIX: usize = 24;
pub const SHIFT: usize = 32 - RADIX;

#[inline(always)]
fn scoring_params() -> (i32, i32) {
    let match_score = (WINDOW as i32 / 10).max(1);
    let mismatch_pen = -((SYNC_S as i32 / 10).max(1));
    (match_score, mismatch_pen)
}

#[derive(Debug, Clone)]
pub struct Index {
    pub offsets: Vec<u32>,
    pub seeds: Vec<u64>,
    pub freq_filter: u32,
    pub genome_size: u64,
    pub ref_seqs: Vec<Vec<u8>>,
    pub ref_names: Vec<String>,
}

pub trait IndexLike {
    fn offsets(&self) -> &[u32];
    fn seeds(&self) -> &[u64];
    fn freq_filter(&self) -> u32;
    fn genome_size(&self) -> u64;
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

        fn read_u32(data: &[u8], pos: &mut usize, label: &str) -> Result<u32> {
            if *pos + 4 > data.len() {
                bail!("{}: unexpected EOF", label);
            }
            let val = u32::from_le_bytes(data[*pos..*pos + 4].try_into().unwrap());
            *pos += 4;
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

        let freq_filter = read_u32(data, &mut pos, "read freq_filter")?;
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

        Ok(Index { offsets, seeds, freq_filter, genome_size, ref_seqs, ref_names })
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

        r.read_exact(&mut buf8).context("read seeds len")?;
        let len = u64::from_le_bytes(buf8) as usize;
        let mut seeds = Vec::with_capacity(len);
        for _ in 0..len {
            r.read_exact(&mut buf8).context("read seed")?;
            seeds.push(u64::from_le_bytes(buf8));
        }

        r.read_exact(&mut buf4).context("read freq_filter")?;
        let freq_filter = u32::from_le_bytes(buf4);

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

        Ok(Index { offsets, seeds, freq_filter, genome_size, ref_seqs, ref_names })
    }

    pub fn to_writer<W: Write>(&self, mut w: W) -> Result<()> {
        w.write_all(&(self.offsets.len() as u64).to_le_bytes())?;
        for v in &self.offsets {
            w.write_all(&v.to_le_bytes())?;
        }

        w.write_all(&(self.seeds.len() as u64).to_le_bytes())?;
        for s in &self.seeds {
            w.write_all(&s.to_le_bytes())?;
        }

        w.write_all(&self.freq_filter.to_le_bytes())?;

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
        Ok(())
    }
}

impl IndexLike for Index {
    fn offsets(&self) -> &[u32] {
        &self.offsets
    }

    fn seeds(&self) -> &[u64] {
        &self.seeds
    }

    fn freq_filter(&self) -> u32 {
        self.freq_filter
    }

    fn genome_size(&self) -> u64 {
        self.genome_size
    }

    fn ref_count(&self) -> usize {
        self.ref_seqs.len()
    }

    fn ref_seq(&self, id: usize) -> &[u8] {
        &self.ref_seqs[id]
    }

    fn ref_name(&self, id: usize) -> &str {
        &self.ref_names[id]
    }
}

enum IndexBuffer {
    Mmap(memmap2::Mmap),
}

impl IndexBuffer {
    #[inline]
    fn as_slice(&self) -> &[u8] {
        match self {
            Self::Mmap(mmap) => mmap.as_ref(),
        }
    }
}

pub struct MemoryIndex {
    buf: IndexBuffer,
    offsets_ptr: *const u32,
    offsets_len: usize,
    seeds_ptr: *const u64,
    seeds_len: usize,
    freq_filter: u32,
    genome_size: u64,
    ref_seq_ranges: Vec<Range<usize>>,
    ref_name_ranges: Vec<Range<usize>>,
}

unsafe impl Send for MemoryIndex {}
unsafe impl Sync for MemoryIndex {}

impl MemoryIndex {
    pub fn from_path<P: AsRef<Path>>(path: P) -> Result<Self> {
        let file = File::open(path)?;
        let mmap = unsafe {
            MmapOptions::new()
                .populate()
                .map(&file)
                .context("mmap index file")?
        };
        #[cfg(all(target_os = "linux", not(target_arch = "wasm32")))]
        unsafe {
            let ptr = mmap.as_ptr() as *mut libc::c_void;
            let len = mmap.len();
            let advice = libc::MADV_RANDOM;
            let _ = libc::madvise(ptr, len, advice);
        }
        Self::from_data(IndexBuffer::Mmap(mmap))
    }

    fn from_data(buf: IndexBuffer) -> Result<Self> {
        let data = buf.as_slice();
        let mut pos = 0usize;

        fn read_u64(data: &[u8], pos: &mut usize, label: &str) -> Result<u64> {
            if *pos + 8 > data.len() {
                bail!("{}: unexpected EOF", label);
            }
            let val = u64::from_le_bytes(data[*pos..*pos + 8].try_into().unwrap());
            *pos += 8;
            Ok(val)
        }

        fn read_u32(data: &[u8], pos: &mut usize, label: &str) -> Result<u32> {
            if *pos + 4 > data.len() {
                bail!("{}: unexpected EOF", label);
            }
            let val = u32::from_le_bytes(data[*pos..*pos + 4].try_into().unwrap());
            *pos += 4;
            Ok(val)
        }

        let offsets_len = read_u64(data, &mut pos, "read offsets len")? as usize;
        let offsets_bytes = offsets_len
            .checked_mul(4)
            .and_then(|n| pos.checked_add(n))
            .filter(|&end| end <= data.len())
            .ok_or_else(|| anyhow::anyhow!("read offsets: unexpected EOF"))?;
        let offsets_range = pos..offsets_bytes;
        let offsets_slice = try_cast_slice::<u8, u32>(&data[offsets_range])
            .map_err(|_| anyhow::anyhow!("read offsets: unaligned data"))?;
        let offsets_ptr = offsets_slice.as_ptr();
        let offsets_len = offsets_slice.len();
        pos = offsets_bytes;

        let seeds_len = read_u64(data, &mut pos, "read seeds len")? as usize;
        let seeds_bytes = seeds_len
            .checked_mul(8)
            .and_then(|n| pos.checked_add(n))
            .filter(|&end| end <= data.len())
            .ok_or_else(|| anyhow::anyhow!("read seeds: unexpected EOF"))?;
        let seeds_range = pos..seeds_bytes;
        let seeds_slice = try_cast_slice::<u8, u64>(&data[seeds_range])
            .map_err(|_| anyhow::anyhow!("read seeds: unaligned data"))?;
        let seeds_ptr = seeds_slice.as_ptr();
        let seeds_len = seeds_slice.len();
        pos = seeds_bytes;

        let freq_filter = read_u32(data, &mut pos, "read freq_filter")?;
        let genome_size = read_u64(data, &mut pos, "read genome_size")?;

        let ref_count = read_u64(data, &mut pos, "read ref seq count")? as usize;
        let mut ref_seq_ranges = Vec::with_capacity(ref_count);
        for _ in 0..ref_count {
            let slen = read_u64(data, &mut pos, "read seq len")? as usize;
            let end = pos
                .checked_add(slen)
                .filter(|&e| e <= data.len())
                .ok_or_else(|| anyhow::anyhow!("read seq bytes: unexpected EOF"))?;
            ref_seq_ranges.push(pos..end);
            pos = end;
        }

        let name_count = read_u64(data, &mut pos, "read name count")? as usize;
        let mut ref_name_ranges = Vec::with_capacity(name_count);
        for _ in 0..name_count {
            let slen = read_u64(data, &mut pos, "read name len")? as usize;
            let end = pos
                .checked_add(slen)
                .filter(|&e| e <= data.len())
                .ok_or_else(|| anyhow::anyhow!("read name bytes: unexpected EOF"))?;
            let name_buf = &data[pos..end];
            std::str::from_utf8(name_buf).context("utf8 ref name")?;
            ref_name_ranges.push(pos..end);
            pos = end;
        }

        Ok(Self {
            buf,
            offsets_ptr,
            offsets_len,
            seeds_ptr,
            seeds_len,
            freq_filter,
            genome_size,
            ref_seq_ranges,
            ref_name_ranges,
        })
    }
}

impl IndexLike for MemoryIndex {
    fn offsets(&self) -> &[u32] {
        unsafe { std::slice::from_raw_parts(self.offsets_ptr, self.offsets_len) }
    }

    fn seeds(&self) -> &[u64] {
        unsafe { std::slice::from_raw_parts(self.seeds_ptr, self.seeds_len) }
    }

    fn freq_filter(&self) -> u32 {
        self.freq_filter
    }

    fn genome_size(&self) -> u64 {
        self.genome_size
    }

    fn ref_count(&self) -> usize {
        self.ref_seq_ranges.len()
    }

    fn ref_seq(&self, id: usize) -> &[u8] {
        let data = self.buf.as_slice();
        &data[self.ref_seq_ranges[id].clone()]
    }

    fn ref_name(&self, id: usize) -> &str {
        let data = self.buf.as_slice();
        let bytes = &data[self.ref_name_ranges[id].clone()];
        unsafe { std::str::from_utf8_unchecked(bytes) }
    }
}

impl<T: IndexLike> IndexLike for std::sync::Arc<T> {
    fn offsets(&self) -> &[u32] {
        self.as_ref().offsets()
    }

    fn seeds(&self) -> &[u64] {
        self.as_ref().seeds()
    }

    fn freq_filter(&self) -> u32 {
        self.as_ref().freq_filter()
    }

    fn genome_size(&self) -> u64 {
        self.as_ref().genome_size()
    }

    fn ref_count(&self) -> usize {
        self.as_ref().ref_count()
    }

    fn ref_seq(&self, id: usize) -> &[u8] {
        self.as_ref().ref_seq(id)
    }

    fn ref_name(&self, id: usize) -> &str {
        self.as_ref().ref_name(id)
    }
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
}

pub fn build_index_from_sequences<I, S>(records: I) -> Result<Index>
where
    I: IntoIterator<Item = (S, Vec<u8>)>,
    S: Into<String>,
{
    let mut ref_seqs = Vec::new();
    let mut ref_names = Vec::new();

    for (_rid, (name, seq)) in records.into_iter().enumerate() {
        ref_names.push(name.into());
        ref_seqs.push(seq);
    }

    let total_bases: u64 = ref_seqs.iter().map(|s| s.len() as u64).sum();
    let freq_filter = (total_bases as f64).powf(1.0 / 3.0) as usize;
    let freq_filter = std::cmp::max(1, freq_filter);
    eprintln!("Dynamic freq_filter set to: {} (Genome size: {} bp)", freq_filter, total_bases);

    
    
    let mut global_counts = vec![0u32; 1 << RADIX];
    let mut mins = Vec::new();
    for seq in &ref_seqs {
        get_syncmers(seq, &mut mins);
        for (h, _) in mins.iter() {
            let bucket = (*h as u64 >> SHIFT) as usize;
            global_counts[bucket] += 1;
        }
    }

    
    let mut kept_counts = vec![0u32; 1 << RADIX];
    let mut kept_buckets = 0u64;
    let mut total_seeds_kept = 0u64;
    for (i, &count) in global_counts.iter().enumerate() {
        if count as usize <= freq_filter {
            kept_counts[i] = count;
            if count > 0 {
                kept_buckets += 1;
                total_seeds_kept += count as u64;
            }
        }
    }

    
    let mut offsets = Vec::with_capacity(kept_counts.len());
    let mut current_offset = 0u32;
    for count in &kept_counts {
        offsets.push(current_offset);
        current_offset += *count;
    }

    
    
    let mut seeds = vec![0u64; total_seeds_kept as usize];
    let mut write_pos = offsets.clone();
    for (rid, seq) in ref_seqs.iter().enumerate() {
        get_syncmers(seq, &mut mins);
        for (h, p) in mins.iter() {
            let hash64 = *h as u64;
            let bucket = (hash64 >> SHIFT) as usize;
            if kept_counts[bucket] == 0 {
                continue;
            }
            let hash16 = (hash64 & 0xFFFF) as u64;
            let packed = (hash16 << 48) | ((rid as u64) << 32) | (*p as u64);
            let pos = write_pos[bucket] as usize;
            seeds[pos] = packed;
            write_pos[bucket] += 1;
        }
    }

    eprintln!("Indexing statistics:");
    eprintln!("  Total seeds (all buckets): {}", global_counts.iter().map(|&c| c as u64).sum::<u64>());
    eprintln!("  Buckets kept: {} / {}", kept_buckets, kept_counts.len());
    eprintln!("  Total seeds kept: {}", seeds.len());

    Ok(Index { offsets, seeds, freq_filter: freq_filter as u32, genome_size: total_bases, ref_seqs, ref_names })
}

#[derive(Clone, Copy, Default, Debug)]
pub struct Hit {
    pub id_strand: u32,
    pub diag: i32,
    pub read_pos: u32,
    pub ref_pos: u32,
}

#[derive(Clone, Copy, Default)]
pub struct ChainRes {
    pub score: i32,
    pub start: usize,
    pub end: usize,
}

pub struct State {
    pub mins: Vec<(u32, u32)>,
    pub candidates: Vec<Hit>,
}

impl State {
    pub fn new() -> Self {
        Self {
            mins: Vec::with_capacity(100),
            candidates: Vec::with_capacity(1000),
        }
    }
}

#[inline(always)]
fn base_to_index(b: u8) -> Option<usize> {
    let v = BASE_LUT[b as usize];
    if v >= 0 { Some(v as usize) } else { None }
}

#[inline(always)]
pub fn get_syncmers(seq: &[u8], out: &mut Vec<(u32, u32)>) {
    out.clear();
    if seq.len() < WINDOW || WINDOW < SYNC_S {
        return;
    }

    const Q_SIZE: usize = 32;
    const Q_MASK: usize = Q_SIZE - 1;
    let mut q_pos = [0u32; Q_SIZE];
    let mut q_hash = [u64::MAX; Q_SIZE];
    let mut head = 0;
    let mut tail = 0;

    let mut h_k = 0u64;
    let mut h_s = 0u64;
    let mut base_buf_k = [0usize; WINDOW];
    let mut base_buf_s = [0usize; SYNC_S];
    let mut ambig_buf_k = [0u8; WINDOW];
    let mut ambig_buf_s = [0u8; SYNC_S];
    let mut ambig_k = 0i32;
    let mut ambig_s = 0i32;

    for i in 0..seq.len() {
        let (base_idx, ambig) = match base_to_index(seq[i]) {
            Some(idx) => (idx, 0u8),
            None => (0usize, 1u8),
        };

        let k_slot = i % WINDOW;
        let prev_idx_k = base_buf_k[k_slot];
        let prev_ambig_k = ambig_buf_k[k_slot];
        base_buf_k[k_slot] = base_idx;
        ambig_buf_k[k_slot] = ambig;
        ambig_k += ambig as i32 - prev_ambig_k as i32;

        let s_slot = i % SYNC_S;
        let prev_idx_s = base_buf_s[s_slot];
        let prev_ambig_s = ambig_buf_s[s_slot];
        base_buf_s[s_slot] = base_idx;
        ambig_buf_s[s_slot] = ambig;
        ambig_s += ambig as i32 - prev_ambig_s as i32;

        if i + 1 <= WINDOW {
            h_k = h_k.rotate_left(ROT) ^ BASES[base_idx];
        } else {
            h_k = h_k.rotate_left(ROT) ^ BASES[base_idx] ^ REMOVE[prev_idx_k];
        }

        if i + 1 <= SYNC_S {
            h_s = h_s.rotate_left(ROT) ^ BASES[base_idx];
        } else {
            h_s = h_s.rotate_left(ROT) ^ BASES[base_idx] ^ REMOVE_S[prev_idx_s];
        }

        if i + 1 < SYNC_S {
            continue;
        }

        let s_pos = i + 1 - SYNC_S;
        if ambig_s == 0 {
            let s_hash = h_s.wrapping_mul(0x517cc1b727220a95);
            let s_pos_u32 = s_pos as u32;

            while tail > head {
                if q_hash[(tail - 1) & Q_MASK] >= s_hash {
                    tail -= 1;
                } else {
                    break;
                }
            }

            q_hash[tail & Q_MASK] = s_hash;
            q_pos[tail & Q_MASK] = s_pos_u32;
            tail += 1;
        }

        if i + 1 >= WINDOW && ambig_k == 0 {
            let k_pos = i + 1 - WINDOW;
            let min_pos = k_pos as u32;
            let max_pos = (k_pos + WINDOW - SYNC_S) as u32;

            while head < tail && q_pos[head & Q_MASK] < min_pos {
                head += 1;
            }

            if head < tail {
                let m_pos = q_pos[head & Q_MASK];
                if m_pos <= max_pos {
                    let k_hash = h_k.wrapping_mul(0x517cc1b727220a95);
                    if out.last().map(|&(_, p)| p) != Some(k_pos as u32) {
                         out.push((k_hash as u32, k_pos as u32));
                    }
                }
            }
        }
    }
}

#[inline(always)]
fn collect_candidates<I: IndexLike>(idx: &I, mins: &[(u32, u32)], is_rev: bool, out: &mut Vec<Hit>) {
    let freq_filter = idx.freq_filter() as usize;
    let offsets = idx.offsets();
    let seeds = idx.seeds();
    let offsets_len = offsets.len();
    let mut last_bucket = usize::MAX;
    let mut start = 0usize;
    let mut end = 0usize;
    for &(h, r_pos) in mins {
        let bucket = (h >> SHIFT) as usize;
        if bucket != last_bucket {
            start = unsafe { *offsets.get_unchecked(bucket) } as usize;
            end = if bucket + 1 < offsets_len {
                (unsafe { *offsets.get_unchecked(bucket + 1) }) as usize
            } else {
                seeds.len()
            };
            last_bucket = bucket;
        }

        let count = end - start;

        if count > freq_filter {
            continue;
        }

        let target_hash = (h & 0xFFFF) as u64;
        let mut i = start;
        while i < end {
            let seed = unsafe { *seeds.get_unchecked(i) };
            if (seed >> 48) == target_hash {
                let rid = ((seed >> 32) & 0xFFFF) as u32;
                let pos = seed as u32;

                let id_strand = (rid << 1) | (is_rev as u32);
                let diag = (pos as i32) - (r_pos as i32);

                out.push(Hit { id_strand, diag, read_pos: r_pos, ref_pos: pos });
            }
            i += 1;
        }
    }
}

#[inline(always)]
fn find_top_chains(candidates: &mut Vec<Hit>) -> (ChainRes, ChainRes) {
    if candidates.is_empty() {
        return (ChainRes::default(), ChainRes::default());
    }

    let bin_width = (WINDOW as i32).max(1);
    candidates.sort_unstable_by(|a, b| {
        let a_key = (a.id_strand, a.diag.div_euclid(bin_width));
        let b_key = (b.id_strand, b.diag.div_euclid(bin_width));
        a_key.cmp(&b_key)
    });

    let mut best = ChainRes::default();
    let mut second = ChainRes::default();

    let mut i = 0usize;
    while i < candidates.len() {
        let key = (candidates[i].id_strand, candidates[i].diag.div_euclid(bin_width));
        let start = i;
        let mut end = i + 1;
        while end < candidates.len() {
            let next_key = (candidates[end].id_strand, candidates[end].diag.div_euclid(bin_width));
            if next_key != key {
                break;
            }
            end += 1;
        }
        let score = (end - start) as i32;
        if score > best.score {
            second = best;
            best = ChainRes { score, start, end };
        } else if score > second.score {
            second = ChainRes { score, start, end };
        }
        i = end;
    }

    (best, second)
}


#[inline(always)]
fn push_cigar(cigar: &mut Vec<u32>, len: u32, op: u32) {
    if len == 0 {
        return;
    }
    if let Some(last) = cigar.last_mut() {
        if (*last & 0xF) == op {
            *last += len << 4;
            return;
        }
    }
    cigar.push((len << 4) | op);
}

#[inline(always)]
pub fn cigar_ref_span(cigar: &[u32]) -> i32 {
    let mut span = 0;
    for &c in cigar {
        let len = (c >> 4) as i32;
        let op = c & 0xF;
        match op {
            0 | 2 | 3 | 7 | 8 => span += len,
            _ => {}
        }
    }
    span
}


fn gen_cigar(
    hits: &[Hit],
    read_seq: &[u8],
    ref_seq: &[u8],
 ) -> (Vec<u32>, i32, usize) {
    if hits.is_empty() {
        return (vec![], 0, 0);
    }

    let mut cigar = Vec::with_capacity(32);

    let diag_sum: i64 = hits
        .iter()
        .map(|h| (h.ref_pos as i64) - (h.read_pos as i64))
        .sum();
    let diag = (diag_sum / hits.len() as i64) as i32;

    let read_len = read_seq.len();
    let ref_len = ref_seq.len();

    let (leading_clip, ref_start) = if diag < 0 {
        ((-diag) as usize, 0usize)
    } else {
        (0usize, diag as usize)
    };

    if leading_clip >= read_len || ref_start >= ref_len {
        push_cigar(&mut cigar, read_len as u32, 4);
        return (cigar, ref_start as i32, 0);
    }

    let max_align_len = read_len - leading_clip;
    let ref_avail = ref_len.saturating_sub(ref_start);
    let aligned_len = max_align_len.min(ref_avail);
    let trailing_clip = read_len - leading_clip - aligned_len;

    if leading_clip > 0 {
        push_cigar(&mut cigar, leading_clip as u32, 4);
    }
    if aligned_len > 0 {
        push_cigar(&mut cigar, aligned_len as u32, 0);
    }
    if trailing_clip > 0 {
        push_cigar(&mut cigar, trailing_clip as u32, 4);
    }

    (cigar, ref_start as i32, aligned_len)
}

fn compute_metrics(read_seq: &[u8], ref_seq: &[u8], start_pos: i32, cigar: &[u32]) -> (i32, String, i32) {
    let mut nm = 0;
    let mut as_score = 0;

    let (match_score, mismatch_pen) = scoring_params();

    let mut r_idx = start_pos as usize;
    let mut q_idx = 0;

    for &c in cigar {
        let len = (c >> 4) as usize;
        let op = c & 0xF;

        match op {
            0 | 7 => {
                for _ in 0..len {
                    if r_idx < ref_seq.len() && q_idx < read_seq.len() {
                        let r = ref_seq[r_idx];
                        let q = read_seq[q_idx];
                        if r.eq_ignore_ascii_case(&q) {
                            as_score += match_score;
                        } else {
                            nm += 1;
                            as_score += mismatch_pen;
                        }
                    } else {
                        nm += 1;
                        as_score += mismatch_pen;
                    }
                    r_idx += 1;
                    q_idx += 1;
                }
            }
            8 => {
                nm += len as i32;
                as_score += len as i32 * mismatch_pen;
                r_idx += len;
                q_idx += len;
            }
            1 => {
                nm += len as i32;
                as_score += len as i32 * mismatch_pen;
                q_idx += len;
            }
            2 => {
                nm += len as i32;
                as_score += len as i32 * mismatch_pen;
                r_idx += len;
            }
            4 => {
                q_idx += len;
            }
            _ => {}
        }
    }

    (nm, String::new(), as_score)
}

pub fn align<I: IndexLike>(seq: &[u8], idx: &I, state: &mut State, rev: &mut Vec<u8>) -> Option<AlignmentResult> {
    let State { mins, candidates } = state;

    candidates.clear();

    get_syncmers(seq, mins);
    mins.sort_unstable_by_key(|(h, _)| (*h >> SHIFT) as u32);
    collect_candidates(idx, mins, false, candidates);

    rev.clear();
    rev.extend(seq.iter().rev().map(|b| match b {
        b'A' => b'T',
        b'C' => b'G',
        b'G' => b'C',
        b'T' => b'A',
        b'a' => b't',
        b'c' => b'g',
        b'g' => b'c',
        b't' => b'a',
        _ => b'N',
    }));
    get_syncmers(rev, mins);
    mins.sort_unstable_by_key(|(h, _)| (*h >> SHIFT) as u32);
    collect_candidates(idx, mins, true, candidates);

    if candidates.is_empty() {
        return None;
    }

    let (c1, c2) = find_top_chains(candidates);
    let min_identity = (SYNC_S as f32) / (WINDOW as f32);
    let evaluate = |chain: &ChainRes,
                    seq: &[u8],
                    rev: &[u8],
                    candidates: &[Hit]|
     -> Option<(AlignmentResult, i32)> {
        if chain.score == 0 {
            return None;
        }
        let hits = &candidates[chain.start..chain.end];
        let ref_id = (hits[0].id_strand >> 1) as i32;
        if ref_id < 0 || (ref_id as usize) >= idx.ref_count() {
            return None;
        }
        let is_rev = (hits[0].id_strand & 1) == 1;
        let target_seq = if is_rev { rev } else { seq };
        let ref_seq = idx.ref_seq(ref_id as usize);

        let (cigar, pos, aligned_len) = gen_cigar(hits, target_seq, ref_seq);
        if pos < 0 || (pos as usize) >= ref_seq.len() {
            return None;
        }

        if aligned_len == 0 {
            return None;
        }

        let (nm, md, as_score) = compute_metrics(target_seq, ref_seq, pos, &cigar);
        let identity = 1.0 - (nm as f32 / aligned_len as f32);

        if identity < min_identity {
            return None;
        }
        Some((
            AlignmentResult { ref_id, pos, is_rev, mapq: 0, cigar, nm, md, as_score },
            chain.score,
        ))
    };

    let r1 = evaluate(&c1, seq, rev, candidates);
    let r2 = if c2.score > 0 { evaluate(&c2, seq, rev, candidates) } else { None };
    let (mut best_res, best_chain_score, second_chain_score) = match (r1, r2) {
        (Some((a1, s1)), Some((a2, s2))) => {
            if s1 >= s2 {
                (a1, s1, s2)
            } else {
                (a2, s2, s1)
            }
        }
        (Some((a1, s1)), None) => (a1, s1, 0),
        (None, Some((a2, s2))) => (a2, s2, 0),
        (None, None) => return None,
    };

    let best_count = best_chain_score.max(0) as u32;
    let second_count = second_chain_score.max(0) as u32;
    if best_count == 0 {
        best_res.mapq = 0;
    } else {
        let diff = best_count.saturating_sub(second_count) as f32;
        let frac = diff / (best_count as f32);
        best_res.mapq = (frac * (u8::MAX as f32)).round() as u8;
    }

    Some(best_res)
}

pub fn write_cigar_string(cigar: &[u32], out: &mut Vec<u8>) {
    if cigar.is_empty() {
        out.push(b'*');
        return;
    }
    for &c in cigar {
        let len = c >> 4;
        let op = match c & 0xF {
            0 => b'M',
            1 => b'I',
            2 => b'D',
            3 => b'N',
            4 => b'S',
            5 => b'H',
            6 => b'P',
            7 => b'=',
            8 => b'X',
            _ => b'?',
        };
        write!(out, "{}{}", len, op as char).unwrap();
    }
}

pub fn reverse_qual(qual: &[u8]) -> Vec<u8> {
    qual.iter().rev().cloned().collect()
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|b| match b {
            b'A' => b'T',
            b'C' => b'G',
            b'G' => b'C',
            b'T' => b'A',
            b'a' => b't',
            b'c' => b'g',
            b'g' => b'c',
            b't' => b'a',
            _ => b'N',
        })
        .collect()
}

pub fn oriented_bases<'a>(seq: &'a [u8], qual: &'a [u8], res: &Option<AlignmentResult>) -> (Cow<'a, [u8]>, Cow<'a, [u8]>) {
    let mut final_seq = Cow::Borrowed(seq);
    let mut final_qual = Cow::Borrowed(qual);
    if let Some(r) = res {
        if r.is_rev {
            final_seq = Cow::Owned(reverse_complement(seq));
            final_qual = Cow::Owned(reverse_qual(qual));
        }
    }
    (final_seq, final_qual)
}
