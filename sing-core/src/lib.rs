use anyhow::{bail, Context, Result};
use bytemuck::try_cast_slice;
use memmap2::MmapOptions;
use std::borrow::Cow;
use std::fs::File;
use std::io::{BufReader, Read, Write};
use std::ops::Range;
use std::path::Path;

#[derive(Clone, Copy)]
pub struct Config {
    pub window: usize,
    pub sync_s: usize,
    pub match_score: i32,
    pub mismatch_pen: i32,
    pub gap_open: i32,
    pub gap_ext: i32,
    pub x_drop: i32,
    pub max_hits: usize,
    pub maxindel: usize,
    pub min_identity: f32,
    pub downsample_threshold: u32,
    pub diag_band: i32,
    pub cluster_window: usize,
}

pub static CONFIG: Config = Config {
    window: 24,
    sync_s: 14,
    match_score: 2,
    mismatch_pen: -4,
    gap_open: -4,
    gap_ext: -2,
    x_drop: 50,
    max_hits: 200,
    maxindel: 20,
    min_identity: 0.8,
    downsample_threshold: 60,
    diag_band: 80,
    cluster_window: 100,
};

const WINDOW: usize = CONFIG.window;
const SYNC_S: usize = CONFIG.sync_s;
const FREQ_FILTER_ROOT: f64 = 2.5;

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
    let freq_filter = (total_bases as f64).powf(1.0 / FREQ_FILTER_ROOT) as usize;
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
    for (i, &count) in global_counts.iter().enumerate() {
        if count as usize <= freq_filter {
            kept_counts[i] = count;
            if count > 0 {
                kept_buckets += 1;
            }
        }
    }

    let downsample_threshold = CONFIG.downsample_threshold;
    let mut downsample_counts = vec![0u32; kept_counts.len()];
    let mut seen = vec![0u32; kept_counts.len()];
    for seq in &ref_seqs {
        get_syncmers(seq, &mut mins);
        for (h, _) in mins.iter() {
            let bucket = (*h as u64 >> SHIFT) as usize;
            if kept_counts[bucket] == 0 {
                continue;
            }
            let count = kept_counts[bucket];
            if count > downsample_threshold {
                let step = (count / downsample_threshold).max(2);
                let s = seen[bucket];
                if s % step == 0 {
                    downsample_counts[bucket] += 1;
                }
                seen[bucket] = s.wrapping_add(1);
            } else {
                downsample_counts[bucket] += 1;
            }
        }
    }

    
    let mut offsets = Vec::with_capacity(downsample_counts.len());
    let mut current_offset = 0u32;
    for count in &downsample_counts {
        offsets.push(current_offset);
        current_offset += *count;
    }

    
    
    let mut seeds = vec![0u64; current_offset as usize];
    let mut write_pos = offsets.clone();
    let mut seen = vec![0u32; kept_counts.len()];
    for (rid, seq) in ref_seqs.iter().enumerate() {
        get_syncmers(seq, &mut mins);
        for (h, p) in mins.iter() {
            let hash64 = *h as u64;
            let bucket = (hash64 >> SHIFT) as usize;
            if kept_counts[bucket] == 0 {
                continue;
            }
            let count = kept_counts[bucket];
            if count > downsample_threshold {
                let step = (count / downsample_threshold).max(2);
                let s = seen[bucket];
                if s % step != 0 {
                    seen[bucket] = s.wrapping_add(1);
                    continue;
                }
                seen[bucket] = s.wrapping_add(1);
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
fn collect_candidates<I: IndexLike>(
    idx: &I,
    mins: &[(u32, u32)],
    is_rev: bool,
    max_hits: usize,
    out: &mut Vec<Hit>,
) -> bool {
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
                if out.len() >= max_hits {
                    return true;
                }
            }
            i += 1;
        }
    }
    false
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

#[inline(always)]
fn find_densest_cluster(hits: &[Hit], band: i32, window: usize) -> (usize, usize, i32) {
    if hits.is_empty() {
        return (0, 0, 0);
    }
    let mut best_start = 0usize;
    let mut best_end = 0usize;
    let mut best_count = 0usize;
    let mut best_id_strand = hits[0].id_strand;
    let mut j = 0usize;
    for i in 0..hits.len() {
        while j < hits.len()
            && hits[j].id_strand == hits[i].id_strand
            && (hits[j].diag - hits[i].diag) < band
        {
            j += 1;
        }
        let cnt = j - i;
        if cnt > best_count || (cnt == best_count && hits[i].id_strand != best_id_strand && cnt >= window) {
            best_count = cnt;
            best_start = i;
            best_end = j;
            best_id_strand = hits[i].id_strand;
        }
    }
    let median_idx = best_start + best_count / 2;
    (best_start, best_end, hits.get(median_idx).map_or(0, |h| h.diag))
}

const CASE_MASK: u64 = 0x2020202020202020;

#[inline(always)]
fn extend_left(
    read: &[u8],
    rseq: &[u8],
    mut rp: usize,
    mut gp: usize,
    cfg: &Config,
) -> (i32, usize, usize, Vec<u32>) {
    let mut score = 0i32;
    let mut max_score = 0i32;
    let mut cigar = Vec::with_capacity(8);
    let mut match_run = 0u32;
    let mut best_rp = rp;
    let mut best_gp = gp;
    let mut best_cigar_len = 0;
    let mut best_match_run = 0;

    while rp > 0 && gp > 0 {
        if rp >= 8 && gp >= 8 {
            unsafe {
                let r_val = std::ptr::read_unaligned(read.as_ptr().add(rp - 8) as *const u64) | CASE_MASK;
                let g_val = std::ptr::read_unaligned(rseq.as_ptr().add(gp - 8) as *const u64) | CASE_MASK;
                if r_val == g_val {
                    score += cfg.match_score * 8;
                    match_run += 8;
                    rp -= 8;
                    gp -= 8;
                    if score > max_score {
                        max_score = score;
                        best_rp = rp;
                        best_gp = gp;
                        best_cigar_len = cigar.len();
                        best_match_run = match_run;
                    }
                    continue;
                }
                
                let x = r_val ^ g_val;
                let match_bytes = (x.leading_zeros() / 8) as usize;
                if match_bytes > 0 {
                    score += cfg.match_score * match_bytes as i32;
                    match_run += match_bytes as u32;
                    rp -= match_bytes;
                    gp -= match_bytes;
                }
            }
        }

        let rb = read[rp - 1];
        let gb = rseq[gp - 1];
        if rb.eq_ignore_ascii_case(&gb) {
            score += cfg.match_score;
            match_run += 1;
            rp -= 1;
            gp -= 1;
        } else {
            if match_run > 0 {
                cigar.push((match_run << 4) | 0);
                match_run = 0;
            }
            let mm_score = score + cfg.mismatch_pen;
            let del_score = score + cfg.gap_open;
            let ins_score = score + cfg.gap_open;
            if mm_score >= del_score && mm_score >= ins_score {
                score = mm_score;
                cigar.push((1 << 4) | 0);
                rp -= 1;
                gp -= 1;
            } else if del_score >= ins_score && gp > 0 {
                score = del_score;
                let mut glen = 1usize;
                while glen < CONFIG.maxindel && gp > glen && score + cfg.gap_ext > max_score - cfg.x_drop {
                    if gp > glen && rp > 0 && read[rp - 1].eq_ignore_ascii_case(&rseq[gp - glen - 1]) {
                        break;
                    }
                    glen += 1;
                    score += cfg.gap_ext;
                }
                cigar.push(((glen as u32) << 4) | 2);
                gp = gp.saturating_sub(glen);
            } else if rp > 0 {
                score = ins_score;
                cigar.push((1 << 4) | 1);
                rp -= 1;
            } else {
                break;
            }
        }
        if score > max_score {
            max_score = score;
            best_rp = rp;
            best_gp = gp;
            best_cigar_len = cigar.len();
            best_match_run = match_run;
        }
        if score < max_score - cfg.x_drop {
            break;
        }
    }
    cigar.truncate(best_cigar_len);
    if best_match_run > 0 {
        cigar.push((best_match_run << 4) | 0);
    }
    cigar.reverse();
    (max_score, best_rp, best_gp, cigar)
}

#[inline(always)]
fn extend_right(
    read: &[u8],
    rseq: &[u8],
    mut rp: usize,
    mut gp: usize,
    cfg: &Config,
) -> (i32, usize, usize, Vec<u32>) {
    let rlen = read.len();
    let glen = rseq.len();
    let mut score = 0i32;
    let mut max_score = 0i32;
    let mut cigar = Vec::with_capacity(8);
    let mut match_run = 0u32;
    let mut best_rp = rp;
    let mut best_gp = gp;
    let mut best_cigar_len = 0;
    let mut best_match_run = 0;

    while rp < rlen && gp < glen {
        if rlen - rp >= 8 && glen - gp >= 8 {
             unsafe {
                let r_val = std::ptr::read_unaligned(read.as_ptr().add(rp) as *const u64) | CASE_MASK;
                let g_val = std::ptr::read_unaligned(rseq.as_ptr().add(gp) as *const u64) | CASE_MASK;
                if r_val == g_val {
                    rp += 8;
                    gp += 8;
                    score += cfg.match_score * 8;
                    match_run += 8;
                    if score > max_score {
                        max_score = score;
                        best_rp = rp;
                        best_gp = gp;
                        best_cigar_len = cigar.len();
                        best_match_run = match_run;
                    }
                    continue;
                }

                let x = r_val ^ g_val;
                let match_bytes = (x.trailing_zeros() / 8) as usize;
                if match_bytes > 0 {
                    score += cfg.match_score * match_bytes as i32;
                    match_run += match_bytes as u32;
                    rp += match_bytes;
                    gp += match_bytes;
                }
             }
        }

        let rb = read[rp];
        let gb = rseq[gp];
        if rb.eq_ignore_ascii_case(&gb) {
            score += cfg.match_score;
            match_run += 1;
            rp += 1;
            gp += 1;
        } else {
            if match_run > 0 {
                cigar.push((match_run << 4) | 0);
                match_run = 0;
            }
            let mm_score = score + cfg.mismatch_pen;
            let del_score = score + cfg.gap_open;
            let ins_score = score + cfg.gap_open;
            if mm_score >= del_score && mm_score >= ins_score {
                score = mm_score;
                cigar.push((1 << 4) | 0);
                rp += 1;
                gp += 1;
            } else if del_score >= ins_score && gp < glen {
                score = del_score;
                let mut dlen = 1usize;
                while dlen < 8 && gp + dlen < glen && score + cfg.gap_ext > max_score - cfg.x_drop {
                    if rp < rlen && read[rp].eq_ignore_ascii_case(&rseq[gp + dlen]) {
                        break;
                    }
                    dlen += 1;
                    score += cfg.gap_ext;
                }
                cigar.push(((dlen as u32) << 4) | 2);
                gp += dlen;
            } else if rp < rlen {
                score = ins_score;
                cigar.push((1 << 4) | 1);
                rp += 1;
            } else {
                break;
            }
        }
        if score > max_score {
            max_score = score;
            best_rp = rp;
            best_gp = gp;
            best_cigar_len = cigar.len();
            best_match_run = match_run;
        }
        if score < max_score - cfg.x_drop {
            break;
        }
    }
    cigar.truncate(best_cigar_len);
    if best_match_run > 0 {
        cigar.push((best_match_run << 4) | 0);
    }
    (max_score, best_rp, best_gp, cigar)
}

#[inline(always)]
fn greedy_extend(
    hits: &[Hit],
    anchor_diag: i32,
    read: &[u8],
    rseq: &[u8],
    cfg: &Config,
) -> (Vec<u32>, i32, i32, usize) {
    let rlen = read.len();
    let anchor_idx = hits.len() / 2;
    let anchor = &hits[anchor_idx.min(hits.len().saturating_sub(1))];
    let a_rp = anchor.read_pos as usize;
    let mut a_gp = if anchor_diag >= 0 {
        anchor_diag as usize + a_rp
    } else {
        a_rp.saturating_sub((-anchor_diag) as usize)
    };
    if a_gp >= rseq.len() {
        a_gp = rseq.len().saturating_sub(1);
    }
    let (l_score, l_rp, l_gp, l_cigar) = extend_left(read, rseq, a_rp, a_gp, cfg);
    let (r_score, r_rp, _, r_cigar) = extend_right(read, rseq, a_rp, a_gp, cfg);
    let mut cigar = Vec::with_capacity(l_cigar.len() + r_cigar.len() + 2);
    if l_rp > 0 {
        cigar.push(((l_rp as u32) << 4) | 4);
    }
    for &c in &l_cigar {
        push_cigar(&mut cigar, c >> 4, c & 0xF);
    }
    for &c in &r_cigar {
        push_cigar(&mut cigar, c >> 4, c & 0xF);
    }
    let trailing = rlen.saturating_sub(r_rp);
    if trailing > 0 {
        push_cigar(&mut cigar, trailing as u32, 4);
    }
    let aligned = r_rp.saturating_sub(l_rp);
    let as_score = l_score + r_score;
    (cigar, l_gp as i32, as_score, aligned)
}

fn compute_nm(read: &[u8], rseq: &[u8], pos: i32, cigar: &[u32]) -> i32 {
    let mut nm = 0i32;
    let mut ri = pos as usize;
    let mut qi = 0usize;
    for &c in cigar {
        let len = (c >> 4) as usize;
        let op = c & 0xF;
        match op {
            0 => {
                for _ in 0..len {
                    if ri < rseq.len() && qi < read.len() && !read[qi].eq_ignore_ascii_case(&rseq[ri]) {
                        nm += 1;
                    }
                    ri += 1;
                    qi += 1;
                }
            }
            1 => { nm += len as i32; qi += len; }
            2 => { nm += len as i32; ri += len; }
            4 => { qi += len; }
            _ => {}
        }
    }
    nm
}

pub fn align<I: IndexLike>(seq: &[u8], idx: &I, state: &mut State, rev: &mut Vec<u8>) -> Vec<AlignmentResult> {
    let State { mins, candidates } = state;
    let cfg = &CONFIG;
    candidates.clear();
    let max_hits = cfg.max_hits;
    let mut repetitive = false;
    
    let offsets = idx.offsets();
    let offsets_len = offsets.len();
    let seeds_len = idx.seeds().len() as u32;

    let sort_mins = |mins: &mut Vec<(u32, u32)>| {
        mins.sort_unstable_by_key(|&(h, _)| {
            let bucket = (h >> SHIFT) as usize;
            let start = unsafe { *offsets.get_unchecked(bucket) };
            let end = if bucket + 1 < offsets_len {
                unsafe { *offsets.get_unchecked(bucket + 1) }
            } else {
                seeds_len
            };
            end.wrapping_sub(start)
        });
    };

    get_syncmers(seq, mins);
    sort_mins(mins);
    if collect_candidates(idx, mins, false, max_hits, candidates) {
        repetitive = true;
    }
    rev.clear();
    rev.extend(seq.iter().rev().map(|b| match b {
        b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A',
        b'a' => b't', b'c' => b'g', b'g' => b'c', b't' => b'a',
        _ => b'N',
    }));
    get_syncmers(rev, mins);
    sort_mins(mins);
    if collect_candidates(idx, mins, true, max_hits, candidates) {
        repetitive = true;
    }
    if candidates.is_empty() {
        return Vec::new();
    }
    candidates.sort_unstable_by(|a, b| (a.id_strand, a.diag).cmp(&(b.id_strand, b.diag)));
    let total_seeds = candidates.len();
    let (c1_start, c1_end, anchor_diag) = find_densest_cluster(candidates, cfg.diag_band, cfg.cluster_window);
    let c1_count = c1_end - c1_start;
    if c1_count < 3 {
        return Vec::new();
    }
    let (pre_start, pre_end, pre_diag) = find_densest_cluster(&candidates[..c1_start], cfg.diag_band, cfg.cluster_window);
    let (post_start, post_end, post_diag) = find_densest_cluster(&candidates[c1_end..], cfg.diag_band, cfg.cluster_window);
    let pre_count = pre_end.saturating_sub(pre_start);
    let post_count = post_end.saturating_sub(post_start);
    let (second_range, second_count) = if pre_count >= post_count {
        if pre_count > 0 {
            ((pre_start, pre_end, pre_diag), pre_count)
        } else {
            ((0, 0, 0), 0)
        }
    } else {
        if post_count > 0 {
            ((c1_end + post_start, c1_end + post_end, post_diag), post_count)
        } else {
            ((0, 0, 0), 0)
        }
    };

    let attempt = |range: (usize, usize, i32)| -> Option<(AlignmentResult, usize)> {
        let (s, e, diag) = range;
        if e <= s {
            return None;
        }
        let hits = &candidates[s..e];
        let ref_id = (hits[0].id_strand >> 1) as i32;
        if (ref_id as usize) >= idx.ref_count() {
            return None;
        }
        let is_rev = (hits[0].id_strand & 1) == 1;
        let target_seq = if is_rev { rev.as_slice() } else { seq };
        let ref_seq = idx.ref_seq(ref_id as usize);
        let (cigar, pos, as_score, aligned_len) = greedy_extend(hits, diag, target_seq, ref_seq, cfg);
        if aligned_len == 0 || pos < 0 || (pos as usize) >= ref_seq.len() {
            return None;
        }
        let nm = compute_nm(target_seq, ref_seq, pos, &cigar);
        let identity = 1.0 - (nm as f32 / aligned_len.max(1) as f32);
        if identity < cfg.min_identity {
            return None;
        }
        Some((AlignmentResult { ref_id, pos, is_rev, mapq: 0, cigar, nm, md: String::new(), as_score }, e - s))
    };

    let mut results = Vec::with_capacity(2);

    let mut other_count = second_count;
    if let Some((mut res, primary_len)) = attempt((c1_start, c1_end, anchor_diag)) {
        let mapq = if repetitive {
            0u8
        } else {
            let frac = (primary_len.saturating_sub(other_count)) as f32 / total_seeds.max(1) as f32;
            (frac * 60.0).round().min(60.0) as u8
        };
        res.mapq = mapq;
        results.push(res);
    }

    if second_count > 0 {
        if let Some((mut res, primary_len)) = attempt(second_range) {
            other_count = c1_count;
            let mapq = if repetitive {
                0u8
            } else {
                let frac = (primary_len.saturating_sub(other_count)) as f32 / total_seeds.max(1) as f32;
                (frac * 60.0).round().min(60.0) as u8
            };
            res.mapq = mapq;
            results.push(res);
        }
    }

    results.sort_by(|a, b| b.as_score.cmp(&a.as_score));
    results
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
