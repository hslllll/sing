use anyhow::{bail, Context, Result};
use bytemuck::try_cast_slice;
use memmap2::MmapOptions;
use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufReader, Read, Write};
use std::ops::Range;
use std::path::Path;

#[derive(Clone, Copy)]
pub struct Config {
    pub minimizer_window: usize,
    pub hash_window: usize,
    pub sync_s: usize,
    pub match_score: i32,
    pub mismatch_pen: i32,
    pub gap_open: i32,
    pub gap_ext: i32,
    pub x_drop: i32,
    pub max_hits: usize,
    pub maxindel: usize,
    pub min_identity: f32,
    pub pair_max_dist: i32,
}

pub static CONFIG: Config = Config {
    minimizer_window: 19,
    hash_window: 16,
    sync_s: 14,
    match_score: 2,
    mismatch_pen: -2,
    gap_open: -2,
    gap_ext: -1,
    x_drop: 15,
    max_hits: 8000,
    pair_max_dist: 1000,
    maxindel: 15,
    min_identity: 0.85,
};

const HASH_WINDOW: usize = CONFIG.hash_window;
const MINIMIZER_WINDOW: usize = CONFIG.minimizer_window;
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

pub const RADIX: usize = 24;
pub const SHIFT: usize = 32 - RADIX;

#[derive(Debug, Clone)]
pub struct Index {
    pub offsets: Vec<u32>,
    pub seeds: Vec<u64>,
    pub genome_size: u64,
    pub ref_seqs: Vec<Vec<u8>>,
    pub ref_names: Vec<String>,
}

pub trait IndexLike {
    fn offsets(&self) -> &[u32];
    fn seeds(&self) -> &[u64];
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

        Ok(Index { offsets, seeds, genome_size, ref_seqs, ref_names })
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

        Ok(Index { offsets, seeds, genome_size, ref_seqs, ref_names })
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
    pub paired: bool,
    pub proper_pair: bool,
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
    let max_hits = CONFIG.max_hits;
    
    eprintln!("Building sort-based index: {} bp, max_hits={}", total_bases, max_hits);

    
    let mut mins = Vec::with_capacity(1024);
    let mut counts: HashMap<u32, u32> = HashMap::new();

    for seq in ref_seqs.iter() {
        get_syncmers(seq, &mut mins);
        for &(h, _p, _is_rev) in mins.iter() {
            let entry = counts.entry(h).or_insert(0);
            if *entry <= max_hits as u32 {
                *entry += 1;
            }
        }
    }

    eprintln!("  Counted {} unique seeds", counts.len());

    let mut kept_seeds = Vec::new();
    let mut kept = 0usize;
    let mut filtered = 0usize;

    for (rid, seq) in ref_seqs.iter().enumerate() {
        get_syncmers(seq, &mut mins);
        for &(h, p, _is_rev) in mins.iter() {
            match counts.get(&h) {
                Some(&c) if c <= max_hits as u32 => {
                    kept_seeds.push((h, rid as u32, p));
                    kept += 1;
                }
                Some(_) => {
                    filtered += 1;
                }
                None => {}
            }
        }
    }

    eprintln!("  Kept {} seeds, filtered {} seeds", kept, filtered);

    
    let mut offsets = vec![0u32; 1 << RADIX];
    let mut final_seeds = Vec::with_capacity(kept_seeds.len());
    
    
    kept_seeds.sort_unstable_by_key(|k| k.0 >> SHIFT);
    
    let mut write_offset = 0u32;
    let mut last_bucket = 0usize;
    
    for &(h, rid, pos) in &kept_seeds {
        let bucket = (h >> SHIFT) as usize;
        
        
        while last_bucket <= bucket {
            offsets[last_bucket] = write_offset;
            last_bucket += 1;
        }
        
        let hash16 = (h & 0xFFFF) as u64;
        let seed = (hash16 << 48) | ((rid as u64) << 32) | (pos as u64);
        final_seeds.push(seed);
        write_offset += 1;
    }
    
    
    while last_bucket < offsets.len() {
        offsets[last_bucket] = write_offset;
        last_bucket += 1;
    }
    
    eprintln!("  Index built: {} seeds", final_seeds.len());

    Ok(Index {
        offsets,
        seeds: final_seeds,
        genome_size: total_bases,
        ref_seqs,
        ref_names,
    })
}

#[derive(Clone, Copy, Default, Debug)]
pub struct Hit {
    pub id_strand: u32,
    pub diag: i32,
    pub read_pos: u32,
    pub ref_pos: u32,
}

pub struct State {
    pub mins: Vec<(u32, u32, bool)>,
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
pub fn get_syncmers(seq: &[u8], out: &mut Vec<(u32, u32, bool)>) {
    out.clear();
    if seq.len() < HASH_WINDOW || MINIMIZER_WINDOW == 0 || SYNC_S > HASH_WINDOW {
        return;
    }

    let mut s_queue: Vec<(u32, u64)> = Vec::new();
    let mut s_head: usize = 0;

    let mut h_k = 0u64;
    let mut h_k_rc = 0u64;
    let mut h_s = 0u64;
    let mut h_s_rc = 0u64;
    let mut base_buf_k = [0usize; HASH_WINDOW];
    let mut base_buf_s = [0usize; SYNC_S];
    let mut ambig_buf_k = [0u8; HASH_WINDOW];
    let mut ambig_buf_s = [0u8; SYNC_S];
    let mut ambig_k = 0i32;
    let mut ambig_s = 0i32;

    for i in 0..seq.len() {
        let (base_idx, ambig) = match base_to_index(seq[i]) {
            Some(idx) => (idx, 0u8),
            None => (0usize, 1u8),
        };
        
        let rc_idx = 3 - base_idx;

        let k_slot = i % HASH_WINDOW;
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

        if i + 1 <= HASH_WINDOW {
            h_k = h_k.rotate_left(ROT) ^ BASES[base_idx];
            h_k_rc = (BASES[rc_idx]).rotate_right(ROT) ^ h_k_rc;
        } else {
            h_k = h_k.rotate_left(ROT) ^ BASES[base_idx] ^ REMOVE[prev_idx_k];
            h_k_rc = (BASES[rc_idx] ^ REMOVE[3 - prev_idx_k]).rotate_right(ROT) ^ h_k_rc;
        }

        if i + 1 <= SYNC_S {
            h_s = h_s.rotate_left(ROT) ^ BASES[base_idx];
            h_s_rc = (BASES[rc_idx]).rotate_right(ROT) ^ h_s_rc;
        } else {
            h_s = h_s.rotate_left(ROT) ^ BASES[base_idx] ^ REMOVE_S[prev_idx_s];
            h_s_rc = (BASES[rc_idx] ^ REMOVE_S[3 - prev_idx_s]).rotate_right(ROT) ^ h_s_rc;
        }

        if i + 1 >= SYNC_S && ambig_s == 0 {
            let s_pos = (i + 1 - SYNC_S) as u32;
            let s_hash = h_s.wrapping_mul(0x517cc1b727220a95);
            let s_hash_rc = h_s_rc.wrapping_mul(0x517cc1b727220a95);
            let canonical_hash = s_hash.min(s_hash_rc);

            if s_head >= s_queue.len() {
                s_queue.clear();
                s_head = 0;
            }

            while let Some(&(_, back_hash)) = s_queue.last() {
                if back_hash >= canonical_hash {
                    s_queue.pop();
                } else {
                    break;
                }
            }
            if s_head >= s_queue.len() {
                s_queue.clear();
                s_head = 0;
            }
            s_queue.push((s_pos, canonical_hash));
        }

        if i + 1 >= HASH_WINDOW {
            if ambig_k != 0 {
                continue;
            }
            let k_pos = (i + 1 - HASH_WINDOW) as u32;

            // Maintain s-mer minima within the current k-mer window.
            while s_head < s_queue.len() {
                let pos = s_queue[s_head].0;
                if pos < k_pos {
                    s_head += 1;
                } else {
                    break;
                }
            }
            let max_pos = k_pos + (HASH_WINDOW - SYNC_S) as u32;
            while s_head < s_queue.len() {
                let pos = s_queue[s_head].0;
                if pos < k_pos {
                    s_head += 1;
                } else if pos > max_pos {
                    break;
                } else {
                    break;
                }
            }

            if s_head >= s_queue.len() {
                s_queue.clear();
                s_head = 0;
                continue;
            }
            let s_min_pos = s_queue[s_head].0;

            if s_head > 1024 && s_head * 2 > s_queue.len() {
                s_queue.drain(0..s_head);
                s_head = 0;
            }

            if (k_pos + 1) >= MINIMIZER_WINDOW as u32 {
                let k_hash = h_k.wrapping_mul(0x517cc1b727220a95);
                let k_hash_rc = h_k_rc.wrapping_mul(0x517cc1b727220a95);
                let is_rev = k_hash_rc < k_hash;
                let canonical_k = k_hash.min(k_hash_rc);

                let is_closed = s_min_pos == k_pos || s_min_pos == max_pos;
                if is_closed {
                    if out.last().map(|&(_, p, _)| p) != Some(k_pos) {
                        out.push((canonical_k as u32, k_pos, is_rev));
                    }
                }
            }
        }
    }
}

#[inline(always)]
fn collect_candidates<I: IndexLike>(
    idx: &I,
    mins: &[(u32, u32, bool)],
    max_hits: usize,
    out: &mut Vec<Hit>,
) -> bool {
    let offsets = idx.offsets();
    let seeds = idx.seeds();
    let mut capped = false;
    
    for &(h, r_pos, is_rev) in mins {
        let bucket = (h >> SHIFT) as usize;
        
        if bucket >= offsets.len() {
            continue;
        }
        
        let start = offsets[bucket] as usize;
        let end = if bucket + 1 < offsets.len() {
            offsets[bucket + 1] as usize
        } else {
            seeds.len()
        };
        
        if start >= end || end > seeds.len() {
            continue;
        }
        
        let target_hash = (h & 0xFFFF) as u64;
        let bucket_seeds = &seeds[start..end];
        
        
        let pos = bucket_seeds.partition_point(|&s| (s >> 48) < target_hash);
        
        if pos >= bucket_seeds.len() {
            continue;
        }
        
        let mut count = 0;
        let mut idx_in_bucket = pos;
        
        while idx_in_bucket < bucket_seeds.len() {
            let seed = bucket_seeds[idx_in_bucket];
            let seed_hash = seed >> 48;
            
            if seed_hash != target_hash {
                break;
            }

            if count >= max_hits {
                capped = true;
                break;
            }
            
            let ref_id = ((seed >> 32) & 0xFFFF) as u32;
            let ref_pos = (seed & 0xFFFFFFFF) as u32;
            let id_strand = (ref_id << 1) | (if is_rev { 1 } else { 0 });
            let diag = (ref_pos as i32).wrapping_sub(r_pos as i32);
            
            out.push(Hit {
                id_strand,
                diag,
                read_pos: r_pos,
                ref_pos,
            });
            
            count += 1;
            idx_in_bucket += 1;
        }
    }
    
    capped
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
fn best_diag_for_group(group: &[Hit]) -> Option<(i32, usize)> {
    if group.is_empty() {
        return None;
    }
    let mut counts: Vec<(i32, usize)> = Vec::new();
    let mut best_diag = 0i32;
    let mut best_count = 0usize;
    for hit in group {
        let mut found = false;
        for (d, c) in counts.iter_mut() {
            if *d == hit.diag {
                *c += 1;
                if *c > best_count {
                    best_count = *c;
                    best_diag = *d;
                }
                found = true;
                break;
            }
        }
        if !found {
            counts.push((hit.diag, 1));
            if best_count == 0 {
                best_count = 1;
                best_diag = hit.diag;
            }
        }
    }
    Some((best_diag, best_count))
}

#[inline(always)]
fn anchor_for_diag(group: &[Hit], diag: i32) -> Option<(i32, u32)> {
    if group.is_empty() {
        return None;
    }
    let mut matches: Vec<usize> = Vec::new();
    for (i, hit) in group.iter().enumerate() {
        if hit.diag == diag {
            matches.push(i);
        }
    }
    if matches.is_empty() {
        return None;
    }
    let mid = matches.len() / 2;
    let hit = &group[matches[mid]];
    Some((diag, hit.read_pos))
}

const CASE_MASK: u64 = 0x2020202020202020;

#[inline(always)]
fn normalize_n(b: u8) -> u8 {
    match b {
        b'N' | b'n' => b'A',
        _ => b,
    }
}

#[inline(always)]
fn eq_for_align(a: u8, b: u8) -> bool {
    normalize_n(a).eq_ignore_ascii_case(&normalize_n(b))
}

#[inline(always)]
fn has_n_u64(v: u64) -> bool {
    v.to_le_bytes().iter().any(|&b| b == b'n')
}

#[inline(always)]
fn extend_left(
    read: &[u8],
    rseq: &[u8],
    mut rp: usize,
    mut gp: usize,
    cfg: &Config,
) -> (i32, usize, usize, Vec<u32>, usize) {
    let mut score = 0i32;
    let mut max_score = 0i32;
    let mut cigar = Vec::with_capacity(8);
    let mut match_run = 0u32;
    let mut best_rp = rp;
    let mut best_gp = gp;
    let mut best_cigar_len = 0;
    let mut best_match_run = 0;
    let mut nm = 0usize;
    let mut best_nm = 0usize;

    while rp > 0 && gp > 0 {
        if rp >= 8 && gp >= 8 {
            unsafe {
                let r_val = std::ptr::read_unaligned(read.as_ptr().add(rp - 8) as *const u64) | CASE_MASK;
                let g_val = std::ptr::read_unaligned(rseq.as_ptr().add(gp - 8) as *const u64) | CASE_MASK;
                if !has_n_u64(r_val) && !has_n_u64(g_val) {
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
                            best_nm = nm;
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
                    
                    let remaining = 8 - match_bytes;
                    if remaining > 0 {
                        let x_masked = x << (match_bytes * 8);
                        for i in 0..remaining {
                            let byte_idx = 7 - i;
                            let byte_mask = 0xFFu64 << (byte_idx * 8);
                            if (x_masked & byte_mask) != 0 {
                                nm += 1;
                            }
                        }
                    }
                }
            }
        }

        let rb = read[rp - 1];
        let gb = rseq[gp - 1];
        if eq_for_align(rb, gb) {
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
                nm += 1;
                rp -= 1;
                gp -= 1;
            } else if del_score >= ins_score && gp > 0 {
                score = del_score;
                let mut glen = 1usize;
                nm += 1;
                while glen < CONFIG.maxindel && gp > glen && score + cfg.gap_ext > max_score - cfg.x_drop {
                    if gp > glen && rp > 0 && eq_for_align(read[rp - 1], rseq[gp - glen - 1]) {
                        break;
                    }
                    glen += 1;
                    nm += 1;
                    score += cfg.gap_ext;
                }
                cigar.push(((glen as u32) << 4) | 2);
                gp = gp.saturating_sub(glen);
            } else if rp > 0 {
                score = ins_score;
                cigar.push((1 << 4) | 1);
                nm += 1;
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
            best_nm = nm;
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
    (max_score, best_rp, best_gp, cigar, best_nm)
}

#[inline(always)]
fn extend_right(
    read: &[u8],
    rseq: &[u8],
    mut rp: usize,
    mut gp: usize,
    cfg: &Config,
) -> (i32, usize, usize, Vec<u32>, usize) {
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
    let mut nm = 0usize;
    let mut best_nm = 0usize;

    while rp < rlen && gp < glen {
        if rlen - rp >= 8 && glen - gp >= 8 {
             unsafe {
                let r_val = std::ptr::read_unaligned(read.as_ptr().add(rp) as *const u64) | CASE_MASK;
                let g_val = std::ptr::read_unaligned(rseq.as_ptr().add(gp) as *const u64) | CASE_MASK;
                if !has_n_u64(r_val) && !has_n_u64(g_val) {
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
                            best_nm = nm;
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
                    
                    let remaining = 8 - match_bytes;
                    if remaining > 0 {
                        let x_masked = x >> (match_bytes * 8);
                        for i in 0..remaining {
                            let byte_mask = 0xFFu64 << (i * 8);
                            if (x_masked & byte_mask) != 0 {
                                nm += 1;
                            }
                        }
                    }
                }
             }
        }

        let rb = read[rp];
        let gb = rseq[gp];
        if eq_for_align(rb, gb) {
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
                nm += 1;
                rp += 1;
                gp += 1;
            } else if del_score >= ins_score && gp < glen {
                score = del_score;
                let mut dlen = 1usize;
                nm += 1;
                while dlen < 8 && gp + dlen < glen && score + cfg.gap_ext > max_score - cfg.x_drop {
                    if rp < rlen && eq_for_align(read[rp], rseq[gp + dlen]) {
                        break;
                    }
                    dlen += 1;
                    nm += 1;
                    score += cfg.gap_ext;
                }
                cigar.push(((dlen as u32) << 4) | 2);
                gp += dlen;
            } else if rp < rlen {
                score = ins_score;
                cigar.push((1 << 4) | 1);
                nm += 1;
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
            best_nm = nm;
        }
        if score < max_score - cfg.x_drop {
            break;
        }
    }
    cigar.truncate(best_cigar_len);
    if best_match_run > 0 {
        cigar.push((best_match_run << 4) | 0);
    }
    (max_score, best_rp, best_gp, cigar, best_nm)
}

#[inline(always)]
fn greedy_extend(
    anchor_read_pos: usize,
    anchor_diag: i32,
    read: &[u8],
    rseq: &[u8],
    cfg: &Config,
) -> (Vec<u32>, i32, i32, usize, usize) {
    let rlen = read.len();
    let a_rp = anchor_read_pos;
    let mut a_gp = if anchor_diag >= 0 {
        anchor_diag as usize + a_rp
    } else {
        a_rp.saturating_sub((-anchor_diag) as usize)
    };
    if a_gp >= rseq.len() {
        a_gp = rseq.len().saturating_sub(1);
    }
    let (l_score, l_rp, l_gp, l_cigar, l_nm) = extend_left(read, rseq, a_rp, a_gp, cfg);
    let (r_score, r_rp, _, r_cigar, r_nm) = extend_right(read, rseq, a_rp, a_gp, cfg);
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
    let nm = l_nm + r_nm;
    (cigar, l_gp as i32, as_score, aligned, nm)
}

#[inline(always)]
fn span_check(hits: &[Hit], cfg: &Config) -> bool {
    if hits.is_empty() {
        return false;
    }
    let mut min_read = hits[0].read_pos as i32;
    let mut max_read = hits[0].read_pos as i32;
    let mut min_ref = hits[0].ref_pos as i32;
    let mut max_ref = hits[0].ref_pos as i32;

    for hit in hits.iter().skip(1) {
        let rp = hit.read_pos as i32;
        let gp = hit.ref_pos as i32;
        if rp < min_read { min_read = rp; }
        if rp > max_read { max_read = rp; }
        if gp < min_ref { min_ref = gp; }
        if gp > max_ref { max_ref = gp; }
    }

    let query_span = max_read - min_read;
    let ref_span = max_ref - min_ref;
    
    let span_diff = (ref_span - query_span).abs();
    span_diff < cfg.maxindel as i32
}

#[inline(always)]
fn is_collinear(hits: &[Hit], cfg: &Config) -> bool {
    if hits.len() < 2 { return true; }

    let first_seed = &hits[0];
    
    let base_diag = (first_seed.ref_pos as i64) - (first_seed.read_pos as i64);
    
    let mut min_diff = 0;
    let mut max_diff = 0;

    for seed in hits.iter().skip(1) {
        let cur_diag = (seed.ref_pos as i64) - (seed.read_pos as i64);
        let diff = cur_diag - base_diag;

        if diff < min_diff { min_diff = diff; }
        if diff > max_diff { max_diff = diff; }

        if max_diff - min_diff > cfg.maxindel as i64 {
            return false;
        }
    }
    true
}

pub fn align<I: IndexLike>(seq: &[u8], idx: &I, state: &mut State, rev: &mut Vec<u8>) -> Vec<AlignmentResult> {
    let State { mins, candidates } = state;
    let cfg = &CONFIG;
    candidates.clear();
    let max_hits = cfg.max_hits;
    let mut capped = false;
    let mut second_failed = false;

    let mut tmp_mins: Vec<(u32, u32, bool)> = Vec::with_capacity(mins.capacity());

    get_syncmers(seq, &mut tmp_mins);
    mins.clear();
    mins.extend(tmp_mins.iter().map(|&(h, p, _)| (h, p, false)));

    let rev_seq = reverse_complement(seq);
    tmp_mins.clear();
    get_syncmers(&rev_seq, &mut tmp_mins);
    mins.extend(tmp_mins.iter().map(|&(h, p, _)| (h, p, true)));
    
    capped |= collect_candidates(idx, mins, max_hits, candidates);
    
    if candidates.is_empty() {
        return Vec::new();
    }

    let total_seeds = candidates.len();
    let mut groups: HashMap<usize, Vec<Hit>> = HashMap::new();
    for hit in candidates.iter().copied() {
        let id = hit.id_strand as usize;
        groups.entry(id).or_default().push(hit);
    }


    let mut group_votes: HashMap<usize, (i32, usize)> = HashMap::new();
    for (id, group) in groups.iter() {
        if let Some((diag, count)) = best_diag_for_group(group) {
            group_votes.insert(*id, (diag, count));
        }
    }

    let mut best_id: Option<usize> = None;
    let mut best_count = 0usize;
    let mut best_diag = 0i32;
    for (id, (diag, count)) in group_votes.iter() {
        if *count > best_count || (*count == best_count && best_id.map_or(true, |bid| *id > bid)) {
            best_id = Some(*id);
            best_count = *count;
            best_diag = *diag;
        }
    }

    let best_id = match best_id {
        Some(id) => id,
        None => return Vec::new(),
    };

    if best_count < 3 {
        return Vec::new();
    }

    let best_hits = match groups.get(&best_id) {
        Some(hits) => hits,
        None => return Vec::new(),
    };
    let mut best_diag_hits: Vec<Hit> = Vec::new();
    for &hit in best_hits {
        if hit.diag == best_diag {
            best_diag_hits.push(hit);
        }
    }
    if best_diag_hits.is_empty() {
        return Vec::new();
    }
    let (anchor_diag, anchor_read_pos) = match anchor_for_diag(best_hits, best_diag) {
        Some(v) => v,
        None => return Vec::new(),
    };

    if !span_check(&best_diag_hits, cfg) {
        return Vec::new();
    }

    if !is_collinear(&best_diag_hits, cfg) {
        return Vec::new();
    }

    let mut second_id: Option<usize> = None;
    let mut second_count = 0usize;
    let mut second_diag = 0i32;
    for (id, (diag, count)) in group_votes.iter() {
        if *id == best_id {
            continue;
        }
        if *count > second_count || (*count == second_count && second_id.map_or(true, |sid| *id > sid)) {
            second_id = Some(*id);
            second_count = *count;
            second_diag = *diag;
        }
    }


    let second_hits = second_id.and_then(|id| groups.get(&id));
    let mut second_diag_hits: Option<Vec<Hit>> = None;
    if let Some(hits) = second_hits {
        let mut filtered: Vec<Hit> = Vec::new();
        for &hit in hits {
            if hit.diag == second_diag {
                filtered.push(hit);
            }
        }
        if !filtered.is_empty() {
            second_diag_hits = Some(filtered);
        }
    }
    let second_anchor_read_pos = second_hits
        .and_then(|group| anchor_for_diag(group, second_diag))
        .map(|(_, rp)| rp);

    let mut attempt = |hits: &[Hit], diag: i32, anchor_rp: u32| -> Option<AlignmentResult> {
        if hits.is_empty() {
            return None;
        }
        
        if !span_check(hits, cfg) {
            return None;
        }
        
        let ref_id = (hits[0].id_strand >> 1) as i32;
        if (ref_id as usize) >= idx.ref_count() {
            return None;
        }
        let is_rev = (hits[0].id_strand & 1) == 1;
        
        let target_seq = if is_rev {
            rev.clear();
            rev.extend(seq.iter().rev().map(|b| match b {
                b'A' => b'T', b'C' => b'G', b'G' => b'C', b'T' => b'A',
                b'a' => b't', b'c' => b'g', b'g' => b'c', b't' => b'a',
                _ => b'N',
            }));
            rev.as_slice()
        } else {
            seq
        };
        
        let ref_seq = idx.ref_seq(ref_id as usize);
        
        let (cigar, pos, as_score, aligned_len, nm) = greedy_extend(anchor_rp as usize, diag, target_seq, ref_seq, cfg);
        
        if aligned_len == 0 || pos < 0 || (pos as usize) >= ref_seq.len() {
            return None;
        }
        let identity = 1.0 - (nm as f32 / aligned_len.max(1) as f32);
        
        if identity < cfg.min_identity {
            return None;
        }
        Some(AlignmentResult { ref_id, pos, is_rev, mapq: 0, cigar, nm: nm as i32, md: String::new(), as_score, paired: false, proper_pair: false })
    };

    let mut results = Vec::with_capacity(2);

    let mut other_count = second_count;
    if let Some(mut res) = attempt(&best_diag_hits, anchor_diag, anchor_read_pos) {
        let mapq = {
            let frac = (best_count.saturating_sub(other_count)) as f32 / total_seeds.max(1) as f32;
            (frac * 60.0).round().min(60.0) as u8
        };
        res.mapq = mapq;
        results.push(res);
    }

    if let (Some(hits), Some(second_anchor_read_pos)) = (second_diag_hits.as_deref(), second_anchor_read_pos) {
        if let Some(mut res) = attempt(hits, second_diag, second_anchor_read_pos) {
            other_count = best_count;
            let mapq = {
                let frac = (second_count.saturating_sub(other_count)) as f32 / total_seeds.max(1) as f32;
                (frac * 60.0).round().min(60.0) as u8
            };
            res.mapq = mapq;
            results.push(res);
        } else {
            second_failed = true;
        }
    }

    if second_failed && !results.is_empty() {
        let primary = &mut results[0];
        primary.mapq = (primary.mapq / 2).max(1);
    }

    if capped {
        for res in &mut results {
            res.mapq = res.mapq.min(30);
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