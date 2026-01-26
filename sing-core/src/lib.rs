use anyhow::{Context, Result};
use block_aligner::scan_block::*;
use block_aligner::scores::*;
use std::borrow::Cow;
use std::collections::VecDeque;
use std::io::{BufReader, Read, Write};

pub const WINDOW: usize = 19;  
pub const MIN_W: usize = WINDOW;
pub const BATCH_SIZE: usize = 10_000;
pub const BATCH_CAP: usize = 64;
pub const FILTERFOR_100BP: f64 = 10_000_000.0;
pub const AMBIGUOUS_SCORE: f32 = 1.001;

const BASES: [u64; 4] = [
    0x243f_6a88_85a3_08d3,
    0xb7e1_5162_8aed_2a6a,
    0x6a09_e667_f3bc_c908,
    0xbb67_ae85_84ca_a73b,
];

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

const SYNC_S: usize = 11;  
const SYNC_WINDOW: usize = WINDOW - SYNC_S + 1; 
const REMOVE_S: [u64; 4] = [
    rot(BASES[0], (SYNC_S as u32 * ROT) % 64),
    rot(BASES[1], (SYNC_S as u32 * ROT) % 64),
    rot(BASES[2], (SYNC_S as u32 * ROT) % 64),
    rot(BASES[3], (SYNC_S as u32 * ROT) % 64),
];

pub const RADIX: usize = 24;
pub const SHIFT: usize = 32 - RADIX;


const BASE_TO_IDX: [usize; 256] = {
    let mut table = [4usize; 256];
    table[b'A' as usize] = 0;
    table[b'a' as usize] = 0;
    table[b'C' as usize] = 1;
    table[b'c' as usize] = 1;
    table[b'G' as usize] = 2;
    table[b'g' as usize] = 2;
    table[b'T' as usize] = 3;
    table[b't' as usize] = 3;
    table
};

#[derive(Clone, Copy)]
pub struct TuningConfig {
    pub chain_max_gap: i32,
    pub seed_weight: u8,
    pub match_score: i32,
    pub mismatch_pen: i32,
    pub gap_open: i32,
    pub gap_ext: i32,
    pub min_identity: f32,
}

pub const CONFIG: TuningConfig = TuningConfig {
    chain_max_gap: 500,
    seed_weight: 1,
    match_score: 2,
    mismatch_pen: -2,
    gap_open: -3,
    gap_ext: -1,
    min_identity: 0.5,
};

#[derive(Debug, Clone)]
pub struct Index {
    pub offsets: Vec<u32>,
    pub seeds: Vec<u64>,
    pub seed_boundaries: Vec<u32>,
    pub freq_filter: u32,
    pub ref_seqs: Vec<Vec<u8>>,
    pub ref_names: Vec<String>,
}

impl Index {
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

        r.read_exact(&mut buf8).context("read seed_boundaries len")?;
        let len = u64::from_le_bytes(buf8) as usize;
        let mut seed_boundaries = Vec::with_capacity(len);
        for _ in 0..len {
            r.read_exact(&mut buf4).context("read seed_boundary")?;
            seed_boundaries.push(u32::from_le_bytes(buf4));
        }

        r.read_exact(&mut buf4).context("read freq_filter")?;
        let freq_filter = u32::from_le_bytes(buf4);

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

        Ok(Index { offsets, seeds, seed_boundaries, freq_filter, ref_seqs, ref_names })
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

        w.write_all(&(self.seed_boundaries.len() as u64).to_le_bytes())?;
        for b in &self.seed_boundaries {
            w.write_all(&b.to_le_bytes())?;
        }

        w.write_all(&self.freq_filter.to_le_bytes())?;

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
    let mut seeds_tmp: Vec<(u64, u64)> = Vec::new();

    for (rid, (name, seq)) in records.into_iter().enumerate() {
        ref_names.push(name.into());

        let mut mins = Vec::new();
        get_syncmers(&seq, &mut mins);
        for (h, p) in mins {
            let hash64 = h as u64;
            let hash16 = (hash64 & 0xFFFF) as u64;
            let packed = (hash16 << 48) | ((rid as u64) << 32) | (p as u64);
            seeds_tmp.push((hash64, packed));
        }
        ref_seqs.push(seq);
    }

    let total_bases: usize = ref_seqs.iter().map(|s| s.len()).sum();
    let freq_filter = ((total_bases as f64 / FILTERFOR_100BP) * 100.0) as usize;
    let freq_filter = std::cmp::max(100, freq_filter);
    eprintln!("Dynamic freq_filter set to: {} (Genome size: {} bp)", freq_filter, total_bases);

    seeds_tmp.sort_unstable_by_key(|k| k.0);

    let mut global_counts = vec![0u32; 1 << RADIX];
    let mut seeds: Vec<u64> = Vec::new();
    let mut i = 0usize;
    let mut uniq_hashes = 0;
    let mut kept_hashes = 0;

    while i < seeds_tmp.len() {
        let mut j = i;
        while j < seeds_tmp.len() && seeds_tmp[j].0 == seeds_tmp[i].0 {
            j += 1;
        }
        uniq_hashes += 1;
        if j - i <= freq_filter {
            kept_hashes += 1;
            for k in i..j {
                let (h, packed) = seeds_tmp[k];
                let bucket = (h >> SHIFT) as usize;
                global_counts[bucket] += 1;
                seeds.push(packed);
            }
        }
        i = j;
    }

    eprintln!("Indexing statistics:");
    eprintln!("  Total hashes: {}", seeds_tmp.len());
    eprintln!("  Unique hashes: {}", uniq_hashes);
    eprintln!("  Unique hashes kept: {}", kept_hashes);
    eprintln!("  Total seeds kept: {}", seeds.len());

    let mut offsets = Vec::with_capacity(global_counts.len());
    let mut seed_boundaries = Vec::with_capacity(global_counts.len());
    let mut current_offset = 0u32;
    
    for (_bucket, count) in global_counts.iter().enumerate() {
        offsets.push(current_offset);
        current_offset += *count;
    }

    let mut prev_hash = 0u64;
    for (i, &seed) in seeds.iter().enumerate() {
        let full_hash = (seed >> 48) | ((((seed >> 32) & 0xFFFF) as u64) << 16);
        if full_hash != prev_hash {
            seed_boundaries.push(i as u32);
            prev_hash = full_hash;
        }
    }
    seed_boundaries.push(seeds.len() as u32);

    Ok(Index { offsets, seeds, seed_boundaries, freq_filter: freq_filter as u32, ref_seqs, ref_names })
}

#[derive(Clone, Copy, Default, Debug)]
pub struct Hit {
    pub id_strand: u32,
    pub diag: i32,
    pub read_pos: u32,
    pub ref_pos: u32,
    pub weight: u8,
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
        Self { mins: Vec::with_capacity(100), candidates: Vec::with_capacity(1000) }
    }
}

#[inline(always)]
fn base_to_index(b: u8) -> Option<usize> {
    let idx = BASE_TO_IDX[b as usize];
    if idx < 4 {
        Some(idx)
    } else {
        None
    }
}

pub fn get_syncmers(seq: &[u8], out: &mut Vec<(u32, u32)>) {
    out.clear();
    if seq.len() < WINDOW {
        return;
    }

    
    let mut h_k = 0u64;
    let mut h_s = 0u64;
    let mut base_buf_k = [0usize; WINDOW];
    let mut base_buf_s = [0usize; SYNC_S];
    let mut ambig_buf_k = [0u8; WINDOW];
    let mut ambig_buf_s = [0u8; SYNC_S];
    let mut ambig_k = 0i32;
    let mut ambig_s = 0i32;

    
    let mut minq: VecDeque<(u64, usize)> = VecDeque::with_capacity(SYNC_WINDOW + 1);
    let mut s_hash_buf = [0u64; SYNC_WINDOW];
    let mut s_valid = [false; SYNC_WINDOW];

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

        
        if i + 1 >= SYNC_S {
            let s_start = i + 1 - SYNC_S; 
            let slot = s_start % SYNC_WINDOW;

            if ambig_s == 0 {
                let s_hash = h_s.wrapping_mul(0x517cc1b727220a95);
                s_hash_buf[slot] = s_hash;
                s_valid[slot] = true;

                while let Some(&(h, _)) = minq.back() {
                    if h >= s_hash {
                        minq.pop_back();
                    } else {
                        break;
                    }
                }
                minq.push_back((s_hash, s_start));
            } else {
                s_valid[slot] = false;
            }
        }

        
        if i + 1 >= WINDOW && ambig_k == 0 {
            let k_pos = i + 1 - WINDOW;

            
            while let Some(&(_, pos)) = minq.front() {
                if pos < k_pos {
                    minq.pop_front();
                } else {
                    break;
                }
            }

            let k_slot = k_pos % SYNC_WINDOW;
            if s_valid[k_slot] {
                if let Some(&(min_hash, min_pos)) = minq.front() {
                    let s_hash_at_k = s_hash_buf[k_slot];
                    if min_pos == k_pos && min_hash == s_hash_at_k {
                        let k_hash = h_k.wrapping_mul(0x517cc1b727220a95);
                        if out.last().map(|&(_, p)| p) != Some(k_pos as u32) {
                            out.push((k_hash as u32, k_pos as u32));
                        }
                    }
                }
            }
        }
    }
}

fn collect_candidates(idx: &Index, mins: &[(u32, u32)], is_rev: bool, out: &mut Vec<Hit>) {
    let freq_filter = idx.freq_filter as usize;
    
    for &(h, r_pos) in mins {
        let bucket = (h >> SHIFT) as usize;
        let start = idx.offsets[bucket] as usize;
        let end = idx
            .offsets
            .get(bucket + 1)
            .copied()
            .unwrap_or(idx.seeds.len() as u32) as usize;

        let count = end - start;
        if count > freq_filter {
            continue;
        }

        let weight = CONFIG.seed_weight;
        let target_hash = (h & 0xFFFF) as u64;

        for i in start..end {
            let seed = idx.seeds[i];
            if (seed >> 48) == target_hash {
                let rid = ((seed >> 32) & 0xFFFF) as u32;
                let pos = seed as u32;

                let id_strand = (rid << 1) | (is_rev as u32);
                let diag = (pos as i32) - (r_pos as i32);

                out.push(Hit { id_strand, diag, read_pos: r_pos, ref_pos: pos, weight });
                
            }
        }
    }
}

fn find_top_chains(candidates: &mut Vec<Hit>) -> (ChainRes, ChainRes) {
    if candidates.is_empty() {
        return (ChainRes::default(), ChainRes::default());
    }

    candidates.sort_unstable_by_key(|h| (h.id_strand, h.read_pos));

    let mut best = ChainRes::default();
    let mut second = ChainRes::default();

    let mut group_start = 0;
    while group_start < candidates.len() {
        let id = candidates[group_start].id_strand;
        let mut group_end = group_start;
        while group_end < candidates.len() && candidates[group_end].id_strand == id {
            group_end += 1;
        }

        
        let n = group_end - group_start;
        let mut dp = vec![0i32; n];
        let mut backtrack = vec![None; n];
        
        for i in 0..n {
            let hit_i = &candidates[group_start + i];
            dp[i] = hit_i.weight as i32;
            
            
            for j in 0..i {
                let hit_j = &candidates[group_start + j];
                
                
                if hit_j.read_pos < hit_i.read_pos && hit_j.ref_pos < hit_i.ref_pos {
                    let gap_read = (hit_i.read_pos - hit_j.read_pos) as i32;
                    let gap_ref = (hit_i.ref_pos - hit_j.ref_pos) as i32;
                    let gap_diff = (gap_read - gap_ref).abs();
                    
                    
                    let gap_penalty = if gap_diff > CONFIG.chain_max_gap {
                        i32::MIN / 2  
                    } else {
                        (gap_diff * CONFIG.gap_ext.abs()) / 10  
                    };
                    
                    let score = dp[j] + hit_i.weight as i32 - gap_penalty;
                    if score > dp[i] {
                        dp[i] = score;
                        backtrack[i] = Some(j);
                    }
                }
            }
        }
        
        
        let mut best_local_idx = 0;
        let mut best_local_score = dp[0];
        for i in 1..n {
            if dp[i] > best_local_score {
                best_local_score = dp[i];
                best_local_idx = i;
            }
        }
        
        
        let mut chain_indices = Vec::new();
        let mut curr = Some(best_local_idx);
        while let Some(idx) = curr {
            chain_indices.push(idx);
            curr = backtrack[idx];
        }
        chain_indices.reverse();
        
        if !chain_indices.is_empty() {
            let chain_start = group_start + chain_indices[0];
            let chain_end = group_start + best_local_idx + 1;
            let chain_score = best_local_score;
            
            if chain_score > best.score {
                second = best;
                best = ChainRes { score: chain_score, start: chain_start, end: chain_end };
            } else if chain_score > second.score {
                second = ChainRes { score: chain_score, start: chain_start, end: chain_end };
            }
        }

        group_start = group_end;
    }

    (best, second)
}

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


fn align_with_block_aligner(read: &[u8], ref_seq: &[u8], cigar: &mut Vec<u32>) {
    if read.is_empty() || ref_seq.is_empty() {
        if !read.is_empty() {
            push_cigar(cigar, read.len() as u32, 1);
        }
        if !ref_seq.is_empty() {
            push_cigar(cigar, ref_seq.len() as u32, 2);
        }
        return;
    }

    
    
    let gaps = Gaps { open: CONFIG.gap_open as i8, extend: CONFIG.gap_ext as i8 };
    
    let scoring = NucMatrix::new_simple(CONFIG.match_score as i8, CONFIG.mismatch_pen as i8);

    
    const BA_MAX: usize = 256;
    let q_len = read.len().min(BA_MAX);
    let r_len = ref_seq.len().min(BA_MAX);
    let query = PaddedBytes::from_bytes::<NucMatrix>(&read[..q_len], 512);
    let reference = PaddedBytes::from_bytes::<NucMatrix>(&ref_seq[..r_len], 512);

    
    let mut block_aligner = Block::<true, false>::new(q_len, r_len, 128);

    block_aligner.align(&query, &reference, &scoring, gaps, 0..=0, 0);
    
    let res = block_aligner.res();
    
    
    
    
    let q_end = res.query_idx;
    let r_end = res.reference_idx;
    
    
    let aligned_len = q_end.min(r_end);
    if aligned_len > 0 {
        push_cigar(cigar, aligned_len as u32, 0);
    }
    
    if q_end > r_end {
        push_cigar(cigar, (q_end - r_end) as u32, 1);
    } else if r_end > q_end {
        push_cigar(cigar, (r_end - q_end) as u32, 2);
    }
}

fn gen_cigar(hits: &[Hit], read_seq: &[u8], ref_seq: &[u8]) -> (Vec<u32>, i32) {
    if hits.is_empty() {
        return (vec![], 0);
    }

    let mut sorted = hits.to_vec();
    sorted.sort_unstable_by_key(|h| h.read_pos);

    let mut filtered = Vec::with_capacity(sorted.len());
    filtered.push(sorted[0]);
    for i in 1..sorted.len() {
        let last = filtered.last().unwrap();
        let curr = sorted[i];
        if curr.ref_pos > last.ref_pos && curr.read_pos > last.read_pos {
            filtered.push(curr);
        }
    }

    let mut cigar = Vec::with_capacity(32);
    let first = filtered[0];


    let r_start = first.read_pos as usize;
    let g_start = first.ref_pos as usize;
    
        
    let left_extend_len = r_start.min(g_start).min(100);
        let final_ref_start = if left_extend_len > 0 && r_start >= left_extend_len && g_start >= left_extend_len {
            let r_left = &read_seq[r_start - left_extend_len..r_start];
            let g_left = &ref_seq[g_start - left_extend_len..g_start];
        
            let left_clip = r_start - left_extend_len;
        if left_clip > 0 {
            push_cigar(&mut cigar, left_clip as u32, 4);
        }
        
            align_with_block_aligner(r_left, g_left, &mut cigar);
            g_start - left_extend_len
    } else {
        if r_start > 0 {
            push_cigar(&mut cigar, r_start as u32, 4);
        }
        g_start
        };

    
    let mut curr_r = r_start;
    let mut curr_g = g_start;
    push_cigar(&mut cigar, WINDOW as u32, 0);
    curr_r += WINDOW;
    curr_g += WINDOW;

    for i in 1..filtered.len() {
        let hit = filtered[i];
        let next_r = hit.read_pos as usize;
        let next_g = hit.ref_pos as usize;

        if next_r < curr_r || next_g < curr_g {
            continue;
        }

        let gap_r = next_r - curr_r;
        let gap_g = next_g - curr_g;

        if gap_r > 0 || gap_g > 0 {
            let r_gap_seq = &read_seq[curr_r..next_r];
            let g_gap_end = next_g.min(ref_seq.len());
            let g_gap_seq = if curr_g < g_gap_end { 
                &ref_seq[curr_g..g_gap_end] 
            } else { 
                &[] 
            };

            align_with_block_aligner(r_gap_seq, g_gap_seq, &mut cigar);
        }

        push_cigar(&mut cigar, WINDOW as u32, 0);
        curr_r = next_r + WINDOW;
        curr_g = next_g + WINDOW;
    }
    
        
    let r_end_idx = curr_r;
    let g_end_idx = curr_g;
    let right_extend_len = (read_seq.len() - r_end_idx).min(ref_seq.len() - g_end_idx).min(100);
    
    if right_extend_len > 0 {
        let r_right = &read_seq[r_end_idx..r_end_idx + right_extend_len];
        let g_right = &ref_seq[g_end_idx..g_end_idx + right_extend_len];
        
            align_with_block_aligner(r_right, g_right, &mut cigar);
    }

    let consumed_r = r_end_idx + right_extend_len;
    let trailing_clip = read_seq.len().saturating_sub(consumed_r);
    if trailing_clip > 0 {
        push_cigar(&mut cigar, trailing_clip as u32, 4);
    }

    (cigar, final_ref_start as i32)
}

fn compute_metrics(read_seq: &[u8], ref_seq: &[u8], start_pos: i32, cigar: &[u32]) -> (i32, String, i32) {
    let mut nm = 0;
    let mut as_score = 0;

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
                            as_score += CONFIG.match_score;
                        } else {
                            nm += 1;
                            as_score += CONFIG.mismatch_pen;
                        }
                    } else {
                        nm += 1;
                        as_score += CONFIG.mismatch_pen;
                    }
                    r_idx += 1;
                    q_idx += 1;
                }
            }
            8 => {
                nm += len as i32;
                as_score += len as i32 * CONFIG.mismatch_pen;
                r_idx += len;
                q_idx += len;
            }
            1 => {
                nm += len as i32;
                as_score += CONFIG.gap_open + ((len as i32 - 1) * CONFIG.gap_ext);
                q_idx += len;
            }
            2 => {
                nm += len as i32;
                as_score += CONFIG.gap_open + ((len as i32 - 1) * CONFIG.gap_ext);
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

pub fn align(seq: &[u8], idx: &Index, state: &mut State, rev: &mut Vec<u8>) -> Option<AlignmentResult> {
    state.candidates.clear();

    get_syncmers(seq, &mut state.mins);
    collect_candidates(idx, &state.mins, false, &mut state.candidates);

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
    get_syncmers(rev, &mut state.mins);
    collect_candidates(idx, &state.mins, true, &mut state.candidates);

    if state.candidates.is_empty() {
        return None;
    }

    let (c1, c2) = find_top_chains(&mut state.candidates);
    let evaluate = |chain: &ChainRes, seq: &[u8], rev: &[u8]| -> Option<(AlignmentResult, i32)> {
        if chain.score == 0 {
            return None;
        }
        let hits = &state.candidates[chain.start..chain.end];
        let ref_id = (hits[0].id_strand >> 1) as i32;
        let is_rev = (hits[0].id_strand & 1) == 1;
        let target_seq = if is_rev { rev } else { seq };
        let ref_seq = &idx.ref_seqs[ref_id as usize];

        let (cigar, pos) = gen_cigar(hits, target_seq, ref_seq);
        if pos < 0 {
            return None;
        }

        let (nm, md, as_score) = compute_metrics(target_seq, ref_seq, pos, &cigar);
        let identity = 1.0 - (nm as f32 / seq.len() as f32);

        if identity < CONFIG.min_identity {
            return None;
        }
        Some((
            AlignmentResult { ref_id, pos, is_rev, mapq: 0, cigar, nm, md, as_score },
            chain.score,
        ))
    };

    let r1 = evaluate(&c1, seq, rev);
    let r2 = if c2.score > 0 { evaluate(&c2, seq, rev) } else { None };
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

    if second_chain_score > 0 {
        let score_ratio = best_chain_score as f32 / second_chain_score as f32;
        if score_ratio < AMBIGUOUS_SCORE {
            return None;
        }
    }

    let diff = if second_chain_score > 0 {
        (best_chain_score - second_chain_score) as f32
    } else {
        best_chain_score as f32
    };

    let read_len_f32 = seq.len() as f32;
    let evidence_weight = (best_chain_score as f32).min(read_len_f32) / read_len_f32;
    let raw_mapq = 1200.0 * diff * evidence_weight;
    best_res.mapq = raw_mapq.min(60.0) as u8;

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
