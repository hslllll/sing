use anyhow::{Context, Result};
use std::borrow::Cow;
use std::io::{BufReader, Read, Write};

const WINDOW: usize = 20;
const SYNC_S: usize = 7;

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

const REMOVE_S: [u64; 4] = [
    rot(BASES[0], (SYNC_S as u32 * ROT) % 64),
    rot(BASES[1], (SYNC_S as u32 * ROT) % 64),
    rot(BASES[2], (SYNC_S as u32 * ROT) % 64),
    rot(BASES[3], (SYNC_S as u32 * ROT) % 64),
];

pub const RADIX: usize = 24;
pub const SHIFT: usize = 32 - RADIX;

#[derive(Clone, Copy)]
pub struct TuningConfig {
    pub chain_max_gap: i32,
    pub seed_weight: u8,
    pub match_score: i32,
    pub mismatch_pen: i32,
    pub gap_open: i32,
    pub gap_ext: i32,
    pub min_identity: f32,
    pub lookup_for_diag: usize,
    pub x_drop: i32,
    pub chain_gap_cost: i32,
    pub chain_min_score: i32,
    pub mapq_scale: f32,
    pub max_mapq: u8,
}

pub const CONFIG: TuningConfig = TuningConfig {
    chain_max_gap: 150,
    seed_weight: 1,
    match_score: 2,
    mismatch_pen: -2,
    gap_open: -3,
    gap_ext: -1,
    min_identity: 0.90,
    lookup_for_diag: 50,
    x_drop: 20,
    chain_gap_cost: 2,
    chain_min_score: 25,
    mapq_scale: 1200.0,
    max_mapq: 60,
};

#[derive(Debug, Clone)]
pub struct Index {
    pub offsets: Vec<u32>,
    pub seeds: Vec<u64>,
    pub freq_filter: u32,
    pub genome_size: u64,
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

    let total_bases: u64 = ref_seqs.iter().map(|s| s.len() as u64).sum();
    let freq_filter = (total_bases as f64).powf(1.0 / 3.0) as usize;
    let freq_filter = std::cmp::max(1, freq_filter);
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
    let mut current_offset = 0u32;
    for count in &global_counts {
        offsets.push(current_offset);
        current_offset += *count;
    }

    Ok(Index { offsets, seeds, freq_filter: freq_filter as u32, genome_size: total_bases, ref_seqs, ref_names })
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
    pub dp_buf: Vec<i32>,
    pub trace_buf: Vec<u8>,
}

impl State {
    pub fn new() -> Self {
        let dp_cap = 1000 * 100;
        Self {
            mins: Vec::with_capacity(100),
            candidates: Vec::with_capacity(1000),
            dp_buf: Vec::with_capacity(dp_cap),
            trace_buf: Vec::with_capacity(dp_cap),
        }
    }
}

#[inline(always)]
fn base_to_index(b: u8) -> Option<usize> {
    match b {
        b'A' | b'a' => Some(0),
        b'C' | b'c' => Some(1),
        b'G' | b'g' => Some(2),
        b'T' | b't' => Some(3),
        _ => None,
    }
}

pub fn get_syncmers(seq: &[u8], out: &mut Vec<(u32, u32)>) {
    out.clear();
    if seq.len() < WINDOW || WINDOW < SYNC_S {
        return;
    }

    const Q_SIZE: usize = 32;
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
                if q_hash[(tail - 1) % Q_SIZE] >= s_hash {
                    tail -= 1;
                } else {
                    break;
                }
            }

            q_hash[tail % Q_SIZE] = s_hash;
            q_pos[tail % Q_SIZE] = s_pos_u32;
            tail += 1;
        }

        if i + 1 >= WINDOW && ambig_k == 0 {
            let k_pos = i + 1 - WINDOW;
            let min_pos = k_pos as u32;
            let max_pos = (k_pos + WINDOW - SYNC_S) as u32;

            while head < tail && q_pos[head % Q_SIZE] < min_pos {
                head += 1;
            }

            if head < tail {
                let m_pos = q_pos[head % Q_SIZE];
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

    candidates.sort_unstable_by(|a, b| {
        a.id_strand.cmp(&b.id_strand)
            .then_with(|| a.ref_pos.cmp(&b.ref_pos))
            .then_with(|| a.read_pos.cmp(&b.read_pos))
    });

    let mut scores = vec![0i32; candidates.len()];

    let mut best_idx = 0;
    let mut max_score = 0;
    
    let group_gap_dist = CONFIG.chain_max_gap;

    let mut group_start = 0;
    while group_start < candidates.len() {
        let id = candidates[group_start].id_strand;
        let mut group_end = group_start;
        while group_end < candidates.len() && candidates[group_end].id_strand == id {
            group_end += 1;
        }

        for i in group_start..group_end {
            let curr = &candidates[i];
            let mut best_s = curr.weight as i32;
            let curr_ref = curr.ref_pos as i64;
            let curr_read = curr.read_pos as i64;
            
            let check_start = if i > CONFIG.lookup_for_diag { i - CONFIG.lookup_for_diag } else { group_start };
            
            for j in (check_start..i).rev() {
                let prev = &candidates[j];
                
                let d_ref = curr_ref - (prev.ref_pos as i64);
                if d_ref <= 0 { continue; }
                if d_ref > group_gap_dist as i64 { break; }

                let d_read = curr_read - (prev.read_pos as i64);
                if d_read <= 0 { continue; } 

                let gap = (d_ref - d_read).abs();
                if gap > CONFIG.chain_max_gap as i64 { continue; }

                let gap_cost = (gap as i32) * CONFIG.chain_gap_cost;
                let score = scores[j] + (curr.weight as i32) - gap_cost;

                if score > best_s {
                    best_s = score;
                }
            }
            scores[i] = best_s;

            if best_s > max_score {
                max_score = best_s;
                best_idx = i;
            }
        }
        group_start = group_end;
    }

    if max_score < CONFIG.chain_min_score {
        return (ChainRes::default(), ChainRes::default()); 
    }

    let target_id = candidates[best_idx].id_strand;
    let center_ref = candidates[best_idx].ref_pos;
    
    let mut start = best_idx;
    let mut end = best_idx;
    
    while start > 0 
        && candidates[start-1].id_strand == target_id 
        && center_ref - candidates[start-1].ref_pos < CONFIG.chain_max_gap as u32
    {
        start -= 1;
    }
    
    while end < candidates.len() - 1 
        && candidates[end+1].id_strand == target_id 
        && candidates[end+1].ref_pos - center_ref < CONFIG.chain_max_gap as u32
    {
        end += 1;
    }

    let best_chain = ChainRes { score: max_score, start, end: end + 1 };
    
    (best_chain, ChainRes::default())
}

fn extend_left(read: &[u8], ref_seq: &[u8], r_start: usize, g_start: usize) -> (usize, usize) {
    let mut mismatches = 0;
    
    let mut score = 0;
    let mut max_score = 0;
    let mut best_len = 0;
    let mut best_mismatches = 0;
    let x_drop = CONFIG.x_drop;

    let len = r_start.min(g_start);
    for i in 1..=len {
        let r_base = read[r_start - i];
        let g_base = ref_seq[g_start - i];
        
        if r_base.eq_ignore_ascii_case(&g_base) {
            score += CONFIG.match_score;
        } else {
            mismatches += 1;
            score += CONFIG.mismatch_pen;
        }

        if score > max_score {
            max_score = score;
            best_len = i;
            best_mismatches = mismatches;
        }

        if max_score - score > x_drop {
            break;
        }
    }
    (best_len, best_mismatches)
}

fn extend_right(read: &[u8], ref_seq: &[u8], r_start: usize, g_start: usize) -> (usize, usize) {
    let mut mismatches = 0;
    
    let mut score = 0;
    let mut max_score = 0;
    let mut best_len = 0;
    let mut best_mismatches = 0;
    let x_drop = CONFIG.x_drop;
    
    let r_len = read.len();
    let g_len = ref_seq.len();
    let mut i = 0;
    
    while r_start + i < r_len && g_start + i < g_len {
        let r_base = read[r_start + i];
        let g_base = ref_seq[g_start + i];
        
        if r_base.eq_ignore_ascii_case(&g_base) {
            score += CONFIG.match_score;
        } else {
            mismatches += 1;
            score += CONFIG.mismatch_pen;
        }

        let curr_len = i + 1;
        if score > max_score {
            max_score = score;
            best_len = curr_len;
            best_mismatches = mismatches;
        }

        if max_score - score > x_drop {
            break;
        }
        i += 1;
    }
    (best_len, best_mismatches)
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


fn align_gap_nw(
    cigar: &mut Vec<u32>,
    read: &[u8],
    ref_seq: &[u8],
    dp_buf: &mut Vec<i32>,
    trace_buf: &mut Vec<u8>,
) {
    let n = read.len();
    let m = ref_seq.len();

    
    
    if n == 0 || m == 0 {
        if n > 0 { push_cigar(cigar, n as u32, 1); } 
        if m > 0 { push_cigar(cigar, m as u32, 2); } 
        return;
    }

    
    
    let match_score = CONFIG.match_score;
    let mismatch_score = CONFIG.mismatch_pen;
    let gap_score = CONFIG.gap_open; 
    
    let cols = m + 1;
    let needed = (n + 1) * cols;
    if dp_buf.len() < needed {
        dp_buf.resize(needed, 0);
    }
    if trace_buf.len() < needed {
        trace_buf.resize(needed, 0);
    }
    let dp = &mut dp_buf[..needed];
    let trace = &mut trace_buf[..needed];
    dp[0] = 0;
    trace[0] = 0;

    for i in 1..=n {
        let idx = i * cols;
        dp[idx] = i as i32 * gap_score;
        trace[idx] = 2;
    }
    for j in 1..=m {
        dp[j] = j as i32 * gap_score;
        trace[j] = 3;
    }

    for i in 1..=n {
        let curr_row = i * cols;
        let prev_row = (i - 1) * cols;
        
        for j in 1..=m {
            let is_match = read[i - 1].eq_ignore_ascii_case(&ref_seq[j - 1]);
            let s_char = if is_match { match_score } else { mismatch_score };

            let score_diag = dp[prev_row + j - 1] + s_char;
            let score_up = dp[prev_row + j] + gap_score;   
            let score_left = dp[curr_row + j - 1] + gap_score; 

            if score_diag >= score_up && score_diag >= score_left {
                dp[curr_row + j] = score_diag;
                trace[curr_row + j] = 1;
            } else if score_up >= score_left {
                dp[curr_row + j] = score_up;
                trace[curr_row + j] = 2;
            } else {
                dp[curr_row + j] = score_left;
                trace[curr_row + j] = 3;
            }
        }
    }

    
    let mut i = n;
    let mut j = m;
    let mut ops_len = 0usize;

    while i > 0 || j > 0 {
        match trace[i * cols + j] {
            1 => { 
                dp[ops_len] = 0;
                ops_len += 1;
                i -= 1;
                j -= 1;
            }
            2 => { 
                dp[ops_len] = 1;
                ops_len += 1;
                i -= 1;
            }
            3 => { 
                dp[ops_len] = 2;
                ops_len += 1;
                j -= 1;
            }
            _ => break, 
        }
    }

    if ops_len == 0 {
        return;
    }

    let mut curr_op = dp[ops_len - 1] as u32;
    let mut curr_len = 0u32;
    for idx in (0..ops_len).rev() {
        let op = dp[idx] as u32;
        if op == curr_op {
            curr_len += 1;
        } else {
            push_cigar(cigar, curr_len, curr_op);
            curr_op = op;
            curr_len = 1;
        }
    }
    push_cigar(cigar, curr_len, curr_op);
}

fn gen_cigar(
    hits: &[Hit],
    read_seq: &[u8],
    ref_seq: &[u8],
    dp_buf: &mut Vec<i32>,
    trace_buf: &mut Vec<u8>,
) -> (Vec<u32>, i32) {
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
    let (match_len_left, _) = extend_left(read_seq, ref_seq, r_start, g_start);
    let leading_clip = r_start - match_len_left;
    let final_ref_start = (g_start - match_len_left) as i32;

    if leading_clip > 0 {
        push_cigar(&mut cigar, leading_clip as u32, 4);
    }
    if match_len_left > 0 {
        push_cigar(&mut cigar, match_len_left as u32, 0);
    }

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
            let g_gap_seq = if curr_g < g_gap_end { &ref_seq[curr_g..g_gap_end] } else { &[] };

            align_gap_nw(&mut cigar, r_gap_seq, g_gap_seq, dp_buf, trace_buf);
        }

        push_cigar(&mut cigar, WINDOW as u32, 0);
        curr_r = next_r + WINDOW;
        curr_g = next_g + WINDOW;
    }
    let r_end_idx = curr_r;
    let g_end_idx = curr_g;
    let (match_len_right, _) = extend_right(read_seq, ref_seq, r_end_idx, g_end_idx);

    if match_len_right > 0 {
        push_cigar(&mut cigar, match_len_right as u32, 0);
    }

    let consumed_r = r_end_idx + match_len_right;
    let trailing_clip = read_seq.len().saturating_sub(consumed_r);
    if trailing_clip > 0 {
        push_cigar(&mut cigar, trailing_clip as u32, 4);
    }

    (cigar, final_ref_start)
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
    let State { mins, candidates, dp_buf, trace_buf } = state;

    candidates.clear();

    get_syncmers(seq, mins);
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
    collect_candidates(idx, mins, true, candidates);

    if candidates.is_empty() {
        return None;
    }

    let (c1, c2) = find_top_chains(candidates);
    let evaluate = |chain: &ChainRes,
                    seq: &[u8],
                    rev: &[u8],
                    candidates: &[Hit],
                    dp_buf: &mut Vec<i32>,
                    trace_buf: &mut Vec<u8>|
     -> Option<(AlignmentResult, i32)> {
        if chain.score == 0 {
            return None;
        }
        let hits = &candidates[chain.start..chain.end];
        let ref_id = (hits[0].id_strand >> 1) as i32;
        let is_rev = (hits[0].id_strand & 1) == 1;
        let target_seq = if is_rev { rev } else { seq };
        let ref_seq = &idx.ref_seqs[ref_id as usize];

        let (cigar, pos) = gen_cigar(hits, target_seq, ref_seq, dp_buf, trace_buf);
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

    let r1 = evaluate(&c1, seq, rev, candidates, dp_buf, trace_buf);
    let r2 = if c2.score > 0 { evaluate(&c2, seq, rev, candidates, dp_buf, trace_buf) } else { None };
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

    let diff = if second_chain_score > 0 {
        (best_chain_score - second_chain_score) as f32
    } else {
        best_chain_score as f32
    };
    let read_len_f32 = seq.len() as f32;
    let evidence_weight = (best_chain_score as f32).min(read_len_f32) / read_len_f32;
    let raw_mapq = CONFIG.mapq_scale * diff * evidence_weight;
    best_res.mapq = raw_mapq.min(CONFIG.max_mapq as f32) as u8;

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
