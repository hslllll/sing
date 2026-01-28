import sys
import re

TOLERANCE = 10

def parse_dwgsim_qname(qname):
    clean_qname = qname.split('/')[0]
    match = re.match(r'^(.*?)_(\d+)_(\d+)_([01])_([01])_', clean_qname)
    if match:
        ref = match.group(1)
        p1 = int(match.group(2))
        p2 = int(match.group(3))
        return ref, p1, p2
    return None

def inspect(filepath):
    print(f"Inspecting {filepath} for mapping errors (Tolerance: {TOLERANCE}bp)...\n")
    print(f"{'Read Name':<45} | {'True':<12} | {'Mapped':<12} | {'CIGAR':<8} | {'NM':<3} | {'AS':<4} | {'MAPQ':<4} | {'Error'}")
    print("-" * 140)

    count_total = 0
    count_fp = 0
    count_fn = 0
    count_tp = 0
    printed_errors = 0
    MAX_ERRORS = 20

    with open(filepath, 'r') as f:
        for line in f:
            if line.startswith('@'): continue
            parts = line.split('\t')
            if len(parts) < 11: continue
            
            qname = parts[0]
            flag = int(parts[1])
            rname = parts[2]
            try:
                pos = int(parts[3])
                mapq = int(parts[4])
            except ValueError:
                continue
            
            cigar = parts[5]
            nm = "-"
            as_score = "-"
            for field in parts[11:]:
                if field.startswith("NM:i:"):
                    nm = field.split(":")[2]
                elif field.startswith("AS:i:"):
                    as_score = field.split(":")[2]

            if (flag & 256) or (flag & 2048): continue
            
            read_num = 2 if (flag & 128) else 1
            parsed = parse_dwgsim_qname(qname)
            if not parsed: continue

            true_ref, t_p1, t_p2 = parsed
            true_pos = t_p2 if read_num == 2 else t_p1
            
            count_total += 1

            if (flag & 4) or mapq == 0: # Unmapped or low confidence
                count_fn += 1
                continue

            error_type = None
            if rname != true_ref:
                count_fp += 1
                error_type = "Wrong Ref"
            elif abs(pos - true_pos) > TOLERANCE:
                count_fp += 1
                diff = pos - true_pos
                error_type = f"Dist {diff}"
            else:
                count_tp += 1
            
            if error_type and printed_errors < MAX_ERRORS:
                short_qname = qname[:43] + ".." if len(qname) > 45 else qname
                print(f"{short_qname:<45} | {true_ref}:{true_pos:<10} | {rname}:{pos:<10} | {cigar:<8} | {nm:<3} | {as_score:<4} | {mapq:<4} | {error_type}")
                printed_errors += 1

    print("-" * 140)
    print(f"Total: {count_total}, TP: {count_tp}, FP: {count_fp}, FN: {count_fn}")
    # if count_total > 0:
    #     print(f"Accuracy: {count_tp/count_total:.2%}")

if __name__ == "__main__":
    if len(sys.argv) < 2: sys.exit(1)
    inspect(sys.argv[1])
