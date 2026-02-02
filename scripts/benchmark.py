#!/usr/bin/env python3
import argparse
import os
import re
import shutil
import subprocess
import sys
import time
from pathlib import Path

TOLERANCE = 10

def run(cmd, cwd=None, env=None, quiet=False):
    if not quiet:
        print("+", " ".join(cmd))
    subprocess.run(cmd, cwd=cwd, env=env, check=True)


def run_shell(cmd, cwd=None):
    print("+", cmd)
    subprocess.run(cmd, cwd=cwd, shell=True, check=True, executable="/bin/bash")


def ensure_dir(path: Path):
    path.mkdir(parents=True, exist_ok=True)


def get_mode_config(mode: str):
    if mode == "y":
        return {
            "species": "yeast",
            "genome_size": 12157105,
            "ref_name": "Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa",
            "ref_url": "http://ftp.ensembl.org/pub/release-110/fasta/saccharomyces_cerevisiae/dna/Saccharomyces_cerevisiae.R64-1-1.dna.toplevel.fa.gz",
            "index_prefix": "yeast",
        }
    if mode == "a":
        return {
            "species": "arabidopsis",
            "genome_size": 135000000,
            "ref_name": "Arabidopsis_thaliana.TAIR10.dna.toplevel.fa",
            "ref_url": "http://ftp.ensemblgenomes.org/pub/plants/release-57/fasta/arabidopsis_thaliana/dna/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa.gz",
            "index_prefix": "arabidopsis",
        }
    if mode == "m":
        return {
            "species": "maize",
            "genome_size": 250000000,
            "ref_name": "Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa",
            "ref_url": "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/current/fasta/zea_mays/dna/Zea_mays.Zm-B73-REFERENCE-NAM-5.0.dna.toplevel.fa.gz",
            "index_prefix": "maize",
        }
    if mode == "h":
        return {
            "species": "human",
            "genome_size": 3.2e9,
            "ref_name": "GRCh38_latest_genomic.fna",
            "ref_url": "https://ftp.ncbi.nlm.nih.gov/refseq/H_sapiens/annotation/GRCh38_latest/refseq_identifiers/GRCh38_latest_genomic.fna.gz",
            "index_prefix": "human",
        }
    if mode == "b":
        return {
            "species": "brassica",
            "genome_size": 400000000,
            "ref_name": "Brassica_rapa.Brapa_1.0.dna.toplevel.fa",
            "ref_url": "https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-57/fasta/brassica_rapa/dna/Brassica_rapa.Brapa_1.0.dna.toplevel.fa.gz",
            "index_prefix": "brassica",
        }
    raise SystemExit("Unknown mode. Use one of: y/a/m/h/b")


def download_and_decompress(ref_gz: Path, ref_url: str, ref_decomp: Path):
    if not ref_gz.exists():
        run(["wget", "-c", "-O", str(ref_gz), ref_url])
    if not ref_decomp.exists():
        run_shell(f"pigz -p 8 -cd '{ref_gz}' > '{ref_decomp}'")


def ensure_sing_binary():
    sing_bin = Path("./target/release/sing")
    if not sing_bin.exists():
        run(["cargo", "build", "--release"])
    return sing_bin

def ensure_sing_index(ref_decomp: Path, index_path: Path):
    if not index_path.exists():
        run(["cargo", "run", "--release", "--", "build", str(ref_decomp), str(index_path)])


def ensure_bwa_index(ref_decomp: Path):
    if not Path(str(ref_decomp) + ".0123").exists():
        run(["bwa-mem2", "index", str(ref_decomp)])


def ensure_bowtie2_index(ref_decomp: Path):
    if not Path(str(ref_decomp) + ".1.bt2").exists():
        run(["bowtie2-build", "-t", "8", str(ref_decomp), str(ref_decomp)])


def ensure_filtered_ref(ref_decomp: Path) -> Path:
    filtered = Path(str(ref_decomp).replace(".fa", ".filtered.fa"))
    if not filtered.exists():
        awk_script = r"""awk '/^>/ {header=$0; lc=tolower(header); skip=0; if (lc ~ /(^|[^a-z])(ct|pt|mt)([^a-z]|$)/) skip=1; if (lc ~ /scaffold/ || lc ~ /contig/) skip=1; if (!skip) print header; next} !skip { print }'"""
        run_shell(f"{awk_script} '{ref_decomp}' > '{filtered}'")
    return filtered


def run_mapper_to_sorted_bam(cmd, out_bam: Path, threads: int):
    # Resolve absolute path for time binary
    time_bin = Path("./time").resolve()
    if not time_bin.exists():
        time_bin = Path("/usr/bin/time")
    
    if not time_bin.exists():
        print(f"Error: {time_bin} not found. installation required.")
        return None, None

    sam_path = out_bam.with_suffix(".sam")
    unsorted = out_bam.with_suffix(".unsorted.bam")
    time_log = out_bam.with_suffix(".time.log")
    elapsed = None
    mem_kb = None
    
    mapper_cmd = [str(time_bin), "-f", "%e %M", "-o", str(time_log)] + cmd
    print("+", " ".join(mapper_cmd), ">", str(sam_path))
    
    try:
        with sam_path.open("w") as f_out:
            subprocess.run(mapper_cmd, stdout=f_out, check=True)
    except Exception as e:
        print(f"Mapper execution error: {e}")
        return None, None

    if time_log.exists():
        log_parts = time_log.read_text().strip().split()
        if len(log_parts) == 2:
            elapsed = float(log_parts[0])
            try:
                mem_kb = int(float(log_parts[1]))
            except ValueError:
                mem_kb = None
    
    # SAM -> BAM
    print("+", f"samtools view -b -o {unsorted} {sam_path}")
    subprocess.run(["samtools", "view", "-b", "-o", str(unsorted), str(sam_path)], check=True)
    sam_path.unlink(missing_ok=True)
    
    sort_threads = max(1, min(threads, 8))
    sort_mem = os.environ.get("SAMTOOLS_SORT_MEM", "512M")
    tmp_prefix = str(out_bam.with_suffix(""))
    run([
        "samtools",
        "sort",
        "-@",
        str(sort_threads),
        "-m",
        sort_mem,
        "-T",
        tmp_prefix,
        "-O",
        "BAM",
        "-o",
        str(out_bam),
        str(unsorted),
    ])
    run(["samtools", "index", "-@", str(sort_threads), str(out_bam)])
    unsorted.unlink(missing_ok=True)
    time_log.unlink(missing_ok=True)
    return elapsed, mem_kb


def parse_dwgsim_qname(qname):
    clean = qname.split("/")[0]
    m = re.match(r"^(.*?)_(\d+)_(\d+)_([01])_([01])_", clean)
    if not m:
        return None
    return clean, m.group(1), int(m.group(2)), int(m.group(3))


def load_bams_and_build_truth(filepaths):
    truth_dict = {}
    tool_results = {name: {} for name in filepaths}

    for tool_name, filepath in filepaths.items():
        if not Path(filepath).exists():
            continue
        proc = subprocess.Popen(["samtools", "view", "-h", str(filepath)], stdout=subprocess.PIPE, text=True)
        for line in proc.stdout:
            if line.startswith("@"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            qname = parts[0]
            flag = int(parts[1])
            rname = parts[2]
            pos = int(parts[3])
            mapq = int(parts[4])

            read_num = 2 if (flag & 128) else 1

            parsed = parse_dwgsim_qname(qname)
            if parsed:
                q_base, true_ref, t_p1, t_p2 = parsed
                if (q_base, 1) not in truth_dict:
                    truth_dict[(q_base, 1)] = (true_ref, t_p1)
                if (q_base, 2) not in truth_dict:
                    truth_dict[(q_base, 2)] = (true_ref, t_p2)

            if (flag & 256) or (flag & 2048):
                continue
            q_base = qname.split("/")[0]

            if (flag & 4) or mapq == 0 or rname == "*":
                tool_results[tool_name][(q_base, read_num)] = None
            else:
                tool_results[tool_name][(q_base, read_num)] = (rname, pos)
        proc.wait()

    return truth_dict, tool_results


def calculate_metrics(truth_dict, tool_map):
    tp = 0
    fp = 0
    fn = 0

    for key, (true_ref, true_pos) in truth_dict.items():
        mapped_val = tool_map.get(key)
        if mapped_val is None:
            fn += 1
        else:
            pred_ref, pred_pos = mapped_val
            if pred_ref == true_ref and abs(pred_pos - true_pos) <= TOLERANCE:
                tp += 1
            else:
                fp += 1

    total_positives = len(truth_dict)
    precision = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall_fraction = tp / total_positives if total_positives > 0 else 0.0
    recall = recall_fraction * 100
    f1 = 2 * (precision * recall_fraction) / (precision + recall_fraction) if (precision + recall_fraction) > 0 else 0.0
    return tp, fp, fn, precision * 100, recall, f1 * 100


def ensure_tools(*tools):
    for t in tools:
        if shutil.which(t) is None:
            raise SystemExit(f"Error: {t} not found. Please install it.")


def benchmark_all(mode, threads_override):
    cfg = get_mode_config(mode)
    ref_gz = Path(cfg["ref_name"] + ".gz")
    ref_decomp = Path(cfg["ref_name"])
    ref_url = cfg["ref_url"]
    index_prefix = Path(cfg["index_prefix"] + ".idx")

    ensure_tools("dwgsim", "minimap2", "bowtie2", "bwa-mem2", "samtools", "wget", "pigz")
    download_and_decompress(ref_gz, ref_url, ref_decomp)
    ensure_bwa_index(ref_decomp)
    ensure_bowtie2_index(ref_decomp)
    ensure_sing_index(ref_decomp, index_prefix)

    filtered_ref = ensure_filtered_ref(ref_decomp)

    coverage_list = [1, 3]
    mut_rates = [0.001, 0.01]
    threads = threads_override

    output_csv = Path(f"benchmark_results_f1.{mode}.csv")
    if not output_csv.exists():
        output_csv.write_text("Exp_ID,Tool,Time_ms,Mem_kb,TotalReads,TP,FP,FN,Precision,Recall,F1,Mode\n")

    workdir = Path(f"bench_all_{mode}")
    ensure_dir(workdir)

    for coverage in coverage_list:
        for mut_rate in mut_rates:
            exp_id = f"Cov_{coverage}_Mut_{mut_rate}.{mode}"
            print("=================================================")
            print(f"Running: {exp_id} (Cov: {coverage}, Mut: {mut_rate})")
            print("=================================================")

            r1 = Path(f"sim_{exp_id}.bwa.read1.fastq.gz")
            r2 = Path(f"sim_{exp_id}.bwa.read2.fastq.gz")

            if not (r1.exists() and r2.exists()):
                run(["dwgsim", "-C", str(coverage), "-1", "150", "-2", "150", "-R", "0", "-X", "0", "-r", str(mut_rate), "-y", "0", "-H", str(filtered_ref), f"sim_{exp_id}"])

            tool_cmds = {
                "Minimap2": ["minimap2", "-t", str(threads), "-ax", "sr", str(ref_decomp), str(r1), str(r2)],
                "BWA-MEM2": ["bwa-mem2", "mem", "-t", str(threads), str(ref_decomp), str(r1), str(r2)],
                "Sing": ["./target/release/sing", "map", "-t", str(threads), str(index_prefix), "-1", str(r1), "-2", str(r2)],
                "Bowtie2": ["bowtie2", "-p", str(threads), "-x", str(ref_decomp), "-1", str(r1), "-2", str(r2)],
            }

            bam_paths = {}
            times_ms = {}
            mem_kb_map = {}
            for tool, cmd in tool_cmds.items():
                bam_path = workdir / f"{tool.lower().replace('-', '')}.bam"
                print(f"Running {tool}...")
                elapsed, mem_kb = run_mapper_to_sorted_bam(cmd, bam_path, threads)
                if elapsed is None:
                    times_ms[tool] = "N/A"
                    mem_kb_map[tool] = "N/A"
                else:
                    times_ms[tool] = int(elapsed * 1000)
                    mem_kb_map[tool] = mem_kb if mem_kb is not None else "N/A"
                bam_paths[tool] = bam_path

            truth_dict, tool_results = load_bams_and_build_truth(bam_paths)
            total_reads = len(truth_dict)

            for tool in ["Minimap2", "BWA-MEM2", "Sing", "Bowtie2"]:
                if tool in tool_results and len(tool_results[tool]) > 0:
                    tp, fp, fn, prec, rec, f1 = calculate_metrics(truth_dict, tool_results[tool])
                    csv_line = f"{exp_id},{tool},{times_ms[tool]},{mem_kb_map[tool]},{total_reads},{tp},{fp},{fn},{prec:.4f},{rec:.4f},{f1:.4f},{mode}"
                else:
                    csv_line = f"{exp_id},{tool},{times_ms[tool]},{mem_kb_map[tool]},{total_reads},0,0,0,0,0,0,{mode}"
                with output_csv.open("a") as f:
                    f.write(csv_line + "\n")

            print("\n--- Summary (ALL) ---")
            print("Tool        | Time_ms | Mem_KB  | Precision | Recall   | F1")
            for tool in ["Minimap2", "BWA-MEM2", "Sing", "Bowtie2"]:
                if tool in tool_results and len(tool_results[tool]) > 0:
                    tp, fp, fn, prec, rec, f1 = calculate_metrics(truth_dict, tool_results[tool])
                    print(f"{tool:<10} | {times_ms[tool]:<7} | {mem_kb_map[tool]:<7} | {prec:6.2f}%   | {rec:6.2f}%   | {f1:6.2f}")
                else:
                    print(f"{tool:<10} | {times_ms[tool]:<7} | {mem_kb_map[tool]:<7} | FAILED")

            for bam in bam_paths.values():
                bam.unlink(missing_ok=True)
                Path(str(bam) + ".bai").unlink(missing_ok=True)

    if output_csv.exists():
        print(f"CSV saved: {output_csv}")
    else:
        print(f"CSV missing: {output_csv}")

    if workdir.exists():
        for leftover in workdir.glob("*.unsorted.bam"):
            leftover.unlink(missing_ok=True)
        try:
            workdir.rmdir()
        except OSError:
            pass


def benchmark_minimap(mode, threads_override):
    cfg = get_mode_config(mode)
    ref_gz = Path(cfg["ref_name"] + ".gz")
    ref_decomp = Path(cfg["ref_name"])
    ref_url = cfg["ref_url"]
    index_prefix = Path(cfg["index_prefix"] + ".idx")

    ensure_tools("dwgsim", "minimap2", "samtools", "wget", "pigz")
    download_and_decompress(ref_gz, ref_url, ref_decomp)
    ensure_sing_index(ref_decomp, index_prefix)
    filtered_ref = ensure_filtered_ref(ref_decomp)

    reads_n = 100000
    mut_rates = [0.001, 0.01]
    threads = threads_override

    output_csv = Path(f"benchmark_results_minimap_comparision.{mode}.csv")
    if not output_csv.exists():
        output_csv.write_text("Exp_ID,Tool,Time_ms,Mem_kb,TotalReads,TP,FP,FN,Precision,Recall,F1,Mode\n")

    workdir = Path(f"bench_minimap_{mode}")
    ensure_dir(workdir)

    for mut_rate in mut_rates:
        exp_id = f"N_{reads_n}_Mut_{mut_rate}.{mode}"
        print("=================================================")
        print(f"Running: {exp_id} (N: {reads_n}, Mut: {mut_rate})")
        print("=================================================")

        r1 = Path(f"sim_{exp_id}.bwa.read1.fastq.gz")
        r2 = Path(f"sim_{exp_id}.bwa.read2.fastq.gz")

        if not (r1.exists() and r2.exists()):
            run(["dwgsim", "-N", str(reads_n), "-1", "150", "-2", "150", "-R", "0", "-X", "0", "-r", str(mut_rate), "-y", "0", "-H", str(filtered_ref), f"sim_{exp_id}"])

        tool_cmds = {
            "Minimap2": ["minimap2", "-t", str(threads), "-ax", "sr", str(ref_decomp), str(r1), str(r2)],
            "Sing": ["./target/release/sing", "map", "-t", str(threads), str(index_prefix), "-1", str(r1), "-2", str(r2)],
        }

        bam_paths = {}
        times_ms = {}
        mem_kb_map = {}
        for tool, cmd in tool_cmds.items():
            bam_path = workdir / f"{tool.lower().replace('-', '')}.bam"
            print(f"Running {tool}...")
            elapsed, mem_kb = run_mapper_to_sorted_bam(cmd, bam_path, threads)
            if elapsed is None:
                times_ms[tool] = "N/A"
                mem_kb_map[tool] = "N/A"
            else:
                times_ms[tool] = int(elapsed * 1000)
                mem_kb_map[tool] = mem_kb if mem_kb is not None else "N/A"
            bam_paths[tool] = bam_path

        truth_dict, tool_results = load_bams_and_build_truth(bam_paths)
        total_reads = len(truth_dict)

        for tool in ["Minimap2", "Sing"]:
            if tool in tool_results and len(tool_results[tool]) > 0:
                tp, fp, fn, prec, rec, f1 = calculate_metrics(truth_dict, tool_results[tool])
                csv_line = f"{exp_id},{tool},{times_ms[tool]},{mem_kb_map[tool]},{total_reads},{tp},{fp},{fn},{prec:.4f},{rec:.4f},{f1:.4f},{mode}"
            else:
                csv_line = f"{exp_id},{tool},{times_ms[tool]},{mem_kb_map[tool]},{total_reads},0,0,0,0,0,0,{mode}"
            with output_csv.open("a") as f:
                f.write(csv_line + "\n")

        print("\n--- Summary (MINIMAP) ---")
        print("Tool        | Time_ms | Mem_KB  | Precision | Recall   | F1")
        for tool in ["Minimap2", "Sing"]:
            if tool in tool_results and len(tool_results[tool]) > 0:
                tp, fp, fn, prec, rec, f1 = calculate_metrics(truth_dict, tool_results[tool])
                print(f"{tool:<10} | {times_ms[tool]:<7} | {mem_kb_map[tool]:<7} | {prec:6.2f}%   | {rec:6.2f}%   | {f1:6.2f}")
            else:
                print(f"{tool:<10} | {times_ms[tool]:<7} | {mem_kb_map[tool]:<7} | FAILED")

        for bam in bam_paths.values():
            bam.unlink(missing_ok=True)
            Path(str(bam) + ".bai").unlink(missing_ok=True)

    if output_csv.exists():
        print(f"CSV saved: {output_csv}")
    else:
        print(f"CSV missing: {output_csv}")

    if workdir.exists():
        for leftover in workdir.glob("*.unsorted.bam"):
            leftover.unlink(missing_ok=True)
        try:
            workdir.rmdir()
        except OSError:
            pass


def benchmark_strobe(mode, threads_override):
    cfg = get_mode_config(mode)
    ref_gz = Path(cfg["ref_name"] + ".gz")
    ref_decomp = Path(cfg["ref_name"])
    ref_url = cfg["ref_url"]
    index_prefix = Path(cfg["index_prefix"] + ".idx")

    ensure_tools("dwgsim", "samtools", "wget", "pigz")
    download_and_decompress(ref_gz, ref_url, ref_decomp)
    ensure_sing_index(ref_decomp, index_prefix)
    filtered_ref = ensure_filtered_ref(ref_decomp)

    reads_n = 100000
    mut_rates = [0.001, 0.01]
    threads = threads_override

    output_csv = Path(f"benchmark_results_custom.{mode}.csv")
    output_csv.write_text("Experiment,Tool,Time_sec,Mem_kb,Precision,Recall,F1\n")

    workdir = Path(f"bench_strobe_{mode}")
    ensure_dir(workdir)

    for mut_rate in mut_rates:
        exp_id = f"N_{reads_n}_Mut_{mut_rate}.{mode}"
        print("=================================================")
        print(f"Running: {exp_id} (N: {reads_n}, Mut: {mut_rate})")
        print("=================================================")

        r1 = Path(f"sim_{exp_id}.bwa.read1.fastq.gz")
        r2 = Path(f"sim_{exp_id}.bwa.read2.fastq.gz")

        if not (r1.exists() and r2.exists()):
            run(["dwgsim", "-N", str(reads_n), "-1", "150", "-2", "150", "-R", "0", "-X", "0", "-r", str(mut_rate), "-y", "0", "-H", str(filtered_ref), f"sim_{exp_id}"])

        tool_cmds = {
            "Sing": ["./target/release/sing", "map", "-t", str(threads), str(index_prefix), "-1", str(r1), "-2", str(r2)],
            "Strobealign": ["strobealign", "-t", str(threads), str(ref_decomp), str(r1), str(r2)],
        }

        bam_paths = {}
        times_sec = {}
        mem_kb_map = {}
        for tool, cmd in tool_cmds.items():
            bam_path = workdir / f"{tool.lower()}.bam"
            print(f"Running {tool}...")
            elapsed, mem_kb = run_mapper_to_sorted_bam(cmd, bam_path, threads)
            if elapsed is None:
                times_sec[tool] = "N/A"
                mem_kb_map[tool] = "N/A"
            else:
                times_sec[tool] = f"{elapsed:.2f}"
                mem_kb_map[tool] = mem_kb if mem_kb is not None else "N/A"
            bam_paths[tool] = bam_path

        truth_dict, tool_results = load_bams_and_build_truth(bam_paths)

        for tool in ["Sing", "Strobealign"]:
            if tool in tool_results and len(tool_results[tool]) > 0:
                tp, fp, fn, prec, rec, f1 = calculate_metrics(truth_dict, tool_results[tool])
                csv_line = f"{exp_id},{tool},{times_sec[tool]},{mem_kb_map[tool]},{prec:.4f},{rec:.4f},{f1:.4f}"
            else:
                csv_line = f"{exp_id},{tool},{times_sec[tool]},{mem_kb_map[tool]},0,0,0"
            with output_csv.open("a") as f:
                f.write(csv_line + "\n")

        print("\n--- Summary (STROBE) ---")
        print("Tool        | Time_sec | Mem_KB  | Precision | Recall   | F1")
        for tool in ["Sing", "Strobealign"]:
            if tool in tool_results and len(tool_results[tool]) > 0:
                tp, fp, fn, prec, rec, f1 = calculate_metrics(truth_dict, tool_results[tool])
                print(f"{tool:<10} | {times_sec[tool]:<8} | {mem_kb_map[tool]:<7} | {prec:6.2f}%   | {rec:6.2f}%   | {f1:6.2f}")
            else:
                print(f"{tool:<10} | {times_sec[tool]:<8} | {mem_kb_map[tool]:<7} | FAILED")

        for bam in bam_paths.values():
            bam.unlink(missing_ok=True)
            Path(str(bam) + ".bai").unlink(missing_ok=True)

    if output_csv.exists():
        print(f"CSV saved: {output_csv}")
    else:
        print(f"CSV missing: {output_csv}")

    if workdir.exists():
        for leftover in workdir.glob("*.unsorted.bam"):
            leftover.unlink(missing_ok=True)
        try:
            workdir.rmdir()
        except OSError:
            pass


def run_timed_cmd(label, cmd_str: str, log_path: Path):
    time_bin = Path("./time").resolve()
    if not time_bin.exists():
        time_bin = Path("/usr/bin/time")

    if not time_bin.exists():
        print(f"Error: {time_bin} not found. installation required.")
        return
    run_shell(f"{time_bin} -f '%e %M' -o '{log_path}' {cmd_str}")


def benchmark_gatk(mode, threads_override):
    cfg = get_mode_config(mode)
    ref_gz = Path(cfg["ref_name"] + ".gz")
    ref_decomp = Path(cfg["ref_name"])
    ref_url = cfg["ref_url"]

    ensure_tools("dwgsim", "samtools", "bwa-mem2", "minimap2", "bcftools", "gatk", "tabix", "wget", "pigz", "bgzip")

    download_and_decompress(ref_gz, ref_url, ref_decomp)

    out_dir = Path(f"{mode}.gatk")
    ensure_dir(out_dir)

    ref = out_dir / ref_decomp.name
    if not ref.exists():
        run_shell(f"sed 's/ .*//' <(pigz -p 8 -cd '{ref_gz}') > '{ref}'")

    r1 = out_dir / "sim_reads.bwa.read1.fastq.gz"
    r2 = out_dir / "sim_reads.bwa.read2.fastq.gz"
    raw_truth = out_dir / "sim_reads.mutations.vcf"

    coverage = 10
    read_len = 150
    mut_rate = 0.001
    indel_frac = 0.05
    err_rate = 0.01
    threads = threads_override

    if not (r1.exists() and r2.exists()):
        run([
            "dwgsim",
            "-C", str(coverage), "-1", str(read_len), "-2", str(read_len),
            "-r", str(mut_rate), "-R", str(indel_frac),
            "-e", str(err_rate), "-E", str(err_rate), "-y", "0",
            str(ref), str(out_dir / "sim_reads")
        ])

    truth_vcf = out_dir / "truth.vcf.gz"
    if not truth_vcf.exists():
        run_shell(f"bgzip -c '{raw_truth}' > '{truth_vcf}'")
        run(["tabix", "-p", "vcf", str(truth_vcf)])

    sing_index = out_dir / "index.sing"
    if not sing_index.exists():
        run(["cargo", "run", "--release", "--", "build", str(ref), str(sing_index)])

    if not Path(str(ref) + ".bwt.2bit.64").exists():
        run(["bwa-mem2", "index", "-p", str(ref), str(ref)])

    tools = ["sing", "bwa", "mini"]
    seconds_map = {}
    mem_map = {}

    def map_and_track(tool, cmd_list):
        bam = out_dir / f"{tool}.bam"
        if bam.exists():
            # If bam exists, try to recover stats from unsorted log if possible or just skip time
            return

        print(f"Running {tool}...")
        elapsed, mem_kb = run_mapper_to_sorted_bam(cmd_list, bam, threads)
        
        if elapsed is not None:
             seconds_map[tool] = f"{elapsed:.2f}"
             mem_map[tool] = str(mem_kb) if mem_kb is not None else "N/A"
        else:
             seconds_map[tool] = "N/A"
             mem_map[tool] = "N/A"

    if not (out_dir / "sing.bam").exists():
        map_and_track("sing", ["./target/release/sing", "map", "-t", str(threads), str(sing_index), "-1", str(r1), "-2", str(r2)])

    if not (out_dir / "bwa.bam").exists():
        # Note: No single quotes around RG string in list mode
        map_and_track("bwa", ["bwa-mem2", "mem", "-t", str(threads), "-R", f"@RG\\tID:bwa\\tSM:{cfg['species']}\\tPL:ILLUMINA", str(ref), str(r1), str(r2)])

    if not (out_dir / "mini.bam").exists():
        map_and_track("mini", ["minimap2", "-t", str(threads), "-ax", "sr", "-R", f"@RG\\tID:mini\\tSM:{cfg['species']}\\tPL:ILLUMINA", str(ref), str(r1), str(r2)])

    dict_path = Path(str(ref).replace(".fa", ".dict"))
    if ref.suffix != ".fa":
        dict_path = Path(str(ref) + ".dict")
    if not dict_path.exists():
        run(["gatk", "CreateSequenceDictionary", "-R", str(ref), "-O", str(dict_path)], quiet=True)
    if not Path(str(ref) + ".fai").exists():
        run(["samtools", "faidx", str(ref)])

    for tool in tools:
        gvcf = out_dir / f"{tool}.g.vcf.gz"
        final_vcf = out_dir / f"{tool}.final.vcf.gz"
        if not gvcf.exists():
            run(["gatk", "--java-options", "-Xmx16g", "HaplotypeCaller", "-R", str(ref), "-I", str(out_dir / f"{tool}.bam"), "-O", str(gvcf), "-ERC", "GVCF"], quiet=True)
        if not final_vcf.exists():
            run(["gatk", "--java-options", "-Xmx16g", "GenotypeGVCFs", "-R", str(ref), "-V", str(gvcf), "-O", str(final_vcf)], quiet=True)

    truth_norm = out_dir / "truth.norm.vcf.gz"
    if not truth_norm.exists():
        run(["bcftools", "norm", "-f", str(ref), "-m", "-any", "-O", "z", "-o", str(truth_norm), str(truth_vcf)])
        run(["tabix", "-p", "vcf", str(truth_norm)])

    tp_map = {}
    fp_map = {}
    fn_map = {}
    prec_map = {}
    rec_map = {}
    f1_map = {}

    for tool in tools:
        eval_dir = out_dir / f"eval_{tool}"
        ensure_dir(eval_dir)
        tool_norm = eval_dir / "tool.norm.vcf.gz"
        if not tool_norm.exists():
            run(["bcftools", "norm", "-f", str(ref), "-m", "-any", "-O", "z", "-o", str(tool_norm), str(out_dir / f"{tool}.final.vcf.gz")])
            run(["tabix", "-p", "vcf", str(tool_norm)])

        run(["bcftools", "isec", "-p", str(eval_dir), "-Oz", str(truth_norm), str(tool_norm)], quiet=True)

        fn = int(subprocess.check_output(["bash", "-lc", f"pigz -p 8 -cd '{eval_dir / '0000.vcf.gz'}' | grep -v '^#' | wc -l"]).strip())
        fp = int(subprocess.check_output(["bash", "-lc", f"pigz -p 8 -cd '{eval_dir / '0001.vcf.gz'}' | grep -v '^#' | wc -l"]).strip())
        tp = int(subprocess.check_output(["bash", "-lc", f"pigz -p 8 -cd '{eval_dir / '0002.vcf.gz'}' | grep -v '^#' | wc -l"]).strip())

        tp_map[tool] = tp
        fp_map[tool] = fp
        fn_map[tool] = fn
        if tp == 0:
            prec, rec, f1 = 0.0, 0.0, 0.0
        else:
            prec = tp / (tp + fp) * 100
            rec = tp / (tp + fn) * 100
            f1 = 2 * tp / (2 * tp + fp + fn) * 100
        prec_map[tool] = f"{prec:.2f}"
        rec_map[tool] = f"{rec:.2f}"
        f1_map[tool] = f"{f1:.2f}"

    output_csv = out_dir / f"benchmark_results_gatk.{mode}.csv"
    if not output_csv.exists():
        output_csv.write_text("Tool,Time_sec,Mem_kb,TP,FP,FN,Precision,Recall,F1\n")
    for tool in tools:
        with output_csv.open("a") as f:
            f.write(
                f"{tool},{seconds_map.get(tool, 'N/A')},{mem_map.get(tool, 'N/A')},{tp_map[tool]},{fp_map[tool]},{fn_map[tool]},{prec_map[tool]},{rec_map[tool]},{f1_map[tool]}\n"
            )

    print("\n##############################################################################")
    print(f"            SIMULATION REPORT [{cfg['species']}] (Clean {coverage}X, Real - No Filter)")
    print("##############################################################################")
    print("\n[1] Absolute Performance")
    print("------------------------------------------------------------------")
    print("Tool       | Time(s)     | Mem(KB)     | Speedup(vs BWA)")
    print("------------------------------------------------------------------")
    bwa_time = float(seconds_map.get("bwa", "0") or "0")
    for tool in tools:
        t = float(seconds_map.get(tool, "0") or "0")
        ratio = (bwa_time / t) if t > 0 else 0
        print(f"{tool:<10} | {seconds_map.get(tool, 'N/A'):<12} | {mem_map.get(tool, 'N/A'):<12} | {ratio:.2f}x")
    print("------------------------------------------------------------------")

    print("\n[2] Accuracy vs Ground Truth (Raw VCF - No Hard Filtering)")
    print("---------------------------------------------------------------------------------------")
    print("Tool       | TP       | FP       | FN       | Recall   | Precis   | F1")
    print("---------------------------------------------------------------------------------------")
    for tool in tools:
        print(f"{tool:<10} | {tp_map[tool]:<8} | {fp_map[tool]:<8} | {fn_map[tool]:<8} | {rec_map[tool]:<8}% | {prec_map[tool]:<8}% | {f1_map[tool]:<8}%")
    print("---------------------------------------------------------------------------------------")
    print(" * Condition: 10X, 0.1% Mut, 10% Indel, 1% Error")
    print(" * Filter: NONE (Raw Output from GenotypeGVCFs)")
    print("##############################################################################")

    if output_csv.exists():
        print(f"CSV saved: {output_csv}")
    else:
        print(f"CSV missing: {output_csv}")

    for tool in tools:
        (out_dir / f"{tool}.sam").unlink(missing_ok=True)
        (out_dir / f"{tool}.unsorted.bam").unlink(missing_ok=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("threads", type=int)
    parser.add_argument("benchmark", choices=["gatk", "all", "minimap", "strobe", "compile"])
    parser.add_argument("mode", nargs="?", choices=["y", "a", "m", "h", "b"])
    args = parser.parse_args()

    if args.mode is None:
        raise SystemExit("mode is required for this benchmark")

    if args.benchmark == "gatk":
        benchmark_gatk(args.mode, args.threads)
    elif args.benchmark == "all":
        benchmark_all(args.mode, args.threads)
    elif args.benchmark == "minimap":
        benchmark_minimap(args.mode, args.threads)
    elif args.benchmark == "strobe":
        benchmark_strobe(args.mode, args.threads)


if __name__ == "__main__":
    main()
