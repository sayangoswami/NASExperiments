"""
bootstrap_ci.py — Bootstrap 95% CIs for F1, precision, recall, specificity.

Usage:
    python bootstrap_ci.py [options] <tsv_file_or_dir> [...]
    python bootstrap_ci.py [options] --list filelist.txt

Options:
    --list FILE       Read input paths from FILE (one per line; # = comment)
    --workers N       Parallel workers (default: 8)
    --nboot N         Total bootstrap iterations (default: 10000)
    --mem-fraction F  Fraction of available RAM to use (default: 0.75)

Chunk size is computed automatically from actual file row counts and
available RAM — no manual tuning needed.

Output: TSV to stdout. Redirect: python bootstrap_ci.py --list f.txt > out.tsv
"""

# Suppress numpy's internal thread pools before importing numpy.
# Without this, each worker spawns as many threads as there are cores,
# causing severe oversubscription when running many workers in parallel.
import os
os.environ["OMP_NUM_THREADS"]       = "1"
os.environ["OPENBLAS_NUM_THREADS"]  = "1"
os.environ["MKL_NUM_THREADS"]       = "1"
os.environ["VECLIB_MAXIMUM_THREADS"]= "1"
os.environ["NUMEXPR_NUM_THREADS"]   = "1"

import sys
import glob
import argparse
import numpy as np
import pandas as pd
from concurrent.futures import ProcessPoolExecutor, as_completed
from multiprocessing import Manager
from tqdm import tqdm

try:
    import psutil
    def available_ram_bytes():
        return psutil.virtual_memory().available
except ImportError:
    def available_ram_bytes():
        # Fallback: read /proc/meminfo (Linux)
        try:
            with open("/proc/meminfo") as f:
                for line in f:
                    if line.startswith("MemAvailable:"):
                        return int(line.split()[1]) * 1024
        except Exception:
            pass
        return 64 * 1024**3  # conservative 64 GB fallback

# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------
def load_tsv(path):
    return pd.read_csv(path, sep="\t", usecols=["decision", "true_label"],
                       dtype={"decision": str, "true_label": np.int8})

def load_input(arg):
    if os.path.isdir(arg):
        files = sorted(glob.glob(os.path.join(arg, "*.tsv")))
        if not files:
            raise FileNotFoundError(f"No .tsv files in {arg}")
        return pd.concat([load_tsv(f) for f in files], ignore_index=True)
    return load_tsv(arg)

def count_rows(arg):
    """Row count via wc -l (minus header line)."""
    import subprocess
    if os.path.isdir(arg):
        files = sorted(glob.glob(os.path.join(arg, "*.tsv")))
        return sum(count_rows(f) for f in files)
    result = subprocess.run(["wc", "-l", arg], capture_output=True, text=True)
    return int(result.stdout.split()[0]) - 1  # subtract header

# ---------------------------------------------------------------------------
# Bootstrap
# ---------------------------------------------------------------------------
def bootstrap(decisions, true_labels, n_boot, chunk_size, label="", position=0):
    rng  = np.random.default_rng(42)
    n    = len(decisions)
    keep = (decisions == "stop_receiving").astype(np.uint8)
    pos  = (true_labels == 1).astype(np.uint8)

    n_chunks    = (n_boot + chunk_size - 1) // chunk_size
    all_metrics = []
    remaining   = n_boot

    bar = tqdm(total=n_chunks, desc=f"  {label}", unit="chunk",
               position=position, leave=False, dynamic_ncols=True)

    while remaining > 0:
        batch = min(chunk_size, remaining)
        idx   = rng.integers(0, n, size=(batch, n), dtype=np.int32)

        k = keep[idx].astype(bool)
        p = pos[idx].astype(bool)

        TP = ( k &  p).sum(axis=1, dtype=np.int32)
        FP = ( k & ~p).sum(axis=1, dtype=np.int32)
        FN = (~k &  p).sum(axis=1, dtype=np.int32)
        TN = (~k & ~p).sum(axis=1, dtype=np.int32)

        prec = np.where(TP+FP > 0, TP/(TP+FP), np.nan)
        rec  = np.where(TP+FN > 0, TP/(TP+FN), np.nan)
        spec = np.where(TN+FP > 0, TN/(TN+FP), np.nan)
        f1   = np.where(prec+rec > 0, 2*prec*rec/(prec+rec), np.nan)

        all_metrics.append(np.stack([prec, rec, spec, f1], axis=1))
        remaining -= batch
        bar.update(1)

    bar.close()
    metrics = np.concatenate(all_metrics)

    # Point estimates on full data
    tp = int(((keep==1) & (pos==1)).sum())
    fp = int(((keep==1) & (pos==0)).sum())
    fn = int(((keep==0) & (pos==1)).sum())
    tn = int(((keep==0) & (pos==0)).sum())
    pr = tp/(tp+fp) if tp+fp else np.nan
    re = tp/(tp+fn) if tp+fn else np.nan
    sp = tn/(tn+fp) if tn+fp else np.nan
    f1 = 2*pr*re/(pr+re) if pr+re else np.nan

    lo = np.nanpercentile(metrics, 2.5,  axis=0)
    hi = np.nanpercentile(metrics, 97.5, axis=0)
    return (pr, re, sp, f1), lo, hi

# ---------------------------------------------------------------------------
# Per-file worker
# ---------------------------------------------------------------------------
def process_file(arg, n_boot, chunk_size, position_queue):
    arg   = arg.rstrip("/")
    label = arg
    position = position_queue.get()   # claim a terminal row
    try:
        df          = load_input(arg)
        decisions   = df["decision"].to_numpy()
        true_labels = df["true_label"].to_numpy()
        point, lo, hi = bootstrap(decisions, true_labels, n_boot, chunk_size,
                                   label=label, position=position)
    finally:
        position_queue.put(position)  # release row for next file
    pr, re, sp, f1 = point

    def fmt(v, l, h):
        return f"{v:.4f} [{l:.4f}, {h:.4f}]"

    return "\t".join([
        label, str(len(df)),
        fmt(pr, lo[0], hi[0]),
        fmt(re, lo[1], hi[1]),
        fmt(sp, lo[2], hi[2]),
        fmt(f1, lo[3], hi[3]),
    ])

# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("inputs", nargs="*")
    p.add_argument("--list",         metavar="FILE")
    p.add_argument("--workers",      type=int,   default=8)
    p.add_argument("--nboot",        type=int,   default=10_000)
    p.add_argument("--mem-fraction", type=float, default=0.75)
    return p.parse_args()

def compute_chunk_size(n_reads, n_workers, mem_fraction):
    """Derive chunk size from actual row count and available RAM.

    Memory per worker per chunk:
      idx matrix : chunk × n × 4 bytes (int32)
      k, p arrays: chunk × n × 1 byte each (uint8)
    Total: chunk × n × 6 bytes per worker.
    """
    available  = available_ram_bytes() * mem_fraction
    per_worker = available / n_workers
    chunk      = int(per_worker / (n_reads * 6))
    return max(1, min(chunk, 10_000))  # clamp: at least 1, at most 10k

def main():
    args = parse_args()

    if args.list:
        paths = []
        with open(args.list) as fh:
            for line in fh:
                line = line.strip()
                if line and not line.startswith("#"):
                    paths.append(line)
    elif args.inputs:
        paths = args.inputs
    else:
        print(__doc__); sys.exit(1)

    # Pre-scan row counts (fast — byte counting, no column parsing)
    tqdm.write("# Scanning file sizes...", file=sys.stderr)
    row_counts = {
        arg: count_rows(arg)
        for arg in tqdm(paths, desc="Scanning", unit="file",
                        dynamic_ncols=True, file=sys.stderr)
    }

    max_n      = max(row_counts.values())
    chunk_size = compute_chunk_size(max_n, args.workers, args.mem_fraction)
    ram_gb     = args.workers * chunk_size * max_n * 6 / 1e9

    tqdm.write(
        f"# {len(paths)} files | {args.workers} workers | "
        f"chunk={chunk_size} (auto) | nboot={args.nboot} | "
        f"max_n={max_n:,} | ~{ram_gb:.1f} GB peak RAM",
        file=sys.stderr
    )

    header = "\t".join(["Input","N_reads","Precision [95% CI]",
                         "Recall [95% CI]","Specificity [95% CI]","F1 [95% CI]"])
    tqdm.write(header)

    # Allocate terminal positions: position 0 = outer bar, 1..n_workers = inner bars
    with Manager() as manager:
        position_queue = manager.Queue()
        for i in range(1, args.workers + 1):
            position_queue.put(i)

        with ProcessPoolExecutor(max_workers=args.workers) as executor:
            futures = {
                executor.submit(process_file, arg, args.nboot,
                                chunk_size, position_queue): arg
                for arg in paths
            }
            bar = tqdm(as_completed(futures), total=len(paths), position=0,
                       desc="Files", unit="file", dynamic_ncols=True)
            for future in bar:
                arg = futures[future]
                bar.set_postfix(done=arg.rstrip("/"))
                try:
                    tqdm.write(future.result())
                except Exception as e:
                    tqdm.write(f"ERROR\t{arg}\t{e}", file=sys.stderr)

if __name__ == "__main__":
    main()