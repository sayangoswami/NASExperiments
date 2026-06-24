#!/usr/bin/env python3
"""
benchmark_aligner.py
====================
Standalone benchmark harness for readfish-compatible aligner plugins.

Queries a BLOW5 file (raw signal) or a FASTA/FASTQ file (sequence) against
any plugin that exposes the readfish Aligner interface, and writes per-read
alignment stats and batch throughput metrics.

Ground truth is parsed from the read ID using the squigulator/seq2squiggle
convention:  <prefix>!<contig>!<start>!<end>!<strand>
contig → true_label is resolved via an optional manifest TSV.

--input accepts any of:
  • A single file      reads.blow5
  • A glob pattern     "reads_d0.2_Comm_0/*.blow5"   (quote to avoid shell expansion)
  • A folder           reads_d0.2_Comm_0/             (all matching files inside)
  • Multiple values    --input a.blow5 --input b.blow5

All inputs must be the same type (all BLOW5 or all FASTA/Q).
Files are streamed in sorted order; batches span file boundaries seamlessly.

Usage
-----
# Single BLOW5 file
python benchmark_aligner.py \
    --input reads.blow5 \
    --plugin pyrawhash \
    --plugin-args "idx=/path/to/index.ind threads=16 x=bacterial" \
    --output results/rawhash_zymo/ \
    --manifest /data/zymo_manifest.tsv \
    --batch-size 512

# Glob pattern (quote it)
python benchmark_aligner.py \
    --input "reads_d0.2_Comm_0/*.blow5" \
    --plugin pyrawhash \
    --plugin-args "idx=/path/to/index.ind threads=16 x=bacterial" \
    --output results/rawhash_zymo/

# Folder
python benchmark_aligner.py \
    --input reads_d0.2_Comm_0/ \
    --plugin pyrawhash \
    --plugin-args "idx=/path/to/index.ind threads=16 x=bacterial" \
    --output results/rawhash_zymo/

# Multiple explicit files
python benchmark_aligner.py \
    --input chunk_000.blow5 --input chunk_001.blow5 --input chunk_002.blow5 \
    --plugin pyrawhash \
    --plugin-args "idx=/path/to/index.ind threads=16 x=bacterial" \
    --output results/rawhash_zymo/

# FASTA/Q input, sequence aligner
python benchmark_aligner.py \
    --input reads.fastq \
    --plugin mappy \
    --plugin-args "fn_idx_in=/path/to/ref.fa n_threads=16" \
    --output results/mappy_zymo/ \
    --manifest /data/zymo_manifest.tsv \
    --batch-size 512
"""

from __future__ import annotations

import argparse
import csv
import importlib
import json
import logging
import re
import sys
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, Iterable, Iterator, List, Optional, Tuple

import numpy as np

import glob

logger = logging.getLogger("benchmark_aligner")

# ---------------------------------------------------------------------------
# Input resolution: file, glob, folder, or multiple --input values
# ---------------------------------------------------------------------------

_BLOW5_EXTS = {".blow5", ".slow5"}
_FASTA_EXTS = {".fa", ".fasta", ".fq", ".fastq"}

def resolve_inputs(raw_inputs: List[str]) -> Tuple[List[Path], str]:
    """
    Accept one or more of:
      - a single file path
      - a glob pattern  (e.g. "reads/*.blow5")
      - a directory     (all BLOW5 or FASTA/Q files inside, sorted)
      - any combination via repeated --input flags

    Returns (sorted_file_list, input_type) where input_type is 'blow5' or 'fasta'.
    Raises if the list is empty or contains mixed types.
    """
    paths: List[Path] = []
    for raw in raw_inputs:
        p = Path(raw)
        if p.is_dir():
            # Collect all recognised files in the directory (non-recursive)
            found = sorted(
                f for f in p.iterdir()
                if f.is_file() and f.suffix.lower() in _BLOW5_EXTS | _FASTA_EXTS
            )
            if not found:
                raise ValueError(f"Directory {p} contains no BLOW5 or FASTA/Q files.")
            paths.extend(found)
        else:
            # Try as glob (also handles plain file paths since glob("exact") == [exact])
            matched = sorted(Path(g) for g in glob.glob(raw))
            if not matched:
                raise FileNotFoundError(f"No files matched: {raw!r}")
            paths.extend(matched)

    if not paths:
        raise ValueError("No input files found.")

    # Determine and validate uniform type
    types = set()
    for p in paths:
        ext = p.suffix.lower()
        if ext in _BLOW5_EXTS:
            types.add("blow5")
        elif ext in _FASTA_EXTS:
            types.add("fasta")
        else:
            raise ValueError(f"Unrecognised file extension: {p.suffix!r} ({p})")

    if len(types) > 1:
        raise ValueError(
            f"Mixed input types: found both BLOW5 and FASTA/Q files. "
            f"All inputs must be the same type."
        )

    input_type = types.pop()
    logger.info("Resolved %d input file(s) [%s]:", len(paths), input_type)
    for p in paths:
        logger.info("  %s", p)
    return paths, input_type


# ---------------------------------------------------------------------------
# Read-ID parsing  (squigulator / seq2squiggle convention)
# ---------------------------------------------------------------------------

def parse_read_id(read_id: str) -> dict:
    """
    Parse a read ID of the form:
        seq2squiggle: <contig>!<start>!<end>!<strand>          (4 fields)
        squigulator:  <prefix>!<contig>!<start>!<end>!<strand> (5 fields)
    parts[-4] is always the contig name.
    """
    parts = read_id.split("!")
    if len(parts) >= 4:
        contig, start, end, strand = parts[-4], parts[-3], parts[-2], parts[-1]
        prefix = "!".join(parts[:-4]) if len(parts) > 4 else ""
        species = re.sub(r"_contig_\d+$", "", contig).replace("_", " ")
        return {
            "prefix":  prefix,
            "contig":  contig,
            "species": species,
            "start":   int(start) if start.lstrip("-").isdigit() else None,
            "end":     int(end)   if end.lstrip("-").isdigit()   else None,
            "strand":  strand,
        }
    return {"prefix": read_id, "contig": None, "species": None,
            "start": None, "end": None, "strand": None}


# ---------------------------------------------------------------------------
# Manifest loading + ground truth resolution
# ---------------------------------------------------------------------------

def load_manifest(path: str) -> Dict[str, int]:
    """Load a contig/species → true_label manifest TSV."""
    mapping: Dict[str, int] = {}
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        key_col = None
        for col in ("contig", "species"):
            if col in (reader.fieldnames or []):
                key_col = col
                break
        if key_col is None:
            raise ValueError(f"Manifest {path} must have a 'contig' or 'species' column.")
        if "true_label" not in (reader.fieldnames or []):
            raise ValueError(f"Manifest {path} must have a 'true_label' column.")
        for row in reader:
            mapping[row[key_col].strip()] = int(row["true_label"])
    logger.info("Manifest: %d entries loaded from %s", len(mapping), path)
    return mapping


def resolve_true_label(contig: Optional[str],
                       manifest: Optional[Dict[str, int]]) -> int:
    """Return 1 (target), 0 (deplete), or -1 (unknown)."""
    if manifest is None or contig is None:
        return -1
    if contig in manifest:
        return manifest[contig]
    species = re.sub(r"_contig_\d+$", "", contig).replace("_", " ")
    if species in manifest:
        return manifest[species]
    if species.replace(" ", "_") in manifest:
        return manifest[species.replace(" ", "_")]
    return -1


# ---------------------------------------------------------------------------
# Plugin loading
# ---------------------------------------------------------------------------

def load_plugin(module_name: str, plugin_args: str, debug_log: Optional[str]):
    """
    Import <module_name> and instantiate <module_name>.Aligner(**kwargs).
    plugin_args is a whitespace-separated string of key=value pairs,
    matching the readfish toml mapper_settings convention.
    """
    kwargs: Dict[str, object] = {}
    for token in plugin_args.split():
        if "=" in token:
            k, v = token.split("=", 1)
            # Coerce booleans and integers
            if v.lower() == "true":
                kwargs[k] = True
            elif v.lower() == "false":
                kwargs[k] = False
            else:
                try:
                    kwargs[k] = int(v)
                except ValueError:
                    try:
                        kwargs[k] = float(v)
                    except ValueError:
                        kwargs[k] = v
        else:
            logger.warning("Ignoring malformed plugin arg token: %s", token)

    mod = importlib.import_module(module_name)
    aligner_cls = getattr(mod, "Aligner")
    logger.info("Instantiating %s.Aligner with kwargs=%s", module_name, kwargs)
    aligner = aligner_cls(debug_log=debug_log, **kwargs)
    aligner.validate()
    return aligner


# ---------------------------------------------------------------------------
# Result dataclass  (mirrors readfish plugin utils)
# ---------------------------------------------------------------------------

@dataclass
class Result:
    channel: int
    read_id: str
    seq: str
    barcode: Optional[str] = None
    basecall_data: Optional[object] = None   # numpy array for signal aligners
    alignment_data: Optional[list] = None


# ---------------------------------------------------------------------------
# Input readers
# ---------------------------------------------------------------------------

def _iter_blow5(path: str) -> Iterator[Tuple[str, np.ndarray]]:
    """
    Yield (read_id, signal_array) from a BLOW5 file using pyslow5.
    signal is returned as float32 (calibrated if available, else raw int16).
    """
    try:
        import pyslow5
    except ImportError:
        logger.error("pyslow5 is required for BLOW5 input: pip install pyslow5")
        sys.exit(1)

    s5 = pyslow5.Open(path, "r")
    for read in s5.seq_reads(aux="all"):
        read_id = read["read_id"]
        raw     = np.array(read["signal"], dtype=np.int16)

        # Apply calibration to produce float32 pA signal
        offset  = read.get("offset",  0.0)
        scaling = read.get("range",   1.0) / read.get("digitisation", 1.0) \
                  if "range" in read and "digitisation" in read else 1.0
        signal  = (raw.astype(np.float32) + float(offset)) * float(scaling)
        yield read_id, signal
    s5.close()


def _iter_fasta(path: str) -> Iterator[Tuple[str, str]]:
    """
    Yield (read_id, sequence) from a FASTA or FASTQ file.
    Uses pysam if available, otherwise a minimal built-in parser.
    """
    path = Path(path)
    suffix = path.suffix.lower()
    is_fastq = suffix in (".fastq", ".fq")

    try:
        import pysam
        with pysam.FastxFile(str(path)) as fh:
            for entry in fh:
                yield entry.name.split()[0], entry.sequence
        return
    except ImportError:
        pass

    # Fallback: plain parser
    with open(path) as fh:
        if is_fastq:
            while True:
                header = fh.readline().strip()
                if not header:
                    break
                seq  = fh.readline().strip()
                fh.readline()   # +
                fh.readline()   # quality
                yield header[1:].split()[0], seq
        else:
            read_id, seq_parts = None, []
            for line in fh:
                line = line.strip()
                if line.startswith(">"):
                    if read_id is not None:
                        yield read_id, "".join(seq_parts)
                    read_id = line[1:].split()[0]
                    seq_parts = []
                elif read_id is not None:
                    seq_parts.append(line)
            if read_id is not None:
                yield read_id, "".join(seq_parts)


def make_result_iter(paths: List[Path],
                     input_type: str,
                     batch_size: int) -> Iterator[List[Result]]:
    """
    Yield batches of Result objects from one or more BLOW5 or FASTA/Q files.
    Files are consumed in order; batches span file boundaries seamlessly.
    channel is a rolling counter mod 512 (dummy pore assignment).
    """
    channel = 0
    batch: List[Result] = []

    def _flush(b):
        if b:
            yield b[:]
            b.clear()

    for path in paths:
        logger.info("Reading %s", path)
        if input_type == "blow5":
            src = _iter_blow5(str(path))
            for read_id, signal in src:
                batch.append(Result(
                    channel=channel % 512,
                    read_id=read_id,
                    seq="N",   # dummy — prevents Readfish no_seq fallback
                    basecall_data=signal.astype(np.float32),
                ))
                channel += 1
                if len(batch) >= batch_size:
                    yield batch
                    batch = []
        else:
            src = _iter_fasta(str(path))
            for read_id, seq in src:
                batch.append(Result(
                    channel=channel % 512,
                    read_id=read_id,
                    seq=seq,
                    basecall_data=None,
                ))
                channel += 1
                if len(batch) >= batch_size:
                    yield batch
                    batch = []

    if batch:
        yield batch


# ---------------------------------------------------------------------------
# Per-read stats
# ---------------------------------------------------------------------------

@dataclass
class ReadStat:
    read_id:     str
    contig:      str
    true_label:  int           # 1=target, 0=deplete, -1=unknown
    mapped:      bool
    mapped_ctg:  str
    mapq:        float         # pres_frac from Alignment, or 0
    seq_len:     int           # len(seq) for sequence aligners, signal len for raw
    latency_s:   float         # wall time for this read's alignment


# ---------------------------------------------------------------------------
# Readfish-compatible decisions TSV writer
# ---------------------------------------------------------------------------

_DECISIONS_TSV_FIELDNAMES = [
    "client_iteration", "read_in_loop", "read_id", "channel",
    "seq_len", "counter", "mode", "decision", "condition",
    "barcode", "previous_action", "action_override", "timestamp",
]


def _result_to_decision(result: Result, input_type: str) -> str:
    """
    Map alignment_data to a Readfish decision string.
      aligned      → stop_receiving  (target hit — keep reading)
      no hit       → unblock         (no match — eject)
      no signal    → no_seq          (skipped before reaching aligner)
    """
    if result.alignment_data is None:
        return "no_seq"
    if result.alignment_data:
        return "stop_receiving"
    return "unblock"


def _seq_len_for_tsv(result: Result, input_type: str) -> int:
    """Return a meaningful seq_len for the TSV: signal samples for BLOW5,
    sequence length for FASTA/Q. Dummy 'N' counts as 0."""
    if input_type == "blow5":
        return len(result.basecall_data) if result.basecall_data is not None else 0
    seq = result.seq or ""
    return len(seq) if seq != "N" else 0


class DecisionTSVWriter:
    """
    Writes per-read decisions in the exact readfish_decisions.tsv column order
    so that analyse_runs.py can consume standalone benchmark output directly.
    """

    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self._fh = open(path, "w", newline="", buffering=1)  # line-buffered
        self._writer = csv.DictWriter(
            self._fh, fieldnames=_DECISIONS_TSV_FIELDNAMES, delimiter="\t"
        )
        self._writer.writeheader()

    def write_batch(
        self,
        results:          List[Result],
        input_type:       str,
        client_iteration: int,
        batch_start_time: float,
        result_times:     Dict[str, float],  # read_id → perf_counter at decision
    ) -> None:
        for read_in_loop, result in enumerate(results, start=1):
            decision = _result_to_decision(result, input_type)
            self._writer.writerow({
                "client_iteration": client_iteration,
                "read_in_loop":     read_in_loop,
                "read_id":          result.read_id,
                "channel":          result.channel,
                "seq_len":          _seq_len_for_tsv(result, input_type),
                "counter":          1,           # no chunk counting in standalone mode
                "mode":             "channel",
                "decision":         decision,
                "condition":        "analysis",
                "barcode":          "None",
                "previous_action":  "None",
                "action_override":  False,
                "timestamp":        result_times.get(result.read_id, batch_start_time),
            })

    def close(self) -> None:
        self._fh.flush()
        self._fh.close()

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.close()


# ---------------------------------------------------------------------------
# Core benchmark loop
# ---------------------------------------------------------------------------

def run_benchmark(
    aligner,
    input_paths:   List[Path],
    input_type:    str,
    batch_size:    int,
    manifest:      Optional[Dict[str, int]],
    outdir:        Path,
    tool:          str,
    dataset:       str,
    decisions_tsv: Optional[Path] = None,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    per_read_path  = outdir / f"{tool}_{dataset}_per_read.tsv"
    batch_log_path = outdir / f"{tool}_{dataset}_batch_stats.tsv"

    per_read_rows:  List[dict] = []
    batch_rows:     List[dict] = []

    total_reads  = 0
    total_mapped = 0

    dec_writer = DecisionTSVWriter(decisions_tsv) if decisions_tsv else None
    if dec_writer:
        logger.info("Writing decisions TSV → %s", decisions_tsv)

    try:
      with open(per_read_path,  "w", newline="") as prf, \
           open(batch_log_path, "w", newline="") as blf:

        pr_writer = csv.DictWriter(prf, fieldnames=[
            "read_id", "contig", "true_label", "mapped",
            "mapped_ctg", "mapq", "seq_len", "latency_s",
        ], delimiter="\t")
        pr_writer.writeheader()

        bl_writer = csv.DictWriter(blf, fieldnames=[
            "batch_idx", "n_reads", "n_mapped",
            "wall_s", "throughput_reads_per_s",
            "mean_latency_s", "median_latency_s",
            "p95_latency_s", "p99_latency_s",
        ], delimiter="\t")
        bl_writer.writeheader()

        for batch_idx, batch in enumerate(
                make_result_iter(input_paths, input_type, batch_size)):

            n_reads  = len(batch)
            t_batch_start = time.perf_counter()
            wall_time_start = time.time()  # real timestamp for TSV

            # Time each read individually inside the batch
            read_times:  List[float] = []
            read_stats:  List[ReadStat] = []

            # Wrap the batch in a generator that records per-read start times
            start_times: Dict[str, float] = {}
            result_timestamps: Dict[str, float] = {}  # read_id → epoch time

            def _timed_iter(results):
                for r in results:
                    start_times[r.read_id] = time.perf_counter()
                    yield r

            results = list(aligner.map_reads(_timed_iter(batch)))

            t_batch_end = time.perf_counter()
            wall_s = t_batch_end - t_batch_start

            # Assign a monotonically increasing epoch timestamp to each result
            # within the batch, interpolated across the batch wall time.
            n_results = len(results)
            for i, result in enumerate(results):
                frac = (i + 1) / n_results if n_results > 0 else 1.0
                result_timestamps[result.read_id] = (
                    wall_time_start + frac * (time.time() - wall_time_start)
                )

            for result in results:
                t_start = start_times.get(result.read_id, t_batch_start)
                latency = t_batch_end - t_start

                parsed     = parse_read_id(result.read_id)
                contig     = parsed.get("contig") or ""
                true_label = resolve_true_label(contig, manifest)
                mapped     = bool(result.alignment_data)
                mapped_ctg = result.alignment_data[0].ctg \
                             if mapped else "*"
                mapq       = result.alignment_data[0].pres_frac \
                             if mapped and hasattr(result.alignment_data[0], "pres_frac") \
                             else 0.0
                seq_len    = _seq_len_for_tsv(result, input_type)

                rs = ReadStat(
                    read_id=result.read_id,
                    contig=contig,
                    true_label=true_label,
                    mapped=mapped,
                    mapped_ctg=mapped_ctg,
                    mapq=mapq,
                    seq_len=seq_len,
                    latency_s=latency,
                )
                read_stats.append(rs)
                read_times.append(latency)

            # Write per-read rows
            for rs in read_stats:
                pr_writer.writerow(asdict(rs))

            # Write readfish-format decisions TSV if requested
            if dec_writer is not None:
                dec_writer.write_batch(
                    results=results,
                    input_type=input_type,
                    client_iteration=batch_idx + 1,
                    batch_start_time=wall_time_start,
                    result_times=result_timestamps,
                )

            # Batch summary
            n_mapped    = sum(1 for rs in read_stats if rs.mapped)
            lat         = np.array(read_times)
            throughput  = n_reads / wall_s if wall_s > 0 else float("nan")

            bl_row = {
                "batch_idx":              batch_idx,
                "n_reads":                n_reads,
                "n_mapped":               n_mapped,
                "wall_s":                 round(wall_s, 6),
                "throughput_reads_per_s": round(throughput, 2),
                "mean_latency_s":         round(float(np.mean(lat)),           6),
                "median_latency_s":       round(float(np.median(lat)),         6),
                "p95_latency_s":          round(float(np.percentile(lat, 95)), 6),
                "p99_latency_s":          round(float(np.percentile(lat, 99)), 6),
            }
            bl_writer.writerow(bl_row)

            total_reads  += n_reads
            total_mapped += n_mapped

            logger.info(
                "Batch %4d | %4d reads | %4d mapped | %.3fs | %.1f reads/s",
                batch_idx, n_reads, n_mapped, wall_s, throughput,
            )

    finally:
        if dec_writer is not None:
            dec_writer.close()

    # Summary JSON
    summary = {
        "tool":          tool,
        "dataset":       dataset,
        "input_type":    input_type,
        "total_reads":   total_reads,
        "total_mapped":  total_mapped,
        "map_rate":      round(total_mapped / total_reads, 4) if total_reads else 0.0,
    }
    summary_path = outdir / f"{tool}_{dataset}_summary.json"
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2)

    logger.info("Done. %d reads, %d mapped (%.1f%%)",
                total_reads, total_mapped,
                100 * summary["map_rate"])
    logger.info("Per-read TSV  → %s", per_read_path)
    logger.info("Batch log     → %s", batch_log_path)
    logger.info("Summary JSON  → %s", summary_path)
    print(json.dumps(summary, indent=2))


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input",       required=True, action="append", dest="inputs",
                   metavar="PATH_OR_GLOB",
                   help="Input file, folder, or glob pattern. "
                        "Repeat to supply multiple sources, e.g. "
                        "--input a.blow5 --input b.blow5 or --input 'reads/*.blow5'. "
                        "All inputs must be the same type (BLOW5 or FASTA/Q).")
    p.add_argument("--plugin",      required=True,
                   help="Python module name of the aligner plugin, e.g. pyrawhash, mappy")
    p.add_argument("--plugin-args", default="",
                   help='Space-separated key=value args passed to Aligner(), '
                        'e.g. "idx=/path/to/index threads=16 x=bacterial"')
    p.add_argument("--output",      required=True,
                   help="Output directory for TSV and JSON results")
    p.add_argument("--tool",        default=None,
                   help="Tool label for output filenames (default: --plugin value)")
    p.add_argument("--dataset",     default="dataset",
                   help="Dataset label for output filenames (default: dataset)")
    p.add_argument("--manifest",    default=None,
                   help="Optional TSV mapping contig/species → true_label")
    p.add_argument("--batch-size",  type=int, default=512,
                   help="Number of reads per batch (default: 512)")
    p.add_argument("--decisions-tsv", default=None,
                   help="If set, write a readfish-compatible readfish_decisions.tsv "
                        "to this path so analyse_runs.py can consume it directly")
    p.add_argument("--log-level",   default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p


def main():
    parser = build_parser()
    args   = parser.parse_args()

    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s  %(message)s",
        datefmt="%H:%M:%S",
    )

    # Resolve all inputs → sorted file list + uniform type
    try:
        input_paths, input_type = resolve_inputs(args.inputs)
    except (ValueError, FileNotFoundError) as e:
        logger.error("%s", e)
        sys.exit(1)

    tool    = args.tool or args.plugin
    outdir  = Path(args.output)
    manifest = load_manifest(args.manifest) if args.manifest else None

    aligner = load_plugin(args.plugin, args.plugin_args, args.debug_log)

    run_benchmark(
        aligner=aligner,
        input_paths=input_paths,
        input_type=input_type,
        batch_size=args.batch_size,
        manifest=manifest,
        outdir=outdir,
        tool=tool,
        dataset=args.dataset,
        decisions_tsv=Path(args.decisions_tsv) if args.decisions_tsv else None,
    )


if __name__ == "__main__":
    main()