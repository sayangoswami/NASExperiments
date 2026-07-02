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

On-the-fly basecalling
----------------------
When a BLOW5 file is supplied but the aligner plugin works on sequences
(not raw signal), pass the Dorado basecall server parameters to basecall
each batch before querying:

    --basecall-address  ipc:///var/lib/minknow/data/.dorado/dorado-basecall-server.sock
    --basecall-config   dna_r10.4.1_e8.2_400bps_fast@v5.2.0||
    [--basecall-timeout  60]

Requires: ont-pybasecall-client-lib (version must match the Dorado server).

Signal / sequence truncation
-----------------------------
To simulate the 0.4-0.8 s MinKNOW API chunk cadence (1600-3200 raw samples
≈ 180-360 bases), optionally truncate all reads to a fixed length before
basecalling or querying:

    --truncate-signals N   keep only the first N raw signal samples
    --truncate-bases   N   keep only the first N sequence bases

Only one of the two flags may be set.  Truncation is applied:
  • BLOW5 input, raw-signal aligner   →  --truncate-signals  (before map_reads)
  • BLOW5 input, basecalling enabled  →  --truncate-signals  (before Dorado)
  • FASTA/Q input                     →  --truncate-bases    (before map_reads)

Usage
-----
# Single BLOW5 file, raw-signal aligner
python benchmark_aligner.py \\
    --input reads.blow5 \\
    --plugin pyrawhash \\
    --plugin-args "idx=/path/to/index.ind threads=16 x=bacterial" \\
    --output results/rawhash_zymo/ \\
    --manifest /data/zymo_manifest.tsv \\
    --truncate-signals 3200

# BLOW5 file, sequence aligner (basecall on-the-fly)
python benchmark_aligner.py \\
    --input reads.blow5 \\
    --plugin mappy \\
    --plugin-args "fn_idx_in=/path/to/ref.fa n_threads=16" \\
    --output results/mappy_zymo/ \\
    --manifest /data/zymo_manifest.tsv \\
    --basecall-address ipc:///var/lib/minknow/data/.dorado/dorado-basecall-server.sock \\
    --basecall-config  dna_r10.4.1_e8.2_400bps_fast@v5.2.0|| \\
    --truncate-signals 3200

# Glob pattern (quote it)
python benchmark_aligner.py \\
    --input "reads_d0.2_Comm_0/*.blow5" \\
    --plugin pyrawhash \\
    --plugin-args "idx=/path/to/index.ind threads=16 x=bacterial" \\
    --output results/rawhash_zymo/

# FASTA/Q input, sequence aligner
python benchmark_aligner.py \\
    --input reads.fastq \\
    --plugin mappy \\
    --plugin-args "fn_idx_in=/path/to/ref.fa n_threads=16" \\
    --output results/mappy_zymo/ \\
    --manifest /data/zymo_manifest.tsv \\
    --truncate-bases 360
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
from concurrent.futures import ThreadPoolExecutor, Future

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
# On-the-fly basecalling
# ---------------------------------------------------------------------------

@dataclass
class BasecallConfig:
    """Parameters for the Dorado basecall server connection."""
    address: str
    config:  str
    timeout: int = 60   # per-batch timeout in seconds


def _connect_basecall_client(cfg: BasecallConfig):
    """
    Import ont-pybasecall-client-lib (deferred so the rest of the script
    works without it when basecalling is not needed) and connect to Dorado.
    """
    try:
        from pybasecall_client_lib.pyclient import PyBasecallClient
    except ImportError:
        logger.error(
            "ont-pybasecall-client-lib is required for on-the-fly basecalling. "
            "Install a version that matches your Dorado server: "
            "pip install ont-pybasecall-client-lib==<server_version>"
        )
        sys.exit(1)

    client = PyBasecallClient(
        address=cfg.address,
        config=cfg.config,
        priority=PyBasecallClient.high_priority,
        client_name="NASExperiments_benchmark",
    )
    client.connect()
    logger.info(
        "Connected to Dorado basecall server at %s (config: %s, priority: high)",
        cfg.address, cfg.config,
    )
    return client


def _basecall_raw_batch(
    client,
    raw_batch: List[Tuple[str, np.ndarray, dict]],  # (read_id, signal, slow5_meta)
    cfg: BasecallConfig,
) -> Tuple[Dict[str, str], float]:
    """
    Submit a batch of raw (int16) signals to the Dorado basecall server and
    collect the resulting sequences.

    Parameters
    ----------
    raw_batch : list of (read_id, float32_signal, slow5_meta)
        slow5_meta must contain at minimum: offset, range, digitisation,
        sampling_rate, start_time (start_time may be 0 if absent).
    cfg : BasecallConfig

    Returns
    -------
    (sequences, basecall_s) where sequences maps read_id → basecalled sequence
    (empty string if basecalling failed for that read), and basecall_s is the
    wall-clock time spent waiting for Dorado.
    """
    try:
        from pybasecall_client_lib.helper_functions import package_read
    except ImportError:
        logger.error("ont-pybasecall-client-lib not available.")
        sys.exit(1)

    if not raw_batch:
        return {}, 0.0

    # Build package_read kwargs for each read.
    # The signal stored in Result.basecall_data is float32 (calibrated pA).
    # Dorado expects raw int16 — we pass the int16 directly from slow5_meta
    # where available; otherwise we reverse-scale from the float32.
    packaged = []
    id_order: List[str] = []
    for read_id, signal_f32, meta in raw_batch:
        # Use raw int16 signal if carried in meta; fall back to rounding float32.
        raw_int16 = meta.get("signal_int16")
        if raw_int16 is None:
            offset   = float(meta.get("offset",  0.0))
            scaling  = (float(meta.get("range", 1.0)) /
                        float(meta.get("digitisation", 1.0)))
            if scaling != 0:
                raw_int16 = ((signal_f32 / scaling) - offset).round().astype(np.int16)
            else:
                raw_int16 = signal_f32.round().astype(np.int16)

        start_time = meta.get("start_time")
        pkg = package_read(
            read_id=read_id,
            raw_data=raw_int16,
            daq_offset=float(meta.get("offset", 0.0)),
            daq_scaling=(float(meta.get("range", 1.0)) /
                         float(meta.get("digitisation", 1.0))),
            sampling_rate=float(meta.get("sampling_rate", 4000.0)),
            start_time=int(start_time) if start_time is not None else 0,
        )
        packaged.append(pkg)
        id_order.append(read_id)

    # Submit with up to 3 retries.
    t0 = time.perf_counter()
    for attempt in range(3):
        if client.pass_reads(packaged):
            break
        time.sleep(0.05)
    else:
        logger.warning(
            "Could not submit basecall batch of %d reads after 3 attempts", len(packaged)
        )
        return {rid: "" for rid in id_order}, time.perf_counter() - t0

    # Collect results.
    sequences: Dict[str, str] = {}
    deadline  = time.monotonic() + cfg.timeout
    collected = 0

    while collected < len(packaged):
        if time.monotonic() > deadline:
            logger.warning(
                "Basecall timeout: %d/%d reads still pending",
                len(packaged) - collected, len(packaged),
            )
            break

        results = client.get_completed_reads()
        if not results:
            time.sleep(client.throttle)
            continue

        for res_batch in results:
            for res in res_batch:
                if res.get("sub_tag", 0) > 0:
                    collected += 1
                    continue
                meta_r = res["metadata"]
                rid    = meta_r.get("read_id", "")
                seq    = res["datasets"].get("sequence", "")
                sequences[rid] = seq
                collected += 1

    basecall_s = time.perf_counter() - t0

    # Fill any reads that never came back.
    for rid in id_order:
        sequences.setdefault(rid, "")

    n_ok = sum(1 for s in sequences.values() if s)
    logger.debug("Basecalled %d/%d reads in %.3fs", n_ok, len(packaged), basecall_s)
    return sequences, basecall_s


# ---------------------------------------------------------------------------
# Input readers
# ---------------------------------------------------------------------------

def _iter_blow5(path: str) -> Iterator[Tuple[str, np.ndarray, dict]]:
    """
    Yield (read_id, float32_pA_signal, slow5_meta) from a BLOW5 file.

    slow5_meta carries the fields needed both for calibration and for
    re-packaging reads to the Dorado basecall server:
        signal_int16, offset, range, digitisation, sampling_rate, start_time
    """
    try:
        import pyslow5
    except ImportError:
        logger.error("pyslow5 is required for BLOW5 input: pip install pyslow5")
        sys.exit(1)

    s5 = pyslow5.Open(path, "r")
    for read in s5.seq_reads(aux="all"):
        read_id   = read["read_id"]
        raw_int16 = np.array(read["signal"], dtype=np.int16)

        offset    = read.get("offset",  0.0)
        digitisation = read.get("digitisation", 1.0)
        rng       = read.get("range", 1.0)
        scaling   = float(rng) / float(digitisation) if float(digitisation) != 0 else 1.0
        signal_f32 = (raw_int16.astype(np.float32) + float(offset)) * scaling

        meta = {
            "signal_int16":  raw_int16,        # kept for Dorado re-packaging
            "offset":        offset,
            "range":         rng,
            "digitisation":  digitisation,
            "sampling_rate": read.get("sampling_rate", 4000.0),
            "start_time":    read.get("start_time"),
        }
        yield read_id, signal_f32, meta
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


# ---------------------------------------------------------------------------
# Batch iterator with optional truncation and on-the-fly basecalling
# ---------------------------------------------------------------------------

def make_result_iter(
    paths:             List[Path],
    input_type:        str,
    batch_size:        int,
    truncate_signals:  Optional[int]       = None,
    truncate_bases:    Optional[int]       = None,
    basecall_cfg:      Optional[BasecallConfig] = None,
    basecall_client=None,
) -> Iterator[List[Result]]:
    """
    Yield batches of Result objects from one or more BLOW5 or FASTA/Q files.

    Truncation
    ----------
    truncate_signals : int or None
        Clip raw signal arrays to this many samples before any further
        processing (basecalling or raw-signal alignment).  Applied to BLOW5
        input only.
    truncate_bases : int or None
        Clip sequence strings to this many bases before alignment.
        Applied to FASTA/Q input, and also to basecalled sequences when
        BLOW5 input is used with on-the-fly basecalling.

    On-the-fly basecalling
    ----------------------
    When basecall_client is not None (and input_type == 'blow5'), each batch
    of raw signals is submitted to the Dorado server before being handed to
    the aligner.  The resulting Result objects have seq=<basecalled_sequence>
    and basecall_data=None, making them compatible with sequence-based aligner
    plugins.
    """
    channel = 0
    batch: List[Result] = []

    # Convenience flag: are we basecalling on the fly?
    do_basecall = (basecall_client is not None) and (input_type == "blow5")

    def _flush_basecall_buffer(buf: list) -> Tuple[List[Result], float]:
        """Basecall a buffer of (read_id, signal, meta) and return (Results, basecall_s)."""
        if not buf:
            return [], 0.0
        sequences, basecall_s = _basecall_raw_batch(basecall_client, buf, basecall_cfg)
        out: List[Result] = []
        for rid, sig_f32, _meta in buf:
            seq = sequences.get(rid, "") or "N"
            if truncate_bases is not None:
                seq = seq[:truncate_bases]
            out.append(Result(
                channel=channel % 512,
                read_id=rid,
                seq=seq or "N",
                basecall_data=None,   # sequence aligner — no signal needed
            ))
        return out, basecall_s

    for path in paths:
        logger.info("Reading %s", path)

        if input_type == "blow5":
            if do_basecall:
                # ── Prefetch: overlap basecalling batch N+1 with alignment of batch N ──
                # We read from the BLOW5 file into a staging buffer. As soon as the
                # buffer is full (batch_size reads), we fire off a background thread
                # to basecall it while the caller's aligner is busy with the previous
                # batch. The next iteration simply waits on the future before yielding.
                with ThreadPoolExecutor(max_workers=1) as pool:
                    staging: List[Tuple[str, np.ndarray, dict]] = []
                    prefetch_future: Optional[Future] = None

                    def _submit_staging(buf):
                        """Copy buf and dispatch basecalling in the thread pool."""
                        return pool.submit(_flush_basecall_buffer, list(buf))

                    for read_id, signal_f32, meta in _iter_blow5(str(path)):
                        if truncate_signals is not None:
                            signal_f32 = signal_f32[:truncate_signals]
                            if meta.get("signal_int16") is not None:
                                meta = dict(meta)
                                meta["signal_int16"] = meta["signal_int16"][:truncate_signals]

                        staging.append((read_id, signal_f32, meta))

                        if len(staging) >= batch_size:
                            # Collect the previous prefetch (if any) before firing the next.
                            if prefetch_future is not None:
                                bc_results, basecall_s = prefetch_future.result()
                                for r in bc_results:
                                    batch.append(r)
                                    channel += 1
                                if batch:
                                    yield batch, basecall_s
                                    batch = []

                            prefetch_future = _submit_staging(staging)
                            staging.clear()

                    # Collect the last in-flight future.
                    if prefetch_future is not None:
                        bc_results, basecall_s = prefetch_future.result()
                        for r in bc_results:
                            batch.append(r)
                            channel += 1

                    # Basecall any leftover reads that didn't fill a full buffer.
                    if staging:
                        bc_results, basecall_s = _flush_basecall_buffer(staging)
                        for r in bc_results:
                            batch.append(r)
                            channel += 1
                        staging.clear()

                    if batch:
                        yield batch, basecall_s
                        batch = []

            else:
                # Raw-signal aligner path — no basecalling, no prefetch needed.
                for read_id, signal_f32, meta in _iter_blow5(str(path)):
                    if truncate_signals is not None:
                        signal_f32 = signal_f32[:truncate_signals]
                        if meta.get("signal_int16") is not None:
                            meta = dict(meta)
                            meta["signal_int16"] = meta["signal_int16"][:truncate_signals]

                    batch.append(Result(
                        channel=channel % 512,
                        read_id=read_id,
                        seq="N",
                        basecall_data=signal_f32,
                    ))
                    channel += 1
                    if len(batch) >= batch_size:
                        yield batch, 0.0
                        batch = []

        else:  # FASTA/Q
            for read_id, seq in _iter_fasta(str(path)):
                if truncate_bases is not None:
                    seq = seq[:truncate_bases]
                batch.append(Result(
                    channel=channel % 512,
                    read_id=read_id,
                    seq=seq or "N",
                    basecall_data=None,
                ))
                channel += 1
                if len(batch) >= batch_size:
                    yield batch, 0.0
                    batch = []

    if batch:
        yield batch, 0.0


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
    if input_type == "blow5" and result.basecall_data is not None:
        return len(result.basecall_data)
    seq = result.seq or ""
    return len(seq) if seq != "N" else 0


class DecisionTSVWriter:
    """
    Writes per-read decisions in the exact readfish_decisions.tsv column order
    so that analyse_runs.py can consume standalone benchmark output directly.
    """

    def __init__(self, path: Path):
        path.parent.mkdir(parents=True, exist_ok=True)
        self._fh = open(path, "w", newline="", buffering=1)
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
        result_times:     Dict[str, float],
    ) -> None:
        for read_in_loop, result in enumerate(results, start=1):
            decision = _result_to_decision(result, input_type)
            self._writer.writerow({
                "client_iteration": client_iteration,
                "read_in_loop":     read_in_loop,
                "read_id":          result.read_id,
                "channel":          result.channel,
                "seq_len":          _seq_len_for_tsv(result, input_type),
                "counter":          1,
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
    input_paths:      List[Path],
    input_type:       str,
    batch_size:       int,
    manifest:         Optional[Dict[str, int]],
    outdir:           Path,
    tool:             str,
    dataset:          str,
    decisions_tsv:    Optional[Path]          = None,  # defaults to outdir/readfish_decisions.tsv
    truncate_signals: Optional[int]           = None,
    truncate_bases:   Optional[int]           = None,
    basecall_cfg:     Optional[BasecallConfig] = None,
) -> None:
    outdir.mkdir(parents=True, exist_ok=True)

    per_read_path  = outdir / f"{tool}_{dataset}_per_read.tsv"
    batch_log_path = outdir / f"{tool}_{dataset}_batch_stats.tsv"

    per_read_rows:  List[dict] = []
    batch_rows:     List[dict] = []

    total_reads  = 0
    total_mapped = 0

    # Connect to the Dorado basecall server if requested.
    basecall_client = None
    if basecall_cfg is not None:
        if input_type != "blow5":
            logger.error(
                "On-the-fly basecalling (--basecall-address / --basecall-config) "
                "is only supported with BLOW5 input, not FASTA/Q."
            )
            sys.exit(1)
        basecall_client = _connect_basecall_client(basecall_cfg)
        logger.info(
            "On-the-fly basecalling enabled: signals will be sent to Dorado "
            "before querying the aligner."
        )
        if truncate_signals is not None:
            logger.info(
                "Signal truncation: keeping first %d samples before basecalling "
                "(≈ %.0f bases at 4 kHz / 450 bps).",
                truncate_signals,
                truncate_signals / 4000 * 450,
            )
        if truncate_bases is not None:
            logger.info(
                "Sequence truncation: keeping first %d bases after basecalling.",
                truncate_bases,
            )
    else:
        if truncate_signals is not None:
            if input_type != "blow5":
                logger.error("--truncate-signals requires BLOW5 input.")
                sys.exit(1)
            logger.info(
                "Signal truncation: keeping first %d samples before alignment "
                "(≈ %.0f bases at 4 kHz / 450 bps).",
                truncate_signals,
                truncate_signals / 4000 * 450,
            )
        if truncate_bases is not None:
            if input_type != "fasta":
                logger.error(
                    "--truncate-bases requires FASTA/Q input (or BLOW5 with "
                    "--basecall-address / --basecall-config for on-the-fly basecalling)."
                )
                sys.exit(1)
            logger.info(
                "Sequence truncation: keeping first %d bases before alignment.",
                truncate_bases,
            )

    decisions_tsv = decisions_tsv or outdir / "readfish_decisions.tsv"
    dec_writer = DecisionTSVWriter(decisions_tsv)
    logger.info("Writing decisions TSV → %s", decisions_tsv)

    # Determine the effective input_type seen by the aligner.
    # When basecalling on the fly, reads arrive at the aligner as sequences.
    effective_input_type = "fasta" if basecall_client is not None else input_type

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
                "basecall_s", "align_s", "total_s",
                "throughput_reads_per_s",
                "mean_latency_s", "median_latency_s",
                "p95_latency_s", "p99_latency_s",
            ], delimiter="\t")
            bl_writer.writeheader()

            for batch_idx, (batch, basecall_s) in enumerate(
                    make_result_iter(
                        paths=input_paths,
                        input_type=input_type,
                        batch_size=batch_size,
                        truncate_signals=truncate_signals,
                        truncate_bases=truncate_bases,
                        basecall_cfg=basecall_cfg,
                        basecall_client=basecall_client,
                    )):

                n_reads = len(batch)
                t_batch_start   = time.perf_counter()
                wall_time_start = time.time()

                start_times:        Dict[str, float] = {}
                result_timestamps:  Dict[str, float] = {}

                def _timed_iter(results):
                    for r in results:
                        start_times[r.read_id] = time.perf_counter()
                        yield r

                results = list(aligner.map_reads(_timed_iter(batch)))

                t_batch_end = time.perf_counter()
                # align_s is the wall time spent purely inside map_reads().
                # basecall_s was already measured upstream (0.0 if not basecalling).
                # When prefetch overlap is active, basecall_s for batch N was
                # spent concurrently with alignment of batch N-1, so total_s
                # reflects the true wall-clock cost of this batch.

                n_results = len(results)
                for i, result in enumerate(results):
                    frac = (i + 1) / n_results if n_results > 0 else 1.0
                    result_timestamps[result.read_id] = (
                        wall_time_start + frac * (time.time() - wall_time_start)
                    )

                read_times: List[float] = []
                read_stats: List[ReadStat] = []

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
                    seq_len    = _seq_len_for_tsv(result, effective_input_type)

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

                for rs in read_stats:
                    pr_writer.writerow(asdict(rs))

                dec_writer.write_batch(
                        results=results,
                        input_type=effective_input_type,
                        client_iteration=batch_idx + 1,
                        batch_start_time=wall_time_start,
                        result_times=result_timestamps,
                    )

                n_mapped   = sum(1 for rs in read_stats if rs.mapped)
                lat        = np.array(read_times)
                align_s    = t_batch_end - t_batch_start
                total_s    = basecall_s + align_s
                throughput = n_reads / total_s if total_s > 0 else float("nan")

                bl_row = {
                    "batch_idx":              batch_idx,
                    "n_reads":                n_reads,
                    "n_mapped":               n_mapped,
                    "basecall_s":             round(basecall_s, 6),
                    "align_s":                round(align_s,    6),
                    "total_s":                round(total_s,    6),
                    "throughput_reads_per_s": round(throughput, 2),
                    "mean_latency_s":         round(float(np.mean(lat)),           6),
                    "median_latency_s":       round(float(np.median(lat)),         6),
                    "p95_latency_s":          round(float(np.percentile(lat, 95)), 6),
                    "p99_latency_s":          round(float(np.percentile(lat, 99)), 6),
                }
                bl_writer.writerow(bl_row)

                total_reads  += n_reads
                total_mapped += n_mapped

                if basecall_s > 0:
                    logger.info(
                        "Batch %4d | %4d reads | %4d mapped | "
                        "bc %.3fs  aln %.3fs  tot %.3fs | %.1f reads/s",
                        batch_idx, n_reads, n_mapped,
                        basecall_s, align_s, total_s, throughput,
                    )
                else:
                    logger.info(
                        "Batch %4d | %4d reads | %4d mapped | %.3fs | %.1f reads/s",
                        batch_idx, n_reads, n_mapped, align_s, throughput,
                    )

    finally:
        dec_writer.close()
        if basecall_client is not None:
            try:
                basecall_client.disconnect()
                logger.info("Disconnected from Dorado basecall server.")
            except Exception:
                pass

    summary = {
        "tool":             tool,
        "dataset":          dataset,
        "input_type":       input_type,
        "effective_input":  effective_input_type,
        "truncate_signals": truncate_signals,
        "truncate_bases":   truncate_bases,
        "total_reads":      total_reads,
        "total_mapped":     total_mapped,
        "map_rate":         round(total_mapped / total_reads, 4) if total_reads else 0.0,
    }
    summary_path = outdir / f"{tool}_{dataset}_summary.json"
    with open(summary_path, "w") as fh:
        json.dump(summary, fh, indent=2)

    logger.info("Done. %d reads, %d mapped (%.1f%%)",
                total_reads, total_mapped, 100 * summary["map_rate"])
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

    # ── Input / output ─────────────────────────────────────────────────────
    p.add_argument("--input",       required=True, action="append", dest="inputs",
                   metavar="PATH_OR_GLOB",
                   help="Input file, folder, or glob pattern. "
                        "Repeat for multiple sources. "
                        "All inputs must be the same type (BLOW5 or FASTA/Q).")
    p.add_argument("--output",      required=True,
                   help="Output directory for TSV and JSON results")
    p.add_argument("--tool",        default=None,
                   help="Tool label for output filenames (default: --plugin value)")
    p.add_argument("--dataset",     default="dataset",
                   help="Dataset label for output filenames (default: dataset)")
    p.add_argument("--manifest",    default=None,
                   help="Optional TSV mapping contig/species → true_label")


    # ── Aligner plugin ─────────────────────────────────────────────────────
    p.add_argument("--plugin",      required=True,
                   help="Python module name of the aligner plugin, e.g. pyrawhash, mappy")
    p.add_argument("--plugin-args", default="",
                   help='Space-separated key=value args passed to Aligner(), '
                        'e.g. "fn_idx_in=/ref.fa n_threads=16"')
    p.add_argument("--batch-size",  type=int, default=512,
                   help="Number of reads per batch (default: 512)")

    # ── On-the-fly basecalling ──────────────────────────────────────────────
    bc = p.add_argument_group(
        "on-the-fly basecalling",
        "Basecall raw BLOW5 signals before sending to a sequence aligner. "
        "Requires ont-pybasecall-client-lib (version must match Dorado server).",
    )
    bc.add_argument("--basecall-address", default=None,
                    help="Dorado basecall server address "
                         "(e.g. ipc:///var/lib/minknow/data/.dorado/"
                         "dorado-basecall-server.sock  or  127.0.0.1:5555). "
                         "Must also set --basecall-config.")
    bc.add_argument("--basecall-config",  default=None,
                    help="Dorado model config string "
                         "(e.g. dna_r10.4.1_e8.2_400bps_fast@v5.2.0||). "
                         "Must also set --basecall-address.")
    bc.add_argument("--basecall-timeout", type=int, default=60,
                    help="Per-batch basecall timeout in seconds (default: 60). "
                         "If the server has not returned all reads within this "
                         "wall-clock budget, remaining reads are skipped.")

    # ── Truncation ─────────────────────────────────────────────────────────
    trunc = p.add_argument_group(
        "truncation",
        "Clip reads to a fixed length to simulate MinKNOW chunk cadence "
        "(0.4–0.8 s ≈ 1 600–3 200 raw samples ≈ 180–360 bases). "
        "At most one of the two flags may be set.",
    )
    trunc.add_argument("--truncate-signals", type=int, default=None,
                       metavar="N",
                       help="Keep only the first N raw signal samples. "
                            "Applied to BLOW5 input before basecalling "
                            "(if --basecall-address is set) or before raw-signal "
                            "alignment. Not compatible with FASTA/Q input.")
    trunc.add_argument("--truncate-bases",   type=int, default=None,
                       metavar="N",
                       help="Keep only the first N sequence bases. "
                            "Applied to FASTA/Q input, or to basecalled sequences "
                            "when BLOW5 input is used with --basecall-address.")

    # ── Misc ───────────────────────────────────────────────────────────────
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

    # Validate truncation flags.
    if args.truncate_signals is not None and args.truncate_bases is not None:
        parser.error("--truncate-signals and --truncate-bases are mutually exclusive.")
    if args.truncate_signals is not None and args.truncate_signals <= 0:
        parser.error("--truncate-signals must be a positive integer.")
    if args.truncate_bases is not None and args.truncate_bases <= 0:
        parser.error("--truncate-bases must be a positive integer.")

    # Validate basecalling flags.
    if bool(args.basecall_address) != bool(args.basecall_config):
        parser.error(
            "--basecall-address and --basecall-config must be set together."
        )

    # Resolve inputs.
    try:
        input_paths, input_type = resolve_inputs(args.inputs)
    except (ValueError, FileNotFoundError) as e:
        logger.error("%s", e)
        sys.exit(1)

    # Further cross-validation now that input_type is known.
    if args.truncate_signals is not None and input_type != "blow5":
        parser.error("--truncate-signals requires BLOW5 input.")
    if args.truncate_bases is not None and input_type == "blow5" \
            and not args.basecall_address:
        parser.error(
            "--truncate-bases with BLOW5 input requires --basecall-address "
            "and --basecall-config for on-the-fly basecalling."
        )

    tool     = args.tool or args.plugin
    outdir   = Path(args.output)
    manifest = load_manifest(args.manifest) if args.manifest else None
    aligner  = load_plugin(args.plugin, args.plugin_args,
                           getattr(args, "debug_log", None))

    basecall_cfg = None
    if args.basecall_address:
        basecall_cfg = BasecallConfig(
            address=args.basecall_address,
            config=args.basecall_config,
            timeout=args.basecall_timeout,
        )

    run_benchmark(
        aligner=aligner,
        input_paths=input_paths,
        input_type=input_type,
        batch_size=args.batch_size,
        manifest=manifest,
        outdir=outdir,
        tool=tool,
        dataset=args.dataset,

        truncate_signals=args.truncate_signals,
        truncate_bases=args.truncate_bases,
        basecall_cfg=basecall_cfg,
    )


if __name__ == "__main__":
    main()