#!/usr/bin/env python3
"""
sample_reads.py
===============
Sample exact substrings (error-free reads) from a collection of reference
FASTA files and write them as a FASTA suitable for seq2squiggle --read-input.

Designed for large metagenomic communities (e.g. 1,383-species gut metagenome)
where seq2squiggle cannot process the full reference in one pass due to runtime
constraints. This script:

  1. Streams each reference FASTA one species at a time — never loads the full
     community into memory simultaneously.
  2. Samples reads at a fixed coverage depth per species (genome_size × coverage
     / mean_read_length reads per species).
  3. Draws read lengths from a user-specified length distribution (exponential
     by default, matching real ONT data; or from a TSV of empirical lengths).
  4. Writes encoded read IDs: <contig>!<start>!<end>!<strand>
     compatible with analyse_runs.py and qc_signals.sh.
  5. Optionally splits output into chunks of N reads so seq2squiggle can be
     run in parallel on each chunk.

Usage
-----
    # Exponential length distribution (mean 1500 bp, min 200 bp)
    python3 sample_reads.py \\
        --input    community1/         \\   # directory of per-species FASTAs
        --output   community1_reads/   \\   # output directory for FASTA chunks
        --coverage 30                  \\   # target coverage per species
        --length-dist exponential      \\
        --mean-length 1500             \\
        --min-length  200              \\
        --chunk-size  50000            \\   # reads per output FASTA chunk
        --seed        42

    # Empirical length distribution from a TSV (one length per line, no header)
    python3 sample_reads.py \\
        --input    community1/         \\
        --output   community1_reads/   \\
        --coverage 30                  \\
        --length-dist empirical        \\
        --lengths-tsv zymo_read_lengths.tsv \\
        --chunk-size  50000            \\
        --seed        42

Then run seq2squiggle on each chunk (see run_seq2squiggle.sh).

Output
------
    community1_reads/
        chunk_0000.fasta    reads 0–49999
        chunk_0001.fasta    reads 50000–99999
        ...
        manifest.tsv        chunk_id, n_reads, fasta_path
        sampling_stats.tsv  per-species: genome_size, n_reads_sampled,
                            actual_coverage, skipped_contigs
"""

from __future__ import annotations

import argparse
import csv
import logging
import math
import os
import random
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Generator, List, Optional, Tuple

import numpy as np

logger = logging.getLogger("sample_reads")

# ---------------------------------------------------------------------------
# Length distribution samplers
# ---------------------------------------------------------------------------

class ExponentialSampler:
    """
    Samples read lengths from an exponential distribution, matching the
    heavy-tailed length distribution of real ONT data.
    Lengths are drawn as: min_length + Exponential(mean - min_length)
    then capped at max_length.
    """
    def __init__(self, mean: int, min_length: int, max_length: int, rng: random.Random):
        self.min_length  = min_length
        self.max_length  = max_length
        self.scale       = max(mean - min_length, 1)
        self.rng         = rng

    def sample(self) -> int:
        length = int(self.min_length + self.rng.expovariate(1.0 / self.scale))
        return min(length, self.max_length)

    @property
    def mean(self) -> float:
        return self.min_length + self.scale


class EmpiricalSampler:
    """
    Samples read lengths by drawing uniformly from a list of observed lengths
    (e.g. extracted from a real or Zymo run via NanoStat or seqkit).

    TSV format: one integer per line, no header.
    """
    def __init__(self, lengths_tsv: str, min_length: int, rng: random.Random):
        lengths = []
        with open(lengths_tsv) as fh:
            for line in fh:
                line = line.strip()
                if line and line.isdigit():
                    l = int(line)
                    if l >= min_length:
                        lengths.append(l)
        if not lengths:
            raise ValueError(f"No valid lengths found in {lengths_tsv} "
                             f"(min_length={min_length})")
        self.lengths    = lengths
        self.rng        = rng
        logger.info("Empirical sampler: %d lengths loaded from %s "
                    "(mean=%.0f, N50 approx)", len(lengths), lengths_tsv,
                    sum(lengths) / len(lengths))

    def sample(self) -> int:
        return self.rng.choice(self.lengths)

    @property
    def mean(self) -> float:
        return sum(self.lengths) / len(self.lengths)


# ---------------------------------------------------------------------------
# FASTA streaming
# ---------------------------------------------------------------------------

def stream_fasta(fasta_path: str) -> Generator[Tuple[str, str], None, None]:
    """
    Stream sequences from a FASTA file one contig at a time.
    Yields (header, sequence) pairs. Does not load the whole file at once.
    """
    header = None
    chunks = []
    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].split()[0]   # use first word of header as name
                chunks = []
            else:
                chunks.append(line.upper())
    if header is not None:
        yield header, "".join(chunks)


# ---------------------------------------------------------------------------
# Complement table for reverse complement
# ---------------------------------------------------------------------------

_COMP = str.maketrans("ACGTN", "TGCAN")

def reverse_complement(seq: str) -> str:
    return seq.translate(_COMP)[::-1]


# ---------------------------------------------------------------------------
# Per-species read sampler
# ---------------------------------------------------------------------------

@dataclass
class SamplingStats:
    fasta_path:       str
    species:          str
    genome_size:      int   = 0
    n_contigs:        int   = 0
    n_contigs_skip:   int   = 0   # contigs too short for even one read
    n_reads_target:   int   = 0
    n_reads_sampled:  int   = 0
    actual_coverage:  float = 0.0


def sample_species(
    fasta_path: str,
    coverage: float,
    length_sampler,
    min_length: int,
    rng: random.Random,
    max_retries: int = 20,
) -> Tuple[List[Tuple[str, str]], SamplingStats]:
    """
    Sample reads from all contigs of one species at a fixed coverage depth.

    Coverage is computed per-contig proportionally — each contig contributes
    reads proportional to its length relative to the total genome size. This
    prevents large contigs from dominating while small contigs get no reads.

    Returns (reads, stats) where reads is a list of (read_id, sequence).
    read_id format: <contig_name>!<start>!<end>!<strand>
    """
    species = Path(fasta_path).stem
    stats   = SamplingStats(fasta_path=fasta_path, species=species)

    # First pass: collect contig names and lengths (no full sequence kept)
    contigs = []
    for name, seq in stream_fasta(fasta_path):
        contigs.append((name, len(seq)))
        stats.genome_size += len(seq)
        stats.n_contigs   += 1

    if stats.genome_size == 0:
        logger.warning("Empty genome for %s, skipping", species)
        return [], stats

    mean_len = length_sampler.mean
    stats.n_reads_target = max(1, int(stats.genome_size * coverage / mean_len))

    # Distribute reads across contigs proportional to length
    contig_read_counts = {}
    total_assigned = 0
    for name, length in contigs:
        if length < min_length:
            stats.n_contigs_skip += 1
            continue
        n = max(1, round(stats.n_reads_target * length / stats.genome_size))
        contig_read_counts[name] = n
        total_assigned += n

    if not contig_read_counts:
        logger.warning("All contigs in %s are shorter than min_length=%d",
                       species, min_length)
        return [], stats

    # Second pass: sample reads from each contig
    reads = []
    for name, seq in stream_fasta(fasta_path):
        n_reads = contig_read_counts.get(name, 0)
        if n_reads == 0:
            continue

        seq_len = len(seq)
        sampled = 0
        retries = 0

        while sampled < n_reads and retries < n_reads * max_retries:
            read_len = length_sampler.sample()
            if read_len > seq_len:
                # Read longer than contig — use whole contig
                read_len = seq_len

            start  = rng.randint(0, seq_len - read_len)
            end    = start + read_len
            strand = rng.choice(["+", "-"])
            subseq = seq[start:end]

            # Skip reads with too many Ns
            n_frac = subseq.count("N") / max(len(subseq), 1)
            if n_frac > 0.1:
                retries += 1
                continue

            if strand == "-":
                subseq = reverse_complement(subseq)

            read_id = f"{name}!{start}!{end}!{strand}"
            reads.append((read_id, subseq))
            sampled += 1

        stats.n_reads_sampled += sampled

    stats.actual_coverage = (
        stats.n_reads_sampled * mean_len / stats.genome_size
        if stats.genome_size > 0 else 0.0
    )

    return reads, stats


# ---------------------------------------------------------------------------
# Chunked FASTA writer
# ---------------------------------------------------------------------------

class ChunkedFastaWriter:
    """
    Writes reads to a series of FASTA files, each containing at most
    chunk_size reads. Maintains a manifest TSV.
    """
    def __init__(self, outdir: Path, chunk_size: int):
        self.outdir      = outdir
        self.chunk_size  = chunk_size
        self.chunk_id    = 0
        self.n_in_chunk  = 0
        self.total       = 0
        self._fh         = None
        self._chunk_path = None
        self.manifest    = []   # (chunk_id, n_reads, path)
        outdir.mkdir(parents=True, exist_ok=True)
        self._open_chunk()

    def _open_chunk(self):
        if self._fh:
            self._fh.close()
            self.manifest.append((self.chunk_id, self.n_in_chunk,
                                  str(self._chunk_path)))
            self.chunk_id += 1
            self.n_in_chunk = 0
        self._chunk_path = self.outdir / f"chunk_{self.chunk_id:04d}.fasta"
        self._fh = open(self._chunk_path, "w")

    def write(self, read_id: str, sequence: str):
        if self.n_in_chunk >= self.chunk_size:
            self._open_chunk()
        self._fh.write(f">{read_id}\n{sequence}\n")
        self.n_in_chunk += 1
        self.total += 1

    def close(self):
        if self._fh:
            self._fh.close()
            if self.n_in_chunk > 0:
                self.manifest.append((self.chunk_id, self.n_in_chunk,
                                      str(self._chunk_path)))
            self._fh = None

    def write_manifest(self, manifest_path: Path):
        with open(manifest_path, "w", newline="") as f:
            w = csv.writer(f, delimiter="\t")
            w.writerow(["chunk_id", "n_reads", "fasta_path"])
            w.writerows(self.manifest)


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--input",   required=True,
                   help="Directory of per-species FASTA files, or a single "
                        "multi-species FASTA. If a directory, every *.fasta / "
                        "*.fa / *.fna file is treated as one species.")
    p.add_argument("--output",  required=True,
                   help="Output directory for chunked FASTA files")
    p.add_argument("--coverage", type=float, default=30.0,
                   help="Target coverage depth per species (default: 30)")
    p.add_argument("--min-length", type=int, default=200,
                   help="Minimum read length in bp (default: 200)")
    p.add_argument("--max-length", type=int, default=100_000,
                   help="Maximum read length in bp (default: 100000)")
    p.add_argument("--chunk-size", type=int, default=50_000,
                   help="Max reads per output FASTA chunk (default: 50000). "
                        "Tune so each chunk finishes in seq2squiggle within "
                        "your time limit.")
    p.add_argument("--seed", type=int, default=42,
                   help="Random seed for reproducibility (default: 42)")

    dist = p.add_argument_group("length distribution (choose one)")
    dist.add_argument("--length-dist",
                      choices=["exponential", "empirical"],
                      default="exponential")
    dist.add_argument("--mean-length", type=int, default=1500,
                      help="Mean read length for exponential distribution "
                           "(default: 1500)")
    dist.add_argument("--lengths-tsv",
                      help="TSV of empirical read lengths (one integer per "
                           "line, no header). Required when --length-dist "
                           "empirical.")

    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    return p.parse_args()


def collect_fasta_files(input_path: str) -> List[str]:
    p = Path(input_path)
    if p.is_file():
        return [str(p)]
    exts = {".fasta", ".fa", ".fna"}
    files = sorted(f for f in p.rglob("*") if f.suffix in exts and f.is_file())
    if not files:
        raise ValueError(f"No FASTA files found in {input_path}")
    return [str(f) for f in files]


def main():
    args = parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s  %(message)s",
        datefmt="%H:%M:%S",
    )

    rng = random.Random(args.seed)

    # Build length sampler
    if args.length_dist == "exponential":
        sampler = ExponentialSampler(
            mean=args.mean_length,
            min_length=args.min_length,
            max_length=args.max_length,
            rng=rng,
        )
        logger.info("Length distribution: exponential (mean=%d, min=%d)",
                    args.mean_length, args.min_length)
    else:
        if not args.lengths_tsv:
            sys.exit("ERROR: --lengths-tsv is required with --length-dist empirical")
        sampler = EmpiricalSampler(args.lengths_tsv, args.min_length, rng)
        logger.info("Length distribution: empirical (mean=%.0f)", sampler.mean)

    # Collect input files
    fasta_files = collect_fasta_files(args.input)
    logger.info("Found %d FASTA file(s) to process", len(fasta_files))

    outdir = Path(args.output)
    writer = ChunkedFastaWriter(outdir, args.chunk_size)
    all_stats: List[SamplingStats] = []

    for i, fasta_path in enumerate(fasta_files):
        species = Path(fasta_path).stem
        logger.info("[%d/%d] Sampling %s ...", i + 1, len(fasta_files),
                    species)
        reads, stats = sample_species(
            fasta_path=fasta_path,
            coverage=args.coverage,
            length_sampler=sampler,
            min_length=args.min_length,
            rng=rng,
        )
        for read_id, seq in reads:
            writer.write(read_id, seq)
        all_stats.append(stats)
        logger.debug(
            "  genome=%d bp, target=%d reads, sampled=%d, "
            "actual_cov=%.1fx, skipped_contigs=%d",
            stats.genome_size, stats.n_reads_target,
            stats.n_reads_sampled, stats.actual_coverage,
            stats.n_contigs_skip,
        )

    writer.close()
    writer.write_manifest(outdir / "manifest.tsv")

    # Write per-species stats
    stats_path = outdir / "sampling_stats.tsv"
    with open(stats_path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["species", "fasta_path", "genome_size", "n_contigs",
                    "n_contigs_skipped", "n_reads_target", "n_reads_sampled",
                    "actual_coverage"])
        for s in all_stats:
            w.writerow([s.species, s.fasta_path, s.genome_size, s.n_contigs,
                        s.n_contigs_skip, s.n_reads_target, s.n_reads_sampled,
                        f"{s.actual_coverage:.2f}"])

    # Summary
    total_reads   = writer.total
    total_chunks  = len(writer.manifest)
    total_genome  = sum(s.genome_size for s in all_stats)
    skipped_contigs = sum(s.n_contigs_skip for s in all_stats)

    print()
    print("=" * 50)
    print(f"  Species processed : {len(all_stats)}")
    print(f"  Total genome size : {total_genome / 1e9:.2f} Gb")
    print(f"  Total reads       : {total_reads:,}")
    print(f"  Output chunks     : {total_chunks}")
    print(f"  Skipped contigs   : {skipped_contigs} (shorter than {args.min_length} bp)")
    print(f"  Output directory  : {outdir}")
    print(f"  Sampling stats    : {stats_path}")
    print(f"  Chunk manifest    : {outdir / 'manifest.tsv'}")
    print("=" * 50)


if __name__ == "__main__":
    main()