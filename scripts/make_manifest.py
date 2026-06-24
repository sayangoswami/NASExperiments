#!/usr/bin/env python3
"""
make_manifest.py
================
Generate a ground-truth manifest TSV from two community FASTA directories.

The manifest maps every contig name to its community label:
    contig        true_label
    Rothia_dentocariosa_contig_1    0
    Neisseria_subflava_contig_1     0
    Prevotella_nigrescens_contig_1  1
    ...

This manifest is consumed by analyse_runs.py --manifest to resolve ground
truth per-read from the contig name encoded in the read ID, enabling
mixed-community runs to be evaluated without separate community replays.

Usage
-----
    python make_manifest.py \\
        --community0  /data/gut/community0/ \\
        --community1  /data/gut/community1/ \\
        --output      /data/gut/gut_manifest.tsv

    # Single FASTA files are also accepted:
    python make_manifest.py \\
        --community0  /data/zymo/community0.fasta \\
        --community1  /data/zymo/community1.fasta \\
        --output      /data/zymo/zymo_manifest.tsv
"""

import argparse
import csv
import sys
from pathlib import Path


FASTA_EXTENSIONS = {".fasta", ".fa", ".fna"}


def iter_contig_names(path: str):
    """Yield contig names (first word of FASTA header) from a file or directory."""
    p = Path(path)
    files = []

    if p.is_file():
        files = [p]
    elif p.is_dir():
        files = sorted(f for f in p.rglob("*")
                       if f.suffix in FASTA_EXTENSIONS and f.is_file())
        if not files:
            print(f"WARNING: No FASTA files found in {path}", file=sys.stderr)
    else:
        print(f"ERROR: {path} is not a file or directory", file=sys.stderr)
        sys.exit(1)

    for fasta in files:
        with open(fasta) as fh:
            for line in fh:
                if line.startswith(">"):
                    contig = line[1:].split()[0].strip()
                    yield contig


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--community0", required=True,
                   help="Directory of FASTA files or single FASTA for "
                        "Community 0 (deplete / true negative)")
    p.add_argument("--community1", required=True,
                   help="Directory of FASTA files or single FASTA for "
                        "Community 1 (target / true positive)")
    p.add_argument("--output", required=True,
                   help="Output manifest TSV path")
    args = p.parse_args()

    rows = []
    seen = {}

    for label, path in [(0, args.community0), (1, args.community1)]:
        for contig in iter_contig_names(path):
            if contig in seen:
                print(
                    f"WARNING: contig '{contig}' appears in both communities "
                    f"(first seen in community {seen[contig]}, "
                    f"now in community {label}) — keeping first.",
                    file=sys.stderr,
                )
                continue
            seen[contig] = label
            rows.append({"contig": contig, "true_label": label})

    if not rows:
        print("ERROR: No contigs found — check your input paths.", file=sys.stderr)
        sys.exit(1)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=["contig", "true_label"], delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    n0 = sum(1 for r in rows if r["true_label"] == 0)
    n1 = sum(1 for r in rows if r["true_label"] == 1)
    print(f"Wrote {len(rows)} contigs to {args.output}")
    print(f"  Community 0 (deplete): {n0} contigs")
    print(f"  Community 1 (target):  {n1} contigs")


if __name__ == "__main__":
    main()