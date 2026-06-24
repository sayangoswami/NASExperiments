#!/usr/bin/env python3
"""
analyse_runs.py
===============
Post-run analysis for the Selective Sequencing benchmark.

Ground truth — three modes (mutually exclusive, checked in this order)
-----------------------------------------------------------------------
1. --manifest <tsv>   A TSV with columns: contig, true_label (0=deplete, 1=target).
                      Generated once from your community FASTA lists (see
                      make_manifest.py). Works for mixed-community runs and for
                      separate community runs alike.

2. --community <0|1>  The original mode: the whole run is one community, so
                      every read is a true positive (1) or true negative (0)
                      by construction. No manifest needed.

3. Read-ID lookup     If neither --manifest nor --community is provided the
                      script attempts to derive truth from the contig name in
                      the encoded read ID by matching against
                      --community1-contigs / --community0-contigs lists.
                      This is a fallback and is less reliable than a manifest.

Read ID formats supported
--------------------------
  seq2squiggle (4 fields): <contig>!<start>!<end>!<strand>
  squigulator   (5 fields): <prefix>!<contig>!<start>!<end>!<strand>
  In both cases parts[-4] is the contig name.

Inputs per run (one run = one community replayed through one tool)
-----------------------------------------------------------------
1. readfish_decisions.tsv  -- the per-read decision TSV that Readfish writes
   Columns (tab-separated, with header):
     client_iteration  read_in_loop  read_id  channel  seq_len  counter
     mode  decision  condition  barcode  previous_action  action_override
     timestamp

2. unblocked_read_ids.txt  -- one read_id per line for every ejected read
   (used as a cross-check / completeness check; not strictly required if the
   TSV is complete)

Exclusions
----------
Reads where action_override == True are excluded from all accuracy metrics.

Directory convention (batch mode)
----------------------------------
Separate communities:
    <results_dir>/<tool>/community{0,1}/readfish_decisions.tsv

Mixed community (single TSV per tool):
    <results_dir>/<tool>/readfish_decisions.tsv

Usage
-----
# Separate community runs (original mode)
python analyse_runs.py single \
    --tsv       results/minimap2/community1/readfish_decisions.tsv \
    --tool      minimap2 --dataset zymo --community 1 \
    --outdir    results/minimap2/community1/

# Mixed-community run with manifest
python analyse_runs.py single \
    --tsv       results/minimap2/readfish_decisions.tsv \
    --tool      minimap2 --dataset zymo \
    --manifest  data/zymo_manifest.tsv \
    --outdir    results/minimap2/ \
    --unblocked /tmp/MinknoApiSimulator/out/unblocked_read_ids.txt # this is optional but can be used to cross-check the TSV

# Batch — separate communities
python analyse_runs.py batch \
    --results-dir results/zymo/ --dataset zymo --outdir figures/zymo/

# Batch — mixed community with manifest
python analyse_runs.py batch \
    --results-dir results/zymo/ --dataset zymo \
    --manifest    data/zymo_manifest.tsv \
        --outdir      figures/zymo/
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

logger = logging.getLogger("analyse_runs")

# ---------------------------------------------------------------------------
# Readfish decision semantics
# ---------------------------------------------------------------------------

# Decisions that mean the read was accepted (sequenced fully / kept)
ACCEPT_DECISIONS = {"stop_receiving"}
# Decisions that mean the read was ejected
REJECT_DECISIONS = {"unblock"}
# Everything else: the aligner couldn't decide yet
PROCEED_DECISIONS = {"proceed", "no_map", "no_seq",
                     "below_min_chunks", "above_max_chunks",
                     "duplex_override", "first_read_override"}

# ---------------------------------------------------------------------------
# Read ID parsing
# ---------------------------------------------------------------------------

def parse_read_id(read_id: str) -> dict:
    """
    Parse a seq2squiggle-style read ID of the form:
        <prefix>!<contig>!<start>!<end>!<strand>
    e.g. S1_64!Neisseria_subflava_contig_1!756462!756936!-

    Returns a dict with keys: prefix, contig, species, start, end, strand.
    'species' is the contig name stripped of the trailing _contig_N suffix,
    with underscores replaced by spaces.
    Falls back gracefully if the format doesn't match.
    """
    # Two formats are supported:
    #   seq2squiggle: <contig>!<start>!<end>!<strand>          (4 fields)
    #   squigulator:  <prefix>!<contig>!<start>!<end>!<strand> (5 fields)
    # The last four fields are always contig, start, end, strand regardless
    # of the tool, so parts[-4:] works for both without format detection.
    parts = read_id.split("!")
    if len(parts) >= 4:
        contig, start, end, strand = parts[-4], parts[-3], parts[-2], parts[-1]
        prefix = "!".join(parts[:-4]) if len(parts) > 4 else ""
        import re
        species = re.sub(r"_contig_\d+$", "", contig).replace("_", " ")
        return {
            "prefix":  prefix,
            "contig":  contig,
            "species": species,
            "start":   int(start) if start.isdigit() else None,
            "end":     int(end)   if end.isdigit()   else None,
            "strand":  strand,
        }
    return {"prefix": read_id, "contig": None, "species": None,
            "start": None, "end": None, "strand": None}


# ---------------------------------------------------------------------------
# TSV loading
# ---------------------------------------------------------------------------

def load_decisions_tsv(tsv_path: str) -> pd.DataFrame:
    """
    Load the Readfish per-read decision TSV.

    Columns used downstream:
        read_id, seq_len, mode, decision, action_override, timestamp
    """
    df = pd.read_csv(tsv_path, sep="\t")
    required = {"read_id", "seq_len", "mode", "decision",
                "action_override", "timestamp"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(
            f"Decision TSV {tsv_path} is missing columns: {missing}\n"
            f"Found: {list(df.columns)}"
        )
    # Normalise types
    df["seq_len"]        = pd.to_numeric(df["seq_len"], errors="coerce").fillna(0).astype(int)
    df["timestamp"]      = pd.to_numeric(df["timestamp"], errors="coerce")
    # action_override may be 'True'/'False' strings or booleans
    df["action_override"] = df["action_override"].apply(
        lambda x: str(x).strip().lower() == "true"
    )
    df["decision"] = df["decision"].str.strip().str.lower()
    logger.info("Loaded %d decision rows from %s", len(df), tsv_path)
    return df


def load_unblocked_ids(txt_path: Optional[str]) -> set:
    """Load unblocked_read_ids.txt into a set (optional cross-check)."""
    if txt_path is None or not Path(txt_path).exists():
        return set()
    ids = set()
    with open(txt_path) as fh:
        for line in fh:
            rid = line.strip()
            if rid:
                ids.add(rid)
    logger.info("Loaded %d unblocked read IDs from %s", len(ids), txt_path)
    return ids


# ---------------------------------------------------------------------------
# Ground truth
# ---------------------------------------------------------------------------

def load_manifest(manifest_path: str) -> Dict[str, int]:
    """
    Load a contig→community manifest TSV.

    Expected columns (tab-separated, with header):
        contig       contig name as it appears in the read ID
        true_label   0 = deplete (Community 0), 1 = target (Community 1)

    A species-level manifest is also accepted:
        species      species name (underscores or spaces)
        true_label   0 or 1

    Generate with make_manifest.py:
        python make_manifest.py \
            --community1-dir /data/gut/community1/ \
            --community0-dir /data/gut/community0/ \
            --output data/gut_manifest.tsv
    """
    df = pd.read_csv(manifest_path, sep="\t")
    # Accept either "contig" or "species" as the key column
    if "contig" in df.columns:
        key_col = "contig"
    elif "species" in df.columns:
        key_col = "species"
    else:
        raise ValueError(
            f"Manifest {manifest_path} must have a 'contig' or 'species' column. "
            f"Found: {list(df.columns)}"
        )
    if "true_label" not in df.columns:
        raise ValueError(f"Manifest {manifest_path} must have a 'true_label' column.")

    mapping = dict(zip(df[key_col].astype(str), df["true_label"].astype(int)))
    logger.info("Manifest: %d entries loaded from %s", len(mapping), manifest_path)
    return mapping


def resolve_true_label(
    contig: Optional[str],
    manifest: Optional[Dict[str, int]],
    community: Optional[int],
) -> int:
    """
    Resolve the ground-truth label (1=target, 0=deplete, -1=unknown) for one
    read using whichever ground-truth source is available.

    Priority:
      1. manifest lookup by contig name
      2. manifest lookup by species name (contig stripped of _contig_N suffix)
      3. --community flag (whole run is one community)
      4. -1 (unknown — read will be excluded from accuracy metrics)
    """
    import re
    if manifest is not None and contig is not None:
        # Try exact contig name first
        if contig in manifest:
            return manifest[contig]
        # Fall back to species name
        species = re.sub(r"_contig_\d+$", "", contig).replace("_", " ")
        if species in manifest:
            return manifest[species]
        # Try with underscores instead of spaces
        species_u = species.replace(" ", "_")
        if species_u in manifest:
            return manifest[species_u]
        return -1   # contig not in manifest
    if community is not None:
        return int(community)
    return -1


# ---------------------------------------------------------------------------
# Core per-run analysis
# ---------------------------------------------------------------------------

@dataclass
class CommunityRunResult:
    """Metrics for a single (tool, dataset) run. community=-1 when mixed."""
    tool:       str
    dataset:    str
    community:  int          # 0, 1, or -1 (mixed-community run)

    n_total:    int = 0      # all reads in TSV
    n_override: int = 0      # excluded: action_override == True
    n_decided:  int = 0      # accepted + rejected (override excluded)
    n_accepted: int = 0      # stop_receiving (override excluded)
    n_rejected: int = 0      # unblock        (override excluded)
    n_undecided:int = 0      # proceed / no_map etc (override excluded)

    # Timing (seconds between consecutive decisions — proxy for per-read latency)
    latency_mean_s:   float = float("nan")
    latency_median_s: float = float("nan")
    latency_p95_s:    float = float("nan")
    latency_p99_s:    float = float("nan")

    # Seq-len stats at decision (basecalled bases used to decide)
    seqlen_mean:   float = float("nan")
    seqlen_median: float = float("nan")
    seqlen_p95:    float = float("nan")

    # Per-species breakdown stored separately (not in the dataclass)


def analyse_community_run(
    tsv_path:   str,
    tool:       str,
    dataset:    str,
    outdir:     str,
    community:  Optional[int] = None,   # set for separate-community runs
    manifest:   Optional[Dict[str, int]] = None,  # set for mixed-community runs
    unblocked_txt: Optional[str] = None,
) -> Tuple[CommunityRunResult, pd.DataFrame]:
    """
    Analyse a single run.

    Ground truth is resolved per-read via resolve_true_label():
      - If manifest is provided: look up each read's contig in the manifest.
      - If community is provided: every read has that label.
      - Otherwise: all reads are labelled -1 (unknown, excluded from metrics).

    Returns (CommunityRunResult, per_read_df).
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_decisions_tsv(tsv_path)
    unblocked_ids = load_unblocked_ids(unblocked_txt)

    # Cross-check unblocked_read_ids.txt against TSV if both available
    if unblocked_ids:
        tsv_unblocked = set(df.loc[df["decision"] == "unblock", "read_id"])
        only_in_txt = unblocked_ids - tsv_unblocked
        only_in_tsv = tsv_unblocked - unblocked_ids
        if only_in_txt or only_in_tsv:
            logger.warning(
                "Unblocked IDs mismatch: %d only in txt, %d only in TSV",
                len(only_in_txt), len(only_in_tsv)
            )

    # Parse species from read_id
    parsed = df["read_id"].apply(parse_read_id).apply(pd.Series)
    df = pd.concat([df, parsed], axis=1)

    # Resolve ground truth per read using manifest, community flag, or read ID
    def _get_label(row):
        parsed = parse_read_id(row["read_id"])
        return resolve_true_label(parsed.get("contig"), manifest, community)

    df["true_label"]  = df.apply(_get_label, axis=1)
    df["is_override"] = df["action_override"]

    # is_correct: True if decision matches ground truth
    #   true_label=1 + stop_receiving → correct (TP)
    #   true_label=1 + unblock        → wrong   (FN)
    #   true_label=0 + unblock        → correct (TN)
    #   true_label=0 + stop_receiving → wrong   (FP)
    #   true_label=-1                 → unknown (excluded)
    def _is_correct(row):
        lbl = row["true_label"]
        dec = row["decision"]
        if lbl == 1:
            return dec in ACCEPT_DECISIONS
        if lbl == 0:
            return dec in REJECT_DECISIONS
        return None   # unknown

    df["is_correct"] = df.apply(_is_correct, axis=1)

    # Save enriched per-read table
    comm_label = community if community is not None else "mixed"

    out_perread = outdir / f"{tool}_{dataset}_community{comm_label}_per_read.tsv"
    df.to_csv(out_perread, sep="\t", index=False)
    logger.info("Saved per-read table → %s", out_perread)

    # ── Metrics (override-excluded) ────────────────────────────────────────
    # Exclude reads with unknown true_label from accuracy metrics
    known_mask    = df["true_label"] != -1
    decided_mask  = known_mask & (~df["is_override"]) & (df["decision"].isin(ACCEPT_DECISIONS | REJECT_DECISIONS))
    override_mask = df["is_override"]
    undecided_mask = known_mask & (~df["is_override"]) & (~df["decision"].isin(ACCEPT_DECISIONS | REJECT_DECISIONS))

    decided_df   = df[decided_mask]
    n_accepted   = int((decided_df["decision"].isin(ACCEPT_DECISIONS)).sum())
    n_rejected   = int((decided_df["decision"].isin(REJECT_DECISIONS)).sum())
    n_unknown    = int((~known_mask).sum())
    if n_unknown > 0:
        logger.warning("%d reads had no ground truth label (not in manifest) "
                       "— excluded from accuracy metrics", n_unknown)

    res = CommunityRunResult(
        tool=tool, dataset=dataset,
        community=community if community is not None else -1,
        n_total=len(df),
        n_override=int(override_mask.sum()),
        n_decided=len(decided_df),
        n_accepted=n_accepted,
        n_rejected=n_rejected,
        n_undecided=int(undecided_mask.sum()),
    )

    # Latency: inter-decision intervals from timestamp column, decided reads only
    ts = decided_df["timestamp"].dropna().sort_values().values
    if len(ts) > 1:
        intervals = np.diff(ts)
        intervals = intervals[(intervals > 0) & (intervals < 60)]  # sanity filter
        if len(intervals) > 0:
            res.latency_mean_s   = float(np.mean(intervals))
            res.latency_median_s = float(np.median(intervals))
            res.latency_p95_s    = float(np.percentile(intervals, 95))
            res.latency_p99_s    = float(np.percentile(intervals, 99))

    # Seq-len at decision
    sl = decided_df["seq_len"].dropna().values
    if len(sl) > 0:
        res.seqlen_mean   = float(np.mean(sl))
        res.seqlen_median = float(np.median(sl))
        res.seqlen_p95    = float(np.percentile(sl, 95))

    out_metrics = outdir / f"{tool}_{dataset}_community{comm_label}_metrics.json"
    with open(out_metrics, "w") as fh:
        json.dump(asdict(res), fh, indent=2)
    logger.info("Saved metrics → %s", out_metrics)

    return res, df


# ---------------------------------------------------------------------------
# Paired analysis: combine community0 and community1 runs → confusion matrix
# ---------------------------------------------------------------------------

@dataclass
class PairedRunMetrics:
    """Combined metrics for a (tool, dataset) pair across both communities."""
    tool:    str
    dataset: str

    # Confusion matrix counts (override-excluded)
    tp: int = 0   # community1 reads correctly accepted
    fn: int = 0   # community1 reads wrongly rejected
    tn: int = 0   # community0 reads correctly rejected
    fp: int = 0   # community0 reads wrongly accepted

    # Derived
    precision:  float = 0.0
    recall:     float = 0.0
    f1:         float = 0.0
    specificity:float = 0.0   # TN / (TN + FP)  — correct rejection rate

    # Per-community latency (seconds)
    c0_latency_mean_s:   float = float("nan")
    c0_latency_median_s: float = float("nan")
    c0_latency_p95_s:    float = float("nan")
    c1_latency_mean_s:   float = float("nan")
    c1_latency_median_s: float = float("nan")
    c1_latency_p95_s:    float = float("nan")

    # Seq-len at decision
    c0_seqlen_median: float = float("nan")
    c1_seqlen_median: float = float("nan")

    # Counts
    c0_n_total:   int = 0
    c0_n_decided: int = 0
    c1_n_total:   int = 0
    c1_n_decided: int = 0

    # Memory (populated separately if mem_log available)
    peak_rss_mb: float = float("nan")


def analyse_mixed_run(
    per_read_df: pd.DataFrame,
    tool: str,
    dataset: str,
) -> PairedRunMetrics:
    """
    Build a PairedRunMetrics from a single mixed-community run where
    true_label is set per-read from the manifest.
    """
    known = per_read_df[per_read_df["true_label"] != -1].copy()
    decided = known[
        (~known["is_override"]) &
        (known["decision"].isin(ACCEPT_DECISIONS | REJECT_DECISIONS))
    ]

    tp = int(((decided["decision"].isin(ACCEPT_DECISIONS)) & (decided["true_label"] == 1)).sum())
    fn = int(((decided["decision"].isin(REJECT_DECISIONS)) & (decided["true_label"] == 1)).sum())
    tn = int(((decided["decision"].isin(REJECT_DECISIONS)) & (decided["true_label"] == 0)).sum())
    fp = int(((decided["decision"].isin(ACCEPT_DECISIONS)) & (decided["true_label"] == 0)).sum())

    precision   = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall      = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1          = (2 * precision * recall / (precision + recall)
                   if (precision + recall) > 0 else 0.0)
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0

    def _latency(df_community: pd.DataFrame):
        """Compute inter-decision intervals for one community's decided reads."""
        ts = (df_community[~df_community["is_override"]]
              .loc[df_community["decision"].isin(ACCEPT_DECISIONS | REJECT_DECISIONS),
                   "timestamp"]
              .dropna()
              .sort_values()
              .values)
        if len(ts) < 2:
            return float("nan"), float("nan"), float("nan")
        iv = np.diff(ts)
        iv = iv[(iv > 0) & (iv < 60)]
        if len(iv) == 0:
            return float("nan"), float("nan"), float("nan")
        return float(np.mean(iv)), float(np.median(iv)), float(np.percentile(iv, 95))

    def _seqlen_median(df_community: pd.DataFrame):
        sl = (df_community[~df_community["is_override"]]
              .loc[df_community["decision"].isin(ACCEPT_DECISIONS | REJECT_DECISIONS),
                   "seq_len"]
              .dropna()
              .values)
        return float(np.median(sl)) if len(sl) > 0 else float("nan")

    c0 = known[known["true_label"] == 0]
    c1 = known[known["true_label"] == 1]

    c0_lat_mean, c0_lat_median, c0_lat_p95 = _latency(c0)
    c1_lat_mean, c1_lat_median, c1_lat_p95 = _latency(c1)

    return PairedRunMetrics(
        tool=tool, dataset=dataset,
        tp=tp, fn=fn, tn=tn, fp=fp,
        precision=round(precision, 4),
        recall=round(recall, 4),
        f1=round(f1, 4),
        specificity=round(specificity, 4),
        c0_latency_mean_s=c0_lat_mean,
        c0_latency_median_s=c0_lat_median,
        c0_latency_p95_s=c0_lat_p95,
        c1_latency_mean_s=c1_lat_mean,
        c1_latency_median_s=c1_lat_median,
        c1_latency_p95_s=c1_lat_p95,
        c0_seqlen_median=_seqlen_median(c0),
        c1_seqlen_median=_seqlen_median(c1),
        c1_n_total=int((known["true_label"] == 1).sum()),
        c1_n_decided=int((decided["true_label"] == 1).sum()),
        c0_n_total=int((known["true_label"] == 0).sum()),
        c0_n_decided=int((decided["true_label"] == 0).sum()),
    )


def combine_community_results(
    c0_result: CommunityRunResult,
    c1_result: CommunityRunResult,
) -> PairedRunMetrics:
    """
    Combine community 0 and community 1 run results into a paired metrics object.
    """
    assert c0_result.tool    == c1_result.tool
    assert c0_result.dataset == c1_result.dataset

    # From community 0 run (deplete):
    #   correct rejection (TN) = reads that were unblocked
    #   wrong acceptance  (FP) = reads that got stop_receiving
    tn = c0_result.n_rejected   # unblock on community0
    fp = c0_result.n_accepted   # stop_receiving on community0

    # From community 1 run (target):
    #   correct acceptance (TP) = reads that got stop_receiving
    #   missed             (FN) = reads that were unblocked
    tp = c1_result.n_accepted   # stop_receiving on community1
    fn = c1_result.n_rejected   # unblock on community1

    precision   = tp / (tp + fp) if (tp + fp) > 0 else 0.0
    recall      = tp / (tp + fn) if (tp + fn) > 0 else 0.0
    f1          = (2 * precision * recall / (precision + recall)
                   if (precision + recall) > 0 else 0.0)
    specificity = tn / (tn + fp) if (tn + fp) > 0 else 0.0

    return PairedRunMetrics(
        tool=c0_result.tool, dataset=c0_result.dataset,
        tp=tp, fn=fn, tn=tn, fp=fp,
        precision=round(precision,   4),
        recall=round(recall,         4),
        f1=round(f1,                 4),
        specificity=round(specificity, 4),
        c0_latency_mean_s=c0_result.latency_mean_s,
        c0_latency_median_s=c0_result.latency_median_s,
        c0_latency_p95_s=c0_result.latency_p95_s,
        c1_latency_mean_s=c1_result.latency_mean_s,
        c1_latency_median_s=c1_result.latency_median_s,
        c1_latency_p95_s=c1_result.latency_p95_s,
        c0_seqlen_median=c0_result.seqlen_median,
        c1_seqlen_median=c1_result.seqlen_median,
        c0_n_total=c0_result.n_total,
        c0_n_decided=c0_result.n_decided,
        c1_n_total=c1_result.n_total,
        c1_n_decided=c1_result.n_decided,
    )


# ---------------------------------------------------------------------------
# Per-species breakdown
# ---------------------------------------------------------------------------

def species_breakdown(per_read_df: pd.DataFrame,
                      community: int,
                      tool: str,
                      dataset: str,
                      outdir: Path) -> pd.DataFrame:
    """
    Compute per-species accuracy for a single community run.
    For community 1: reports TP rate (recall) per species.
    For community 0: reports TN rate (specificity) per species.
    """
    df = per_read_df[~per_read_df["is_override"]].copy()
    df = df[df["decision"].isin(ACCEPT_DECISIONS | REJECT_DECISIONS)]
    if "species" not in df.columns or df["species"].isna().all():
        logger.warning("No species information available for breakdown.")
        return pd.DataFrame()

    breakdown = (df.groupby("species")
                   .agg(
                       n_reads=("read_id", "count"),
                       n_correct=("is_correct", "sum"),
                       mean_seq_len=("seq_len", "mean"),
                   )
                   .assign(correct_rate=lambda x: x["n_correct"] / x["n_reads"])
                   .sort_values("correct_rate")
                   .reset_index())

    label = "recall" if community == 1 else "specificity"
    breakdown = breakdown.rename(columns={"correct_rate": label})

    comm_label = community if community is not None else "mixed"
    out = outdir / f"{tool}_{dataset}_community{comm_label}_species_breakdown.tsv"
    breakdown.to_csv(out, sep="\t", index=False)
    logger.info("Saved species breakdown → %s", out)
    return breakdown


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------

def _mpl():
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns
        sns.set_theme(style="whitegrid", font_scale=1.1)
        return plt, sns
    except ImportError:
        logger.error("matplotlib/seaborn not installed — skipping figures.")
        return None, None


def plot_precision_recall(df: pd.DataFrame, outdir: Path) -> None:
    plt, sns = _mpl()
    if plt is None or df.empty:
        return
    fig, ax = plt.subplots(figsize=(7, 6))
    markers = {"zymo": "o", "gut": "s"}
    colors  = {"zymo": "#1f77b4", "gut": "#d62728"}
    for dataset, grp in df.groupby("dataset"):
        ax.scatter(grp["recall"], grp["precision"],
                   marker=markers.get(dataset, "^"),
                   color=colors.get(dataset, "grey"),
                   s=90, zorder=3, label=dataset)
        for _, row in grp.iterrows():
            ax.annotate(row["tool"],
                        (row["recall"], row["precision"]),
                        textcoords="offset points", xytext=(5, 3), fontsize=8)
    ax.set_xlabel("Recall (sensitivity)", fontsize=12)
    ax.set_ylabel("Precision", fontsize=12)
    ax.set_xlim(-0.02, 1.05)
    ax.set_ylim(-0.02, 1.05)
    ax.plot([0, 1], [0, 1], ls="--", color="grey", lw=0.8)
    ax.legend(title="Dataset")
    ax.set_title("Precision vs Recall per tool")
    plt.tight_layout()
    out = outdir / "precision_recall.pdf"
    plt.savefig(out); plt.close()
    logger.info("→ %s", out)


def plot_f1_and_specificity(df: pd.DataFrame, outdir: Path) -> None:
    plt, sns = _mpl()
    if plt is None or df.empty:
        return
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=False)
    tool_order = sorted(df["tool"].unique())
    for ax, metric, title in zip(
        axes,
        ["f1", "specificity"],
        ["F1 score", "Specificity (correct rejection rate)"]
    ):
        sns.barplot(data=df, x="tool", y=metric, hue="dataset",
                    order=tool_order, ax=ax)
        ax.set_ylim(0, 1.05)
        ax.set_ylabel(title, fontsize=11)
        ax.set_xlabel("Tool", fontsize=11)
        ax.legend(title="Dataset")
    plt.suptitle("Classification performance per tool and dataset")
    plt.tight_layout()
    out = outdir / "f1_and_specificity.pdf"
    plt.savefig(out); plt.close()
    logger.info("→ %s", out)


def plot_latency_violin(all_per_read: pd.DataFrame, outdir: Path) -> None:
    """Violin of inter-decision intervals (latency proxy), decided reads only."""
    plt, sns = _mpl()
    if plt is None or all_per_read.empty:
        return
    if "timestamp" not in all_per_read.columns:
        return

    records = []
    for (tool, dataset, community), grp in all_per_read.groupby(
            ["tool", "dataset", "community"]):
        decided = grp[
            (~grp["is_override"]) &
            (grp["decision"].isin(ACCEPT_DECISIONS | REJECT_DECISIONS))
        ]
        ts = decided["timestamp"].dropna().sort_values().values
        if len(ts) < 2:
            continue
        intervals = np.diff(ts)
        intervals = intervals[(intervals > 0) & (intervals < 60)]
        for iv in intervals:
            records.append({
                "tool": tool,
                "dataset": dataset,
                "community": f"c{community}",
                "latency_ms": iv * 1000,
            })

    if not records:
        logger.warning("No latency intervals — skipping violin.")
        return

    lat_df = pd.DataFrame(records)
    tool_order = sorted(lat_df["tool"].unique())
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
    for ax, dataset in zip(axes, ["zymo", "gut"]):
        sub = lat_df[lat_df["dataset"] == dataset]
        if sub.empty:
            ax.set_title(f"{dataset} (no data)")
            continue
        sns.violinplot(data=sub, x="tool", y="latency_ms", hue="community",
                       split=True, inner="quartile", ax=ax, order=tool_order,
                       cut=0)
        ax.set_yscale("log")
        ax.set_ylabel("Latency (ms, log scale)" if ax == axes[0] else "")
        ax.set_xlabel("Tool")
        ax.set_title(f"Dataset: {dataset}")
        ax.legend(title="Community")
    plt.suptitle("Decision latency distribution per tool")
    plt.tight_layout()
    out = outdir / "latency_violin.pdf"
    plt.savefig(out); plt.close()
    logger.info("→ %s", out)


def plot_seqlen_at_decision(all_per_read: pd.DataFrame, outdir: Path) -> None:
    """CDF of basecalled sequence length at the moment of decision."""
    plt, sns = _mpl()
    if plt is None or all_per_read.empty:
        return
    decided = all_per_read[
        (~all_per_read["is_override"]) &
        (all_per_read["decision"].isin(ACCEPT_DECISIONS | REJECT_DECISIONS))
    ]
    if decided.empty:
        return
    tool_order = sorted(decided["tool"].unique())
    fig, axes = plt.subplots(1, 2, figsize=(13, 5), sharey=True)
    for ax, dataset in zip(axes, ["zymo", "gut"]):
        sub = decided[decided["dataset"] == dataset]
        if sub.empty:
            continue
        for tool in tool_order:
            t = sub[sub["tool"] == tool]["seq_len"].dropna().sort_values().values
            if len(t) == 0:
                continue
            cdf = np.arange(1, len(t) + 1) / len(t)
            ax.plot(t, cdf, label=tool)
        ax.set_xlabel("Basecalled bases at decision", fontsize=11)
        ax.set_ylabel("CDF" if ax == axes[0] else "")
        ax.set_title(f"Dataset: {dataset}")
        ax.legend(fontsize=8)
        ax.axvline(450, color="grey", lw=0.8, ls="--", label="450 bp")
    plt.suptitle("Sequence length used for decision (CDF)")
    plt.tight_layout()
    out = outdir / "seqlen_at_decision_cdf.pdf"
    plt.savefig(out); plt.close()
    logger.info("→ %s", out)


def plot_species_heatmap(species_dfs: list, outdir: Path) -> None:
    """
    Heatmap of per-species correct-decision rate across tools,
    one heatmap per (dataset, community).
    species_dfs is a list of dicts:
      {tool, dataset, community, df}  where df has columns species + recall|specificity.
    """
    plt, sns = _mpl()
    if plt is None:
        return
    from itertools import product as iproduct

    for dataset in set(d["dataset"] for d in species_dfs):
        for community in [0, 1]:
            subset = [d for d in species_dfs
                      if d["dataset"] == dataset and d["community"] == community]
            if not subset:
                continue
            metric_col = "recall" if community == 1 else "specificity"
            # Build matrix: rows = species, cols = tools
            combined = {}
            for entry in subset:
                df = entry["df"]
                if df.empty or metric_col not in df.columns:
                    continue
                combined[entry["tool"]] = df.set_index("species")[metric_col]
            if not combined:
                continue
            matrix = pd.DataFrame(combined).fillna(float("nan"))
            if matrix.empty:
                continue
            fig, ax = plt.subplots(figsize=(max(6, len(matrix.columns) * 1.4),
                                             max(4, len(matrix) * 0.35)))
            sns.heatmap(matrix, annot=True, fmt=".2f", cmap="RdYlGn",
                        vmin=0, vmax=1, ax=ax,
                        linewidths=0.3, linecolor="grey")
            label = "Recall" if community == 1 else "Specificity"
            ax.set_title(f"{dataset.capitalize()} Community {community} — "
                         f"per-species {label.lower()} per tool")
            ax.set_ylabel("Species")
            ax.set_xlabel("Tool")
            plt.tight_layout()
            out = outdir / f"species_heatmap_{dataset}_c{community}.pdf"
            plt.savefig(out, bbox_inches="tight"); plt.close()
            logger.info("→ %s", out)


def plot_summary_table_image(df: pd.DataFrame, outdir: Path) -> None:
    """Render the summary metrics table as a PDF/PNG for inclusion in paper."""
    plt, _ = _mpl()
    if plt is None or df.empty:
        return
    cols = ["tool", "dataset", "precision", "recall", "f1", "specificity",
            "c1_latency_median_s", "c0_latency_median_s",
            "c1_seqlen_median", "c0_seqlen_median"]
    cols = [c for c in cols if c in df.columns]
    table = df[cols].sort_values(["dataset", "f1"], ascending=[True, False])
    # Rename for display
    rename = {
        "c1_latency_median_s": "lat_c1_med_s",
        "c0_latency_median_s": "lat_c0_med_s",
        "c1_seqlen_median":    "seqlen_c1_med",
        "c0_seqlen_median":    "seqlen_c0_med",
    }
    table = table.rename(columns=rename)
    fig, ax = plt.subplots(figsize=(14, max(2, len(table) * 0.45 + 1)))
    ax.axis("off")
    t = ax.table(cellText=table.round(4).values,
                 colLabels=table.columns,
                 cellLoc="center", loc="center")
    t.auto_set_font_size(False)
    t.set_fontsize(8)
    t.auto_set_column_width(col=list(range(len(table.columns))))
    ax.set_title("Summary metrics across all tools and datasets", fontsize=10, pad=10)
    plt.tight_layout()
    out = outdir / "summary_table.pdf"
    plt.savefig(out, bbox_inches="tight"); plt.close()
    logger.info("→ %s", out)



# ---------------------------------------------------------------------------
# Memory log loading
# ---------------------------------------------------------------------------

def load_mem_log(mem_log_path: str) -> float:
    """
    Load a readfish_mem_<timestamp>.tsv written by simulate_run.sh and return
    peak RSS in MB.

    Expected format (tab-separated, with header):
        timestamp_s    rss_mb
        1234567890.1   512.34
    """
    try:
        df = pd.read_csv(mem_log_path, sep="\t")
        if "rss_mb" not in df.columns:
            logger.warning("mem log %s has no \'rss_mb\' column -- skipping", mem_log_path)
            return float("nan")
        peak = float(df["rss_mb"].max())
        logger.info("Peak RSS %.1f MB from %s", peak, mem_log_path)
        return peak
    except Exception as e:
        logger.warning("Could not load mem log %s: %s", mem_log_path, e)
        return float("nan")


def _find_mem_log(tool_dir: Path) -> Optional[str]:
    """
    Auto-discover a readfish_mem_*.tsv in tool_dir or its logs/ sibling.
    Returns the path of the most recently modified match, or None.
    """
    candidates = sorted(
        list(tool_dir.glob("readfish_mem_*.tsv")) +
        list(tool_dir.glob("logs/readfish_mem_*.tsv")),
        key=lambda p: p.stat().st_mtime,
        reverse=True,
    )
    return str(candidates[0]) if candidates else None


# ---------------------------------------------------------------------------
# Batch discovery
# ---------------------------------------------------------------------------

def discover_runs(results_dir: str, dataset: str,
                  mixed: bool = False) -> List[dict]:
    """
    Discover runs under results_dir.

    Separate-community layout (mixed=False):
        results_dir/<tool>/community0/readfish_decisions.tsv
        results_dir/<tool>/community1/readfish_decisions.tsv

    Mixed-community layout (mixed=True):
        results_dir/<tool>/readfish_decisions.tsv
    """
    runs = []
    root = Path(results_dir)
    for tool_dir in sorted(root.iterdir()):
        if not tool_dir.is_dir():
            continue
        tool = tool_dir.name

        if mixed:
            tsv = tool_dir / "readfish_decisions.tsv"
            txt = tool_dir / "unblocked_read_ids.txt"
            if not tsv.exists():
                logger.warning("Missing %s — skipping %s", tsv, tool)
                continue
            runs.append({
                "tool":      tool,
                "dataset":   dataset,
                "community": None,          # resolved per-read from manifest
                "tsv":       str(tsv),
                "unblocked": str(txt) if txt.exists() else None,
                "outdir":    str(tool_dir),
            })
        else:
            for community in [0, 1]:
                comm_dir = tool_dir / f"community{community}"
                tsv = comm_dir / "readfish_decisions.tsv"
                txt = comm_dir / "unblocked_read_ids.txt"
                if not tsv.exists():
                    logger.warning("Missing %s — skipping %s/community%d",
                                tsv, tool, community)
                    continue
                runs.append({
                        "tool":      tool,
                        "dataset":   dataset,
                        "community": community,
                        "tsv":       str(tsv),
                        "unblocked": str(txt) if txt.exists() else None,
                        "outdir":    str(comm_dir),
                })

    logger.info("Discovered %d run(s) under %s (mixed=%s)",
                len(runs), results_dir, mixed)
    return runs


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def cmd_single(args: argparse.Namespace) -> None:
    manifest = load_manifest(args.manifest) if args.manifest else None
    community = getattr(args, "community", None)

    if manifest is None and community is None:
        logger.error("Provide either --manifest or --community.")
        sys.exit(1)

    result, per_read_df = analyse_community_run(
        tsv_path=args.tsv,
        tool=args.tool,
        dataset=args.dataset,
        outdir=args.outdir,
        community=community,
        manifest=manifest,
        unblocked_txt=args.unblocked,
    )
    per_read_df["tool"]    = args.tool
    per_read_df["dataset"] = args.dataset

    # Mixed run: build paired metrics directly from per_read_df
    if manifest is not None:
        pm = analyse_mixed_run(per_read_df, args.tool, args.dataset)
        if args.mem_log:
            pm.peak_rss_mb = load_mem_log(args.mem_log)
        out = Path(args.outdir) / f"{args.tool}_{args.dataset}_paired_metrics.json"
        with open(out, "w") as fh:
            json.dump(asdict(pm), fh, indent=2)
        print(json.dumps(asdict(pm), indent=2))
    else:
        sp_df = species_breakdown(per_read_df, community,
                               args.tool, args.dataset, Path(args.outdir))
    print(json.dumps(asdict(result), indent=2))


def cmd_batch(args: argparse.Namespace) -> None:
    manifest = load_manifest(args.manifest) if args.manifest else None
    runs = discover_runs(args.results_dir, args.dataset, mixed=(manifest is not None))
    if not runs:
        logger.error("No runs found.")
        sys.exit(1)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    # Analyse each community run
    community_results: Dict[str, Dict[int, CommunityRunResult]] = {}
    all_per_read_frames = []
    species_dfs_list = []

    for run in runs:
        tool, community = run["tool"], run["community"]
        logger.info("── %s / community%s ──", tool, community if community is not None else "mixed")
        try:
            result, per_read_df = analyse_community_run(
                tsv_path=run["tsv"],
                tool=tool,
                dataset=args.dataset,
                outdir=run["outdir"],
                community=run.get("community"),
                manifest=manifest,
                unblocked_txt=run["unblocked"],
            )
        except Exception as exc:
            logger.error("Failed: %s / community%s: %s", tool,
                         community if community is not None else "mixed", exc)
            continue

        per_read_df["tool"]      = tool
        per_read_df["dataset"]   = args.dataset
        per_read_df["community"] = community
        all_per_read_frames.append(per_read_df)

        sp_df = species_breakdown(per_read_df, community, tool,
                                   args.dataset, Path(run["outdir"]))
        if not sp_df.empty:
            species_dfs_list.append({
                "tool": tool, "dataset": args.dataset,
                "community": community, "df": sp_df
            })

        if tool not in community_results:
            community_results[tool] = {}
        community_results[tool][community] = result

    # Build paired metrics — either from mixed per_read_df or from paired community runs
    paired_metrics: List[PairedRunMetrics] = []
    if manifest is not None:
        all_per_read_for_metrics = pd.concat(all_per_read_frames, ignore_index=True) \
                                   if all_per_read_frames else pd.DataFrame()
        for tool in all_per_read_for_metrics["tool"].unique():
            tool_df = all_per_read_for_metrics[all_per_read_for_metrics["tool"] == tool]
            pm = analyse_mixed_run(tool_df, tool, args.dataset)
            # Populate peak RSS — explicit log dir takes priority, else auto-discover
            mem_log = None
            if args.mem_log_dir:
                tool_dir = Path(args.mem_log_dir) / tool
                mem_log = _find_mem_log(tool_dir)
            else:
                tool_dir = Path(args.results_dir) / tool
                mem_log = _find_mem_log(tool_dir)
            if mem_log:
                pm.peak_rss_mb = load_mem_log(mem_log)
            paired_metrics.append(pm)
            out = outdir / f"{tool}_{args.dataset}_paired_metrics.json"
            with open(out, "w") as fh:
                json.dump(asdict(pm), fh, indent=2)
            logger.info("Paired metrics → %s", out)
    else:
        # Separate community runs: pair community0 and community1 results
        for tool, comm_map in community_results.items():
            if 0 not in comm_map or 1 not in comm_map:
                logger.warning("Tool %s missing one community run — skipping.", tool)
                continue
            pm = combine_community_results(comm_map[0], comm_map[1])
            mem_log = None
            if args.mem_log_dir:
                tool_dir = Path(args.mem_log_dir) / tool
                mem_log = _find_mem_log(tool_dir)
            else:
                tool_dir = Path(args.results_dir) / tool
                mem_log = _find_mem_log(tool_dir)
            if mem_log:
                pm.peak_rss_mb = load_mem_log(mem_log)
            paired_metrics.append(pm)
            out = outdir / f"{tool}_{args.dataset}_paired_metrics.json"
            with open(out, "w") as fh:
                json.dump(asdict(pm), fh, indent=2)
            logger.info("Paired metrics → %s", out)

    if not paired_metrics:
        logger.error("No paired metrics produced.")
        sys.exit(1)

    # Aggregate into DataFrames
    metrics_df   = pd.DataFrame([asdict(pm) for pm in paired_metrics])
    all_per_read = pd.concat(all_per_read_frames, ignore_index=True) \
                   if all_per_read_frames else pd.DataFrame()

    # Save master summary CSV
    summary_csv = outdir / f"{args.dataset}_metrics_summary.csv"
    metrics_df.to_csv(summary_csv, index=False)
    logger.info("Summary CSV → %s", summary_csv)
    print(metrics_df.to_string(index=False))

    # Generate all figures
    plot_precision_recall(metrics_df, outdir)
    plot_f1_and_specificity(metrics_df, outdir)
    if not all_per_read.empty:
        plot_latency_violin(all_per_read, outdir)
        plot_seqlen_at_decision(all_per_read, outdir)
    plot_species_heatmap(species_dfs_list, outdir)
    plot_summary_table_image(metrics_df, outdir)


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    p.add_argument("--log-level", default="INFO",
                   choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    sub = p.add_subparsers(dest="command", required=True)

    # ── single ──────────────────────────────────────────────────────────────
    s = sub.add_parser("single", help="Analyse one community run")
    s.add_argument("--tsv",       required=True, help="readfish_decisions.tsv")
    s.add_argument("--tool",      required=True, help="Tool name")
    s.add_argument("--dataset",   required=True, help="Dataset name (zymo / gut)")
    gt = s.add_mutually_exclusive_group(required=True)
    gt.add_argument("--community", type=int, choices=[0, 1],
                    help="Whole run is one community (0=deplete, 1=target)")
    gt.add_argument("--manifest",
                    help="TSV mapping contig→true_label for mixed-community runs")
    s.add_argument("--unblocked", default=None,
                   help="unblocked_read_ids.txt (optional cross-check)")
    s.add_argument("--mem-log",   default=None,
                   help="readfish_mem_<timestamp>.tsv from simulate_run.sh "
                        "(populates peak_rss_mb in paired metrics)")
    s.add_argument("--outdir",    default="./output")
    s.set_defaults(func=cmd_single)

    # ── batch ───────────────────────────────────────────────────────────────
    b = sub.add_parser("batch",
                       help="Analyse all tools under a results directory")
    b.add_argument("--results-dir", required=True,
                   help="Root dir: <results-dir>/<tool>/community{0,1}/")
    b.add_argument("--dataset",     required=True,
                   help="Dataset label (zymo / gut)")
    b.add_argument("--manifest",
                   help="TSV mapping contig→true_label (enables mixed-community mode)")
    b.add_argument("--mem-log-dir", default=None,
                   help="Directory containing per-tool memory logs, laid out as "
                        "<mem-log-dir>/<tool>/readfish_mem_*.tsv. "
                        "If omitted, auto-discovery looks under --results-dir/<tool>/. "
                        "Populates peak_rss_mb in paired metrics.")
    b.add_argument("--outdir",      default="./figures",
                   help="Where to write aggregate figures and CSVs")
    b.set_defaults(func=cmd_batch)

    return p


def main():
    parser = build_parser()
    args = parser.parse_args()
    logging.basicConfig(
        level=getattr(logging, args.log_level),
        format="%(asctime)s %(levelname)s  %(message)s",
        datefmt="%H:%M:%S",
    )
    args.func(args)


if __name__ == "__main__":
    main()