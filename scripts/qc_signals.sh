#!/usr/bin/env bash
# =============================================================================
# qc_signals.sh
# =============================================================================
# Quality control for simulated nanopore signals (BLOW5/SLOW5 format).
#
# Pipeline:
#   BLOW5 → basecall_blow5.py (ont-pybasecall-client-lib → Dorado server)
#   → align (Minimap2) → per-read stats → summary TSV
#
# Basecalling is done by talking directly to your already-running Dorado
# basecall server via ont-pybasecall-client-lib. This avoids the missing-
# metadata crash in slow5-dorado standalone and is guaranteed compatible
# with whichever server version you have (versions must match).
#
# Additionally verifies that each read maps back to its source contig using
# the encoded read ID format: <contig_name>!<start>!<end>!<strand>
#
# Usage:
#   bash qc_signals.sh \
#       --blow5          /path/to/community1.blow5 \
#       --ref            /path/to/community1.fasta \
#       --outdir         /path/to/qc_output/ \
#       --label          zymo_community1 \
#       [--threads       16] \
#       [--dorado-config dna_r10.4.1_e8.2_400bps_fast@v5.2.0||] \
#       [--dorado-address ipc:///var/lib/minknow/data/.dorado/dorado-basecall-server.sock]
#
# Requires basecall_blow5.py in the same directory as this script.
#
# Outputs (all written to --outdir):
#   <label>_basecalled.fastq       Raw basecalled reads
#   <label>_aligned.bam            Minimap2 alignment (sorted + indexed)
#   <label>_per_read.tsv           Per-read: read_id, q_score, read_len,
#                                  mapped, identity, source_contig,
#                                  aligned_contig, contig_correct
#   <label>_summary.tsv            Single-row summary for aggregation across runs
#   <label>_nanostat/              NanoStat report on the FASTQ
# =============================================================================
set -euo pipefail

# Directory containing this script — used to locate basecall_blow5.py
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
BLOW5=""
REF=""
OUTDIR=""
LABEL="signals"
THREADS=16
DORADO_CONFIG="dna_r10.4.1_e8.2_400bps_fast@v5.2.0||"
DORADO_ADDRESS="ipc:///var/lib/minknow/data/.dorado/dorado-basecall-server.sock"
SKIP_BASECALL=false

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
usage() {
    grep '^#' "$0" | grep -v '#!/' | sed 's/^# \?//'
    exit 1
}

while [[ $# -gt 0 ]]; do
    case $1 in
        --blow5)             BLOW5="$2";            shift 2 ;;
        --ref)               REF="$2";              shift 2 ;;
        --outdir)            OUTDIR="$2";           shift 2 ;;
        --label)             LABEL="$2";            shift 2 ;;
        --threads)           THREADS="$2";          shift 2 ;;
        --dorado-config)     DORADO_CONFIG="$2";    shift 2 ;;
        --dorado-address)    DORADO_ADDRESS="$2";   shift 2 ;;
        --skip-basecall)     SKIP_BASECALL=true;     shift ;;
        -h|--help)           usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------
[[ -z "$BLOW5"  ]] && { echo "ERROR: --blow5 is required";  exit 1; }
[[ -z "$REF"   ]] && { echo "ERROR: --ref is required";    exit 1; }
[[ -z "$OUTDIR" ]] && { echo "ERROR: --outdir is required"; exit 1; }
[[ -f "$BLOW5" ]]  || { echo "ERROR: BLOW5 file not found: $BLOW5"; exit 1; }
[[ -f "$REF"   ]]  || { echo "ERROR: Reference FASTA not found: $REF"; exit 1; }

for tool in minimap2 samtools NanoStat python3; do
    command -v "$tool" &>/dev/null || { echo "ERROR: $tool not found on PATH"; exit 1; }
done
if [[ "$SKIP_BASECALL" == false ]]; then
    python3 -c "import pybasecall_client_lib" 2>/dev/null || { echo "ERROR: ont-pybasecall-client-lib not installed. Run: pip install ont-pybasecall-client-lib==7.13.6"; exit 1; }
fi
[[ -f "$SCRIPT_DIR/basecall_blow5.py" ]] || { echo "ERROR: basecall_blow5.py not found in $SCRIPT_DIR"; exit 1; }

mkdir -p "$OUTDIR"

# ---------------------------------------------------------------------------
# File paths
# ---------------------------------------------------------------------------
FASTQ="$OUTDIR/${LABEL}_basecalled.fastq"
BAM="$OUTDIR/${LABEL}_aligned.bam"
PER_READ_TSV="$OUTDIR/${LABEL}_per_read.tsv"
SUMMARY_TSV="$OUTDIR/${LABEL}_summary.tsv"
NANOSTAT_DIR="$OUTDIR/${LABEL}_nanostat"
LOG="$OUTDIR/${LABEL}_qc.log"

echo "========================================"
echo " Signal QC: $LABEL"
echo " BLOW5:     $BLOW5"
echo " Reference: $REF"
echo " Output:    $OUTDIR"
echo " Threads:   $THREADS"
echo " Dorado:    $DORADO_ADDRESS ($DORADO_CONFIG)"
echo "========================================"
echo ""

# ---------------------------------------------------------------------------
# Step 1 — Basecall via Dorado basecall server (ont-pybasecall-client-lib)
#
# slow5-dorado standalone crashes on seq2squiggle BLOW5 files because they
# lack metadata fields (protocol_run_id, sample_id, protocol_start_time)
# that a real MinKNOW run would populate. The basecall server is tolerant of
# this because it is designed for live sequencing where incomplete run
# metadata is normal.
#
# Reads are loaded from BLOW5 via pyslow5 and submitted directly to the
# running Dorado basecall server. ont-pybasecall-client-lib version must
# exactly match the server version (both 7.13.6 here).
# ---------------------------------------------------------------------------
if [[ "$SKIP_BASECALL" == true ]]; then
    echo "[1/5] Skipping basecalling (--skip-basecall set)."
    [[ -f "$FASTQ" ]] || { echo "ERROR: --skip-basecall set but FASTQ not found: $FASTQ"; exit 1; }
elif [[ -f "$FASTQ" ]]; then
    echo "[1/5] Skipping basecalling — FASTQ already exists: $FASTQ"
else
    echo "[1/5] Basecalling via Dorado server ($DORADO_ADDRESS)..."
    python3 "$SCRIPT_DIR/basecall_blow5.py" \
        --blow5    "$BLOW5" \
        --fastq    "$FASTQ" \
        --config   "$DORADO_CONFIG" \
        --address  "$DORADO_ADDRESS" \
        --log-file "$LOG"
fi

N_READS_BASECALLED=$(grep -c "^@" "$FASTQ" || true)
echo "      Done: $N_READS_BASECALLED reads basecalled → $FASTQ"
printf "\a"  # beep

# ---------------------------------------------------------------------------
# Step 3 — Align to reference with Minimap2
# --secondary=no: suppress secondary alignments so each read has at most one
#   reported alignment, simplifying per-read identity calculation.
# -Y: use soft-clipping for supplementary alignments.
# ---------------------------------------------------------------------------
echo "[2/5] Aligning to reference..."
minimap2 \
    -ax map-ont \
    -t "$THREADS" \
    --secondary=no \
    -Y \
    "$REF" \
    "$FASTQ" \
    2>>"$LOG" \
    | samtools sort -@ "$THREADS" -o "$BAM"
printf "\a" # beep    

samtools index "$BAM"
printf "\a" # beep

N_TOTAL=$(samtools view -c "$BAM")
N_MAPPED=$(samtools view -c -F 4 "$BAM")
echo "      Done: $N_MAPPED / $N_TOTAL reads mapped → $BAM"
printf "\a" # beep

# ---------------------------------------------------------------------------
# Step 4 — NanoStat on the FASTQ
# Gives Q-score distribution, length distribution, N50.
# ---------------------------------------------------------------------------
echo "[3/5] Running NanoStat..."
mkdir -p "$NANOSTAT_DIR"
NanoStat \
    --fastq "$FASTQ" \
    --outdir "$NANOSTAT_DIR" \
    --name "$LABEL" \
    --threads "$THREADS" \
    >>"$LOG" 2>&1
echo "      Done: NanoStat report in $NANOSTAT_DIR"
printf "\a" # beep

# ---------------------------------------------------------------------------
# Step 5 — Per-read stats
# For each mapped read, extract:
#   read_id, q_score, read_len, identity, source_contig (from encoded ID),
#   aligned_contig (from BAM), contig_correct (bool)
#
# Identity = (aligned_length - NM) / aligned_length  where NM = edit distance.
# We use aligned_length (CIGAR-based) rather than read length to avoid
# penalising soft-clipped ends in the identity calculation.
# ---------------------------------------------------------------------------
echo "[4/5] Computing per-read statistics..."

# Extract Q-scores from FASTQ into a lookup table (read_id → mean_q)
# Mean Q from Phred: mean_q = -10 * log10(mean_error_prob)
python3 - "$FASTQ" "$BAM" "$PER_READ_TSV" << 'PYEOF'
import sys
import math
import pysam

fastq_path = sys.argv[1]
bam_path   = sys.argv[2]
out_path   = sys.argv[3]

# ── Build read_id → mean_q_score from FASTQ ────────────────────────────────
q_scores = {}
with pysam.FastxFile(fastq_path) as fq:
    for entry in fq:
        quals = entry.get_quality_array()
        if quals is not None and len(quals) > 0:
            # Convert Phred to error probs, take mean, back to Phred
            mean_error = sum(10 ** (-q / 10) for q in quals) / len(quals)
            mean_q = -10 * math.log10(mean_error) if mean_error > 0 else 0
        else:
            mean_q = float("nan")
        q_scores[entry.name] = round(mean_q, 3)

# ── Parse BAM for alignment stats ──────────────────────────────────────────
rows = []
with pysam.AlignmentFile(bam_path, "rb") as bam:
    for read in bam.fetch():
        read_id      = read.query_name
        mapped       = not read.is_unmapped
        read_len     = read.query_length or 0
        mean_q       = q_scores.get(read_id, float("nan"))

        # Parse encoded read ID. Two formats are supported:
        #   seq2squiggle: <contig>!<start>!<end>!<strand>          (4 fields)
        #   squigulator:  <prefix>!<contig>!<start>!<end>!<strand> (5 fields)
        # The contig name is always 4th from the end, so parts[-4] works for
        # both formats without needing to know which tool generated the file.
        parts = read_id.split("!")
        source_contig = parts[-4] if len(parts) >= 4 else "unknown"

        if mapped:
            aligned_contig  = read.reference_name
            nm              = read.get_tag("NM") if read.has_tag("NM") else None
            aligned_len     = read.query_alignment_length or 0
            if nm is not None and aligned_len > 0:
                identity = round((aligned_len - nm) / aligned_len * 100, 3)
            else:
                identity = float("nan")
            contig_correct = (aligned_contig == source_contig)
        else:
            aligned_contig = "*"
            identity       = float("nan")
            contig_correct = False

        rows.append({
            "read_id":        read_id,
            "mean_q":         mean_q,
            "read_len":       read_len,
            "mapped":         int(mapped),
            "identity":       identity,
            "source_contig":  source_contig,
            "aligned_contig": aligned_contig,
            "contig_correct": int(contig_correct),
        })

# ── Write TSV ──────────────────────────────────────────────────────────────
cols = ["read_id","mean_q","read_len","mapped","identity",
        "source_contig","aligned_contig","contig_correct"]

with open(out_path, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    for row in rows:
        fh.write("\t".join(str(row[c]) for c in cols) + "\n")

print(f"      Wrote {len(rows)} rows to {out_path}")
PYEOF

echo "      Done: per-read stats → $PER_READ_TSV"
printf "\a" # beep

# ---------------------------------------------------------------------------
# Step 6 — Summary TSV
# One row per run, suitable for aggregating across communities/tools/datasets.
# ---------------------------------------------------------------------------
echo "[5/5] Computing summary statistics..."

python3 - "$PER_READ_TSV" "$LABEL" "$BLOW5" "$SUMMARY_TSV" << 'PYEOF'
import sys
import math
import os
import numpy as np

per_read_path = sys.argv[1]
label         = sys.argv[2]
blow5_path    = sys.argv[3]
out_path      = sys.argv[4]

rows = []
with open(per_read_path) as fh:
    header = fh.readline().strip().split("\t")
    for line in fh:
        vals = line.strip().split("\t")
        rows.append(dict(zip(header, vals)))

n_total    = len(rows)
n_mapped   = sum(1 for r in rows if r["mapped"] == "1")
n_correct  = sum(1 for r in rows if r["contig_correct"] == "1")

mapping_rate   = round(n_mapped / n_total * 100, 2)   if n_total   > 0 else float("nan")
correct_rate   = round(n_correct / n_mapped * 100, 2) if n_mapped  > 0 else float("nan")

def safe_floats(rows, key):
    vals = []
    for r in rows:
        try:
            v = float(r[key])
            if not math.isnan(v):
                vals.append(v)
        except (ValueError, KeyError):
            pass
    return vals

identities = safe_floats([r for r in rows if r["mapped"] == "1"], "identity")
q_scores   = safe_floats(rows, "mean_q")
lengths    = safe_floats(rows, "read_len")

def pct(vals, p):
    return round(float(np.percentile(vals, p)), 3) if vals else float("nan")

summary = {
    "label":              label,
    "blow5_file":         os.path.basename(blow5_path),
    "n_reads_total":      n_total,
    "n_reads_mapped":     n_mapped,
    "n_reads_correct_contig": n_correct,
    "mapping_rate_pct":   mapping_rate,
    "correct_contig_rate_pct": correct_rate,
    "identity_mean":      round(float(np.mean(identities)),   3) if identities else float("nan"),
    "identity_median":    pct(identities, 50),
    "identity_p10":       pct(identities, 10),   # lower tail matters most
    "identity_p90":       pct(identities, 90),
    "q_mean":             round(float(np.mean(q_scores)),     3) if q_scores   else float("nan"),
    "q_median":           pct(q_scores,   50),
    "readlen_mean":       round(float(np.mean(lengths)),      1) if lengths     else float("nan"),
    "readlen_median":     pct(lengths,    50),
    "readlen_n50":        None,  # computed below
}

# N50
if lengths:
    sorted_lens = sorted(lengths, reverse=True)
    half_total  = sum(sorted_lens) / 2
    cumsum = 0
    n50 = sorted_lens[-1]
    for l in sorted_lens:
        cumsum += l
        if cumsum >= half_total:
            n50 = l
            break
    summary["readlen_n50"] = round(n50, 1)

cols = list(summary.keys())
with open(out_path, "w") as fh:
    fh.write("\t".join(cols) + "\n")
    fh.write("\t".join(str(summary[c]) for c in cols) + "\n")

print(f"      Wrote summary → {out_path}")
print()
print("  ┌─────────────────────────────────────────┐")
print(f"  │  Mapping rate:        {summary['mapping_rate_pct']:>6.1f}%           │")
print(f"  │  Correct contig rate: {summary['correct_contig_rate_pct']:>6.1f}%           │")
print(f"  │  Mean identity:       {summary['identity_mean']:>6.2f}%           │")
print(f"  │  Median Q-score:      {summary['q_median']:>6.2f}            │")
print(f"  │  Read length N50:     {summary['readlen_n50']:>6.0f} bp         │")
print("  └─────────────────────────────────────────┘")
PYEOF

echo ""
echo "✅ QC complete. Summary at: $SUMMARY_TSV"
echo "   Per-read stats at:       $PER_READ_TSV"
echo "   NanoStat report at:      $NANOSTAT_DIR/"

# No intermediate files to clean up.

echo ""
echo "Done."
printf "\a" # beep