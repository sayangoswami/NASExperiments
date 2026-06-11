#!/usr/bin/env bash
# =============================================================================
# run_seq2squiggle.sh
# =============================================================================
# Run seq2squiggle --read-input on each FASTA chunk produced by sample_reads.py
# and merge the resulting BLOW5 files into a single output.
#
# Each chunk is processed independently, so the job can be parallelised with
# --parallel-jobs. Failed chunks are retried once and logged; the script does
# not abort on a single chunk failure so the merge still proceeds.
#
# Usage:
#   bash run_seq2squiggle.sh \
#       --manifest  community1_reads/manifest.tsv \
#       --output    community1.blow5 \
#       --profile   dna-r10-prom \
#       [--parallel-jobs 8] \
#       [--seq2squiggle-args "--noise-std 1.0 --sample-rate 4000"] \
#       [--keep-chunks]
#
# Dependencies: seq2squiggle, slow5tools (for merge), GNU parallel or xargs
# =============================================================================
set -euo pipefail

# ---------------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------------
MANIFEST=""
OUTPUT=""
PROFILE="dna-r10-prom"
PARALLEL_JOBS=4
EXTRA_ARGS=""
KEEP_CHUNKS=false

# ---------------------------------------------------------------------------
# Argument parsing
# ---------------------------------------------------------------------------
usage() { grep '^#' "$0" | grep -v '#!/' | sed 's/^# \?//'; exit 1; }

while [[ $# -gt 0 ]]; do
    case $1 in
        --manifest)           MANIFEST="$2";        shift 2 ;;
        --output)             OUTPUT="$2";           shift 2 ;;
        --profile)            PROFILE="$2";          shift 2 ;;
        --parallel-jobs)      PARALLEL_JOBS="$2";    shift 2 ;;
        --seq2squiggle-args)  EXTRA_ARGS="$2";       shift 2 ;;
        --keep-chunks)        KEEP_CHUNKS=true;      shift ;;
        -h|--help)            usage ;;
        *) echo "Unknown argument: $1"; usage ;;
    esac
done

[[ -z "$MANIFEST" ]] && { echo "ERROR: --manifest is required"; exit 1; }
[[ -z "$OUTPUT"   ]] && { echo "ERROR: --output is required";   exit 1; }
[[ -f "$MANIFEST" ]] || { echo "ERROR: manifest not found: $MANIFEST"; exit 1; }

for tool in seq2squiggle slow5tools; do
    command -v "$tool" &>/dev/null || { echo "ERROR: $tool not found on PATH"; exit 1; }
done

CHUNK_DIR="$(dirname "$MANIFEST")"
BLOW5_DIR="$CHUNK_DIR/blow5_chunks"
LOG_DIR="$CHUNK_DIR/logs"
mkdir -p "$BLOW5_DIR" "$LOG_DIR"

# ---------------------------------------------------------------------------
# Read manifest: chunk_id, n_reads, fasta_path
# ---------------------------------------------------------------------------
mapfile -t FASTA_PATHS < <(tail -n +2 "$MANIFEST" | tr -d '\r' | awk -F'\t' '{print $3}')
N_CHUNKS="${#FASTA_PATHS[@]}"

echo "========================================"
echo " seq2squiggle chunk runner"
echo " Chunks    : $N_CHUNKS"
echo " Profile   : $PROFILE"
echo " Parallel  : $PARALLEL_JOBS jobs"
echo " Output    : $OUTPUT"
echo "========================================"
echo ""

# ---------------------------------------------------------------------------
# Function: process one chunk
# ---------------------------------------------------------------------------
process_chunk() {
    local fasta_path="$1"
    local profile="$2"
    local extra_args="$3"
    local blow5_dir="$4"
    local log_dir="$5"

    local chunk_name
    chunk_name="$(basename "$fasta_path" .fasta)"
    local blow5_out="$blow5_dir/${chunk_name}.blow5"
    local log_out="$log_dir/${chunk_name}.log"

    # Skip if already done (allows resuming after partial failure)
    if [[ -f "$blow5_out" ]]; then
        echo "  [SKIP] $chunk_name — BLOW5 already exists"
        return 0
    fi

    echo "  [RUN ] $chunk_name ($(wc -l < "$fasta_path" | awk '{print $1/2}') reads)"

    # shellcheck disable=SC2086
    if seq2squiggle predict \
            --read-input \
            --profile "$profile" \
            --preserve-read-ids \
            $extra_args \
            "$fasta_path" \
            -o "$blow5_out" \
            >"$log_out" 2>&1; then
        echo "  [OK  ] $chunk_name → $blow5_out"
    else
        echo "  [FAIL] $chunk_name — see $log_out"
        # Retry once
        echo "  [RTRY] $chunk_name ..."
        if seq2squiggle predict \
                --read-input \
                --profile "$profile" \
                --preserve-read-ids \
                $extra_args \
                "$fasta_path" \
                -o "$blow5_out" \
                >>"$log_out" 2>&1; then
            echo "  [OK  ] $chunk_name (retry succeeded)"
        else
            echo "  [ERR ] $chunk_name failed after retry — skipping"
            return 1
        fi
    fi
}

export -f process_chunk

# ---------------------------------------------------------------------------
# Run chunks — use GNU parallel if available, otherwise xargs
# ---------------------------------------------------------------------------
echo "Processing $N_CHUNKS chunks with $PARALLEL_JOBS parallel jobs..."
echo ""

FAILED=0

if command -v parallel &>/dev/null; then
    # GNU parallel: print progress, keep going on failure
    printf '%s\n' "${FASTA_PATHS[@]}" \
        | parallel --jobs "$PARALLEL_JOBS" --bar --joblog "$LOG_DIR/parallel.log" \
            process_chunk {} "$PROFILE" "$EXTRA_ARGS" "$BLOW5_DIR" "$LOG_DIR" \
        || FAILED=1
else
    # Fallback: xargs (less informative but no extra dependency)
    printf '%s\n' "${FASTA_PATHS[@]}" \
        | xargs -P "$PARALLEL_JOBS" -I{} \
            bash -c 'process_chunk "$@"' _ {} "$PROFILE" "$EXTRA_ARGS" "$BLOW5_DIR" "$LOG_DIR" \
        || FAILED=1
fi

# ---------------------------------------------------------------------------
# Check all chunks produced a BLOW5
# ---------------------------------------------------------------------------
echo ""
echo "Checking chunk outputs..."
MISSING=()
for fasta_path in "${FASTA_PATHS[@]}"; do
    chunk_name="$(basename "$fasta_path" .fasta)"
    blow5_out="$BLOW5_DIR/${chunk_name}.blow5"
    if [[ ! -f "$blow5_out" ]]; then
        MISSING+=("$chunk_name")
    fi
done

if [[ ${#MISSING[@]} -gt 0 ]]; then
    echo "WARNING: ${#MISSING[@]} chunk(s) failed to produce a BLOW5:"
    printf '  %s\n' "${MISSING[@]}"
    echo "These will be absent from the merged output."
fi

# ---------------------------------------------------------------------------
# Merge all chunk BLOW5s into the final output
# slow5tools merge takes a directory or list of files
# ---------------------------------------------------------------------------
echo ""
echo "Merging $(ls "$BLOW5_DIR"/*.blow5 2>/dev/null | wc -l) BLOW5 chunks → $OUTPUT ..."

# Collect all chunk BLOW5s in manifest order for deterministic output
BLOW5_LIST=()
for fasta_path in "${FASTA_PATHS[@]}"; do
    chunk_name="$(basename "$fasta_path" .fasta)"
    blow5_out="$BLOW5_DIR/${chunk_name}.blow5"
    [[ -f "$blow5_out" ]] && BLOW5_LIST+=("$blow5_out")
done

if [[ ${#BLOW5_LIST[@]} -eq 0 ]]; then
    echo "ERROR: No BLOW5 chunks to merge — all chunks failed."
    exit 1
fi

slow5tools merge \
    --output "$OUTPUT" \
    --threads "$(nproc)" \
    "${BLOW5_LIST[@]}"

# Verify output
N_READS_OUT=$(slow5tools stats "$OUTPUT" | grep "^number of reads" | awk '{print $NF}')
echo ""
echo "========================================"
echo " Merge complete"
echo " Output          : $OUTPUT"
echo " Reads in output : $N_READS_OUT"
echo " Chunks failed   : ${#MISSING[@]}"
echo "========================================"

# ---------------------------------------------------------------------------
# Cleanup chunk BLOW5s unless --keep-chunks
# ---------------------------------------------------------------------------
if [[ "$KEEP_CHUNKS" == false ]]; then
    echo ""
    echo "Removing chunk BLOW5s (use --keep-chunks to retain)..."
    rm -f "${BLOW5_LIST[@]}"
    echo "Done."
fi

[[ "$FAILED" -eq 1 ]] && exit 1 || exit 0