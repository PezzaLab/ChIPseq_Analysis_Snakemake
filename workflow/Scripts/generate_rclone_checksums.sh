#!/usr/bin/env bash
# ==============================================================================
# generate_rclone_checksums.sh
#
# Generates a checksum manifest (checksums.tsv) for all deliverables under
# Results/ to be exported to Dropbox via rclone, conforming to the
# data-tracking specifications from dialog 32cfd667-47a8-4121-88cf-f7cc23f687fd.
#
# Format:
# relative_path\thash\tsize_bytes\tmodified_date
#
# Compatibility: macOS BSD and Linux GNU (Bash 3.2+)
# ==============================================================================

set -euo pipefail

TARGET_DIR="${1:-.}"
OUTPUT_FILE="${2:-${TARGET_DIR}/Results/checksums.tsv}"

# Resolve target directory to absolute path
TARGET_DIR="$(cd "$TARGET_DIR" && pwd)"

# If output file is relative, make it relative to current working dir
if [[ "$OUTPUT_FILE" != /* ]]; then
    OUTPUT_FILE="$(pwd)/${OUTPUT_FILE}"
fi

OUTPUT_DIR="$(dirname "$OUTPUT_FILE")"
mkdir -p "$OUTPUT_DIR"

# Verify rclone is available
if ! command -v rclone >/dev/null 2>&1; then
    echo "[ERROR] 'rclone' command not found in PATH or environment modules." >&2
    echo "  Please ensure rclone is loaded (e.g. 'module load rclone' or 'ml rclone')." >&2
    exit 1
fi

RESULTS_DIR="${TARGET_DIR}/Results"
if [[ ! -d "$RESULTS_DIR" ]]; then
    echo "[WARN] Results directory not found at: $RESULTS_DIR" >&2
    echo "  Creating an empty checksum manifest." >&2
    printf "relative_path\thash\tsize_bytes\tmodified_date\n" > "$OUTPUT_FILE"
    touch "$OUTPUT_FILE"
    exit 0
fi

TMP_OUT="$(mktemp "${OUTPUT_DIR}/checksums.tmp.XXXXXX")"
printf "relative_path\thash\tsize_bytes\tmodified_date\n" > "$TMP_OUT"

echo "Generating checksums for exported deliverables under Results/..."

# Collect and hash deliverables strictly under Results/:
# 1. General files <= 50MB matching the rclone backup filters from commands.sh
# 2. All *.bw BigWig coverage files under Results/ (regardless of size)
(
    # Set 1: General deliverables <= 50MB under Results/
    rclone lsf "$TARGET_DIR" -R --files-only --format "phst" --hash md5 --separator $'\t' \
        --max-size 50M \
        --filter '- .git/**' \
        --filter '- .*' \
        --filter '- .*/' \
        --filter '- ~*' \
        --filter '- *fastp.html' \
        --filter '- not_copy*' \
        --filter '- *.bai' \
        --filter '- *.bam' \
        --filter '- *.bw' \
        --filter '- *.filename' \
        --filter '- *.matrix' \
        --filter '- *.fq.gz' \
        --filter '- *.fastq.gz' \
        --filter '- *.sra' \
        --filter '- *.homer_anotated.tsv' \
        --filter '- *.gappedPeak' \
        --filter '- Intersect_HSs_plus_minus_2000_bp/' \
        --filter '+ Results/**blacklist**' \
        --filter '- Peaks/**' \
        --filter '- *checksums*' \
        --filter '+ Results/**' \
        --filter '- *'

    # Set 2: All BigWig coverage files under Results/
    rclone lsf "$TARGET_DIR" -R --files-only --format "phst" --hash md5 --separator $'\t' \
        --filter '+ Results/**.bw' \
        --filter '- *'
) | sort -t$'\t' -u -k1,1 >> "$TMP_OUT"

# Atomically replace final output
mv "$TMP_OUT" "$OUTPUT_FILE"

# Ensure mtime is current (satisfies Stale Checksum Guard from dialog 32cfd667-47a8-4121-88cf-f7cc23f687fd)
touch "$OUTPUT_FILE"

TOTAL_COUNT=$(($(wc -l < "$OUTPUT_FILE" | tr -d ' ') - 1))
echo "[SUCCESS] Generated $OUTPUT_FILE ($TOTAL_COUNT deliverables indexed)"
