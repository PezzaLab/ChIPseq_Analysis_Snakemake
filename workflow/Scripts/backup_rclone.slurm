#!/usr/bin/env bash
#SBATCH --job-name=rclone_backup
#SBATCH --account=pezzar-lab
#SBATCH --partition=serial
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=04:00:00
#SBATCH --output=.snakemake/slurm_logs/rclone_backup_%j.log
#SBATCH --error=.snakemake/slurm_logs/rclone_backup_%j.log
#SBATCH --mail-type=FAIL

# ==============================================================================
# backup_rclone.slurm
#
# Slurm batch script to transfer library deliverables and BigWigs to Dropbox
# via rclone.
#
# Defaults:
#   Source:      Current directory ($PWD)
#   Destination: dropboxOMRF:Bioinformatics/Libraries/<library_name>
#
# Options:
#   --non-library <path>  Specify an alternative subpath after 'Bioinformatics/'
#                         instead of the default 'Libraries'.
#                         Example: --non-library OtherProjects/ProjectA
#   --dest <remote>       Explicitly specify the full destination remote path.
#   --dry-run             Run rclone in dry-run mode (no files copied).
#   -h, --help            Show usage information.
# ==============================================================================

set -euo pipefail

SRC_DIR=""
DEST_REMOTE=""
NON_LIBRARY=""
DRY_RUN_FLAG=""

# Parse command line options
while [[ $# -gt 0 ]]; do
    case "$1" in
        --non-library)
            if [[ -z "${2:-}" || "${2:-}" == --* ]]; then
                echo "[ERROR] '--non-library' requires a path argument (path after Bioinformatics/)." >&2
                exit 1
            fi
            NON_LIBRARY="$2"
            shift 2
            ;;
        --non-library=*)
            NON_LIBRARY="${1#*=}"
            shift 1
            ;;
        --dest)
            if [[ -z "${2:-}" || "${2:-}" == --* ]]; then
                echo "[ERROR] '--dest' requires a remote path argument." >&2
                exit 1
            fi
            DEST_REMOTE="$2"
            shift 2
            ;;
        --dest=*)
            DEST_REMOTE="${1#*=}"
            shift 1
            ;;
        --dry-run)
            DRY_RUN_FLAG="--dry-run"
            shift 1
            ;;
        -h|--help)
            echo "Usage: sbatch workflow/Scripts/backup_rclone.slurm [options] [SRC_DIR] [DEST_REMOTE]"
            echo ""
            echo "Options:"
            echo "  --non-library <path>  Specify alternate subpath after 'Bioinformatics/'"
            echo "                        (Default: 'Libraries', targeting 'dropboxOMRF:Bioinformatics/Libraries/<lib_name>')"
            echo "  --dest <remote>       Explicitly set full destination remote path (overrides --non-library)"
            echo "  --dry-run             Perform trial run without transferring files"
            echo "  -h, --help            Show this help message and exit"
            echo ""
            echo "Positional Arguments:"
            echo "  SRC_DIR               Source directory (default: current working directory)"
            echo "  DEST_REMOTE           Full destination remote (optional, overrides default/--non-library)"
            exit 0
            ;;
        *)
            if [[ -z "$SRC_DIR" ]]; then
                SRC_DIR="$1"
            elif [[ -z "$DEST_REMOTE" ]]; then
                DEST_REMOTE="$1"
            else
                echo "[ERROR] Unexpected positional argument: $1" >&2
                exit 1
            fi
            shift 1
            ;;
    esac
done

# Resolve source directory
if [[ -z "$SRC_DIR" ]]; then
    SRC_DIR="$(pwd)"
fi

if [[ ! -d "$SRC_DIR" ]]; then
    echo "[ERROR] Source directory does not exist: $SRC_DIR" >&2
    exit 1
fi

SRC_DIR="$(cd "$SRC_DIR" && pwd)"
LIB_NAME="$(basename "$SRC_DIR")"

# Resolve destination remote
if [[ -z "$DEST_REMOTE" ]]; then
    if [[ -n "$NON_LIBRARY" ]]; then
        # Strip leading and trailing slashes
        CLEAN_SUBPATH="$(echo "$NON_LIBRARY" | sed -e 's|^/*||' -e 's|/*$||')"
        if [[ "$(basename "$CLEAN_SUBPATH")" == "$LIB_NAME" ]]; then
            DEST_REMOTE="dropboxOMRF:Bioinformatics/${CLEAN_SUBPATH}"
        else
            DEST_REMOTE="dropboxOMRF:Bioinformatics/${CLEAN_SUBPATH}/${LIB_NAME}"
        fi
    else
        DEST_REMOTE="dropboxOMRF:Bioinformatics/Libraries/${LIB_NAME}"
    fi
fi

echo "================================================================================"
echo ">>> SLURM RCLONE BACKUP"
echo "  Date:        $(date)"
echo "  Source:      ${SRC_DIR}"
echo "  Destination: ${DEST_REMOTE}"
if [[ -n "$NON_LIBRARY" ]]; then
    echo "  Subpath:     ${NON_LIBRARY} (custom via --non-library)"
else
    echo "  Subpath:     Libraries (default)"
fi
if [[ -n "$DRY_RUN_FLAG" ]]; then
    echo "  Mode:        DRY RUN (simulation only)"
fi
echo "================================================================================"

# Verify and load rclone module
if ! command -v rclone >/dev/null 2>&1; then
    if [ -f /etc/profile.d/modules.sh ]; then
        source /etc/profile.d/modules.sh
    elif [ -f /etc/profile.d/lmod.sh ]; then
        source /etc/profile.d/lmod.sh
    fi
    if command -v ml >/dev/null 2>&1; then
        echo ">>> Loading rclone module via ml..."
        module purge 2>/dev/null || true
        ml rclone
    elif command -v module >/dev/null 2>&1; then
        echo ">>> Loading rclone module via module load..."
        module purge 2>/dev/null || true
        module load rclone
    fi
fi

if ! command -v rclone >/dev/null 2>&1; then
    echo "[ERROR] 'rclone' command not found in PATH or environment modules." >&2
    exit 1
fi

# Ensure target slurm log directory exists
mkdir -p "${SRC_DIR}/.snakemake/slurm_logs"

# 1) Backup deliverables and metadata (<= 50M)
echo ">>> [1/2] Backing up deliverables (<= 50M)..."
rclone copy ${DRY_RUN_FLAG} "${SRC_DIR}" "${DEST_REMOTE}" \
    --max-size 50M \
    --stats 1m \
    --filter '+ .snakemake/slurm_logs/**' \
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
    --filter '+ **blacklist**' \
    --filter '- Peaks/**'

# 2) Backup BigWig coverage files
echo ">>> [2/2] Backing up BigWig (*.bw) files..."
rclone copy ${DRY_RUN_FLAG} "${SRC_DIR}" "${DEST_REMOTE}" \
    --include '*.bw' \
    --stats 1m

echo "================================================================================"
echo ">>> Rclone backup completed successfully at: $(date)"
echo "================================================================================"
