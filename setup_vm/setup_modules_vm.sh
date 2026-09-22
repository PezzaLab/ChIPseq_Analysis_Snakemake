#!/usr/bin/env bash
set -euo pipefail

# ==============================================================================
# Lmod Module Generator for Linux VM
# Creates modulefiles that mimic HPC environment modules for Snakemake
# ==============================================================================

MODULE_BASE_DIR="${HOME}/modulefiles"
mkdir -p "${MODULE_BASE_DIR}"

echo ">>> Creating module definitions in ${MODULE_BASE_DIR}..."

# List of modules defined in Config/config.yaml:
# bamutil: "bamutil/1.0.15"
# bedtools: "bedtools/2.30.0"
# bioconductor: "bioconductor/3.14"
# bwa: "bwa/0.7.15"
# cutadapt: "cutadapt/3.7"
# deeptools: "deeptools/3.4.3"
# fastp: "fastp/0.23.2"
# homer: "homer/5.1"
# macs2: "macs2/2.2.7.1"
# multiqc: "multiqc/1.15"
# phantompeakqualtools: "phantompeakqualtools/1.2"
# picard: "picard/2.21.2"
# R: "R/4.1.2-mkl"
# samtools: "samtools/1.14"
# sambamba: "sambamba/0.8.2"

create_module() {
    local mod_name="$1"
    local mod_version="$2"
    local bin_path="$3"
    local target_dir="${MODULE_BASE_DIR}/${mod_name}"
    mkdir -p "${target_dir}"

    cat <<MOD_EOF > "${target_dir}/${mod_version}.lua"
-- Lmod modulefile for ${mod_name}/${mod_version}
whatis("Loads ${mod_name} version ${mod_version} for local pipeline emulation")
prepend_path("PATH", "${bin_path}")
MOD_EOF
    echo "  [+] Created module ${mod_name}/${mod_version}"
}

# Point to conda environment or system bin directory
# Default to /usr/bin or conda env path
DEFAULT_BIN="/usr/local/bin"

create_module "bamutil" "1.0.15" "${DEFAULT_BIN}"
create_module "bedtools" "2.30.0" "${DEFAULT_BIN}"
create_module "bioconductor" "3.14" "${DEFAULT_BIN}"
create_module "bwa" "0.7.15" "${DEFAULT_BIN}"
create_module "cutadapt" "3.7" "${DEFAULT_BIN}"
create_module "deeptools" "3.4.3" "${DEFAULT_BIN}"
create_module "fastp" "0.23.2" "${DEFAULT_BIN}"
create_module "homer" "5.1" "${DEFAULT_BIN}"
create_module "macs2" "2.2.7.1" "${DEFAULT_BIN}"
create_module "multiqc" "1.15" "${DEFAULT_BIN}"
create_module "phantompeakqualtools" "1.2" "${DEFAULT_BIN}"
create_module "picard" "2.21.2" "${DEFAULT_BIN}"
create_module "R" "4.1.2-mkl" "${DEFAULT_BIN}"
create_module "samtools" "1.14" "${DEFAULT_BIN}"
create_module "sambamba" "0.8.2" "${DEFAULT_BIN}"

echo ""
echo "=============================================================================="
echo "Modulefiles created successfully in: ${MODULE_BASE_DIR}"
echo "To make them available in Lmod, add this line to your ~/.bashrc or run:"
echo "    module use ${MODULE_BASE_DIR}"
echo "You can then verify with: module avail"
echo "=============================================================================="
