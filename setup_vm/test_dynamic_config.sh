#!/usr/bin/env bash
set -euo pipefail

REPO_DIR="${HOME}/projects/ChIPseq_snakemake_pipeline"
cd "${REPO_DIR}"
SNAKEMAKE="${HOME}/miniforge3/envs/snakemake_target/bin/snakemake"

echo "=============================================================================="
echo ">>> DYNAMIC CONFIGURATION & FILENAME VERIFICATION SUITE"
echo "=============================================================================="

# Backup original files
cp Config/config.yaml Config/config.yaml.test_backup
cp Config/samples.csv Config/samples.csv.test_backup

cleanup() {
    echo ""
    echo ">>> Restoring original config and samples files..."
    mv Config/config.yaml.test_backup Config/config.yaml
    mv Config/samples.csv.test_backup Config/samples.csv
}
trap cleanup EXIT

# ------------------------------------------------------------------------------
# TEST 1: Baseline Check
# ------------------------------------------------------------------------------
echo ""
echo "--- [TEST 1] Baseline Configuration Verification ---"
${SNAKEMAKE} --profile Config/Profiles/slurm_quio_repeat_3 -np full > /tmp/test1_dryrun.txt 2>&1

check_test1() {
    local pass=true
    if grep -q "CPM_bs4_sm1_ex250" /tmp/test1_dryrun.txt; then
        echo "  [PASS] Found default coverage directory: CPM_bs4_sm1_ex250"
    else
        echo "  [FAIL] Missing coverage directory CPM_bs4_sm1_ex250"; pass=false
    fi

    if grep -q "qv_05__no_input" /tmp/test1_dryrun.txt && grep -q "bco_1_qv_05__no_input" /tmp/test1_dryrun.txt; then
        echo "  [PASS] Found default MACS2 directories: qv_05__no_input & bco_1_qv_05__no_input"
    else
        echo "  [FAIL] Missing default MACS2 directories"; pass=false
    fi

    if grep -q "multiqc_report_Test_library.html" /tmp/test1_dryrun.txt; then
        echo "  [PASS] Found default MultiQC report: multiqc_report_Test_library.html"
    else
        echo "  [FAIL] Missing MultiQC report for Test_library"; pass=false
    fi

    if grep -q "Test_library_aggregate_profiles_data.smoothed.RData" /tmp/test1_dryrun.txt; then
        echo "  [PASS] Found default aggregate profile: Test_library_aggregate_profiles_data.smoothed.RData"
    else
        echo "  [FAIL] Missing default aggregate profile"; pass=false
    fi
    $pass
}
check_test1

# ------------------------------------------------------------------------------
# TEST 2: Coverage & Library Name Mutation
# Mutate:
#   library.name -> "AlphaRun"
#   coverage.normalization -> "RPKM"
#   coverage.bin_size -> "10"
#   coverage.smooth -> "3"
#   coverage.extend_reads -> "180"
# ------------------------------------------------------------------------------
echo ""
echo "--- [TEST 2] Mutating Library Name & Coverage Parameters ---"
sed -i 's/name: "Test_library"/name: "AlphaRun"/' Config/config.yaml
sed -i 's/normalization: "CPM"/normalization: "RPKM"/' Config/config.yaml
sed -i 's/bin_size: "4"/bin_size: "10"/' Config/config.yaml
sed -i 's/smooth: "1"/smooth: "3"/' Config/config.yaml
sed -i 's/extend_reads: "250"/extend_reads: "180"/' Config/config.yaml

${SNAKEMAKE} --profile Config/Profiles/slurm_quio_repeat_3 -np full > /tmp/test2_dryrun.txt 2>&1

check_test2() {
    local pass=true
    if grep -q "RPKM_bs10_sm3_ex180" /tmp/test2_dryrun.txt; then
        echo "  [PASS] BigWigs dynamically moved to directory: RPKM_bs10_sm3_ex180"
    else
        echo "  [FAIL] BigWigs did NOT move to RPKM_bs10_sm3_ex180"; pass=false
    fi

    if grep -q "RPKM_bs10_sm3_ex180/matrixes" /tmp/test2_dryrun.txt; then
        echo "  [PASS] Deeptools matrix inputs dynamically adapted to: RPKM_bs10_sm3_ex180"
    else
        echo "  [FAIL] Deeptools matrix path did NOT adapt"; pass=false
    fi

    if grep -q "multiqc_report_AlphaRun.html" /tmp/test2_dryrun.txt; then
        echo "  [PASS] MultiQC report filename adapted to: multiqc_report_AlphaRun.html"
    else
        echo "  [FAIL] MultiQC report filename did NOT adapt"; pass=false
    fi

    if grep -q "AlphaRun_aggregate_profiles_data.smoothed.RData" /tmp/test2_dryrun.txt; then
        echo "  [PASS] Aggregate profile RData adapted to: AlphaRun_aggregate_profiles_data.smoothed.RData"
    else
        echo "  [FAIL] Aggregate profile RData did NOT adapt"; pass=false
    fi

    if grep -q "AlphaRun.mm10.smoothed.html" /tmp/test2_dryrun.txt; then
        echo "  [PASS] Final HTML report adapted to: AlphaRun.mm10.smoothed.html"
    else
        echo "  [FAIL] Final HTML report did NOT adapt"; pass=false
    fi
    $pass
}
check_test2

# ------------------------------------------------------------------------------
# TEST 3: MACS2 Q-value & Broad Cutoff Mutation
# Mutate:
#   MACS2.qvalue -> "0.01"  (directory tag qv_01)
#   MACS2.broad_cutoff -> "0.05" (directory tag bco_05)
# ------------------------------------------------------------------------------
echo ""
echo "--- [TEST 3] Mutating MACS2 Peak Calling Cutoffs ---"
sed -i 's/qvalue: "0.05"/qvalue: "0.01"/' Config/config.yaml
sed -i 's/broad_cutoff: "0.1"/broad_cutoff: "0.05"/' Config/config.yaml

${SNAKEMAKE} --profile Config/Profiles/slurm_quio_repeat_3 -np full > /tmp/test3_dryrun.txt 2>&1

check_test3() {
    local pass=true
    if grep -q "Peaks/MACS2/narrow/qv_01__no_input" /tmp/test3_dryrun.txt; then
        echo "  [PASS] Narrow peaks directory dynamically shifted to: qv_01__no_input"
    else
        echo "  [FAIL] Narrow peaks did NOT shift to qv_01__no_input"; pass=false
    fi

    if grep -q "Peaks/MACS2/broad/bco_05_qv_01__no_input" /tmp/test3_dryrun.txt; then
        echo "  [PASS] Broad peaks directory dynamically shifted to: bco_05_qv_01__no_input"
    else
        echo "  [FAIL] Broad peaks did NOT shift to bco_05_qv_01__no_input"; pass=false
    fi

    if grep -q "qv_01__no_input/blacklist_filtered" /tmp/test3_dryrun.txt; then
        echo "  [PASS] Downstream blacklist filtering correctly tracks new qv_01 directory"
    else
        echo "  [FAIL] Downstream blacklist filtering failed to track new directory"; pass=false
    fi
    $pass
}
check_test3

# ------------------------------------------------------------------------------
# TEST 4: Library Technology Mutation (Adaptase Trimming Verification)
# Mutate:
#   Bloom_10_PE_no_dros -> library_technology: adaptase
# ------------------------------------------------------------------------------
echo ""
echo "--- [TEST 4] Mutating Library Technology (regular -> adaptase) ---"
# Replace library_technology for Bloom_10_PE_no_dros from regular to adaptase
sed -i 's/Bloom_10_PE_no_dros,Resources\/test_FASTQs\/Bloom_10_HS_peaks_R1.no_dros.R1.fastq.gz,Resources\/test_FASTQs\/Bloom_10_HS_peaks_R1.no_dros.R2.fastq.gz,True,regular,/Bloom_10_PE_no_dros,Resources\/test_FASTQs\/Bloom_10_HS_peaks_R1.no_dros.R1.fastq.gz,Resources\/test_FASTQs\/Bloom_10_HS_peaks_R1.no_dros.R2.fastq.gz,True,adaptase,/' Config/samples.csv

${SNAKEMAKE} --profile Config/Profiles/slurm_quio_repeat_3 -np full > /tmp/test4_dryrun.txt 2>&1

check_test4() {
    local pass=true
    if grep -q "rule trim_adaptase:" /tmp/test4_dryrun.txt; then
        echo "  [PASS] Rule 'trim_adaptase' was dynamically triggered by adaptase technology"
    else
        echo "  [FAIL] Rule 'trim_adaptase' was NOT triggered"; pass=false
    fi

    if grep -A 5 "cutadapt" /tmp/test4_dryrun.txt | grep -q -- "-u 10"; then
        echo "  [PASS] Cutadapt command correctly planned with -u 10 for adaptase trimming"
    else
        echo "  [FAIL] Cutadapt command missing or incorrect"; pass=false
    fi

    if grep -q "Bloom_10_PE_no_dros.adap_trimmed.adaptase_trimmed.R1.PE.fq.gz" /tmp/test4_dryrun.txt; then
        echo "  [PASS] Alignment rule correctly consumes .adaptase_trimmed. FASTQ input"
    else
        echo "  [FAIL] Alignment rule did not consume adaptase_trimmed input"; pass=false
    fi
    $pass
}
check_test4

# ------------------------------------------------------------------------------
# TEST 5: Drosophila Spike-in Common Scale Verification
# Verify all spike-in samples normalize to 'drosophila_100K_reads'
# ------------------------------------------------------------------------------
echo ""
echo "--- [TEST 5] Drosophila Spike-in Common Scale Verification ---"
${SNAKEMAKE} --profile Config/Profiles/slurm_quio_repeat_3 -np full > /tmp/test5_dryrun.txt 2>&1

check_test5() {
    local pass=true
    if grep -q "drosophila_100K_reads/drosophila_equalization_report.tsv" /tmp/test5_dryrun.txt; then
        echo "  [PASS] Drosophila report path correctly uses unified: drosophila_100K_reads"
    else
        echo "  [FAIL] Drosophila report path missing drosophila_100K_reads"; pass=false
    fi

    if grep -q "drosophila_100K_reads/Bloom_10_PE_with_dros" /tmp/test5_dryrun.txt; then
        echo "  [PASS] Drosophila BigWig coverage correctly routed to: drosophila_100K_reads"
    else
        echo "  [FAIL] Drosophila BigWig missing drosophila_100K_reads"; pass=false
    fi
    $pass
}
check_test5

echo ""
echo "=============================================================================="
echo ">>> ALL 5 DYNAMIC CONFIGURATION TESTS PASSED CLEANLY!"
echo "=============================================================================="
