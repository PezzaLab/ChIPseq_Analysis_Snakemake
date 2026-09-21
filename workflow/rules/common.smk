import os
import re
import sys
import pandas as pd
import numpy as np
import smk_functions as smkf
import suffixes as sfxs

#########################################
## Load sample table and check format ###
#########################################
samples_table_raw = pd.read_csv(
    config.get("samples", "Config/samples.csv"),
    comment='#',
    dtype=str,
    na_values={
        "merge_with": "-",
        "fastq2": "-",
        "peak_ctrl_file_alias": "-",
    },
    keep_default_na=True
)
samples_table_raw.columns = samples_table_raw.columns.str.strip()

samples_table_no_merged_samples = smkf.expand_sample_table(samples_table_raw)

smkf.config = config
smkf.check_sample_table_format(samples_table_no_merged_samples)

# Pre formed strings
cov_config_params_string=(
    f"{config['coverage']['normalization']}_"
    f"bs{config['coverage']['bin_size']}_sm{config['coverage']['smooth']}"
    f"ex{config['coverage']['extend_reads']}"
)
dros_norm_cov_config_params_string=(
    f"drosNormalized_bs{config['coverage']['bin_size']}_sm"
    f"{config['coverage']['smooth']}"
    f"ex{config['coverage']['extend_reads']}"
)

############################################################
########## Add merged samples to samples table #############
############################################################
new_rows = samples_table_no_merged_samples.loc[
    pd.notna(samples_table_no_merged_samples['merge_with'])
].copy()

if len(new_rows) > 0:
    new_rows_names = new_rows['sample_name'] + "_MERGED"
    new_rows['sample_name'] = new_rows_names
    new_rows['merge_with'] = np.nan
    new_rows_index = new_rows['sample_name'] + "." + new_rows['reference_genome']
    new_rows.index = new_rows_index
    new_rows['sample_entry_id'] = new_rows_index
    samples_table = pd.concat([samples_table_no_merged_samples, new_rows], verify_integrity=True)
else:
    samples_table = samples_table_no_merged_samples.copy()

####################################################
##### Add samples names with suffixes to table #####
####################################################
samples_table_2 = sfxs.generate_samples_table_2(samples_table, config)
samples_table_2.to_csv("Config/samples_table_processed.csv", index=False)

############################################
#####   Add variables to smkf module   #####
############################################
smkf.samples_table_no_merged_samples = samples_table_no_merged_samples
smkf.samples_table = samples_table
smkf.config = config
smkf.samples_table_2 = samples_table_2

# Determine if running in basic or full mode
# If running through basic.smk, or if 'basic' is targeted on CLI without 'full'
is_basic_mode = (
    config.get("mode") == "basic" or 
    ("basic" in sys.argv and "full" not in sys.argv)
)
smkf.include_hotspots = not is_basic_mode

#####################################################
#####        Target output collections          #####
#####################################################
# MultiQC report
multiqc = []
if config['qctrl']:
    multiqc = expand(
        "Results/{genomes}/Qctrl/multiqc_report_{lib_name}.html",
        genomes = samples_table["reference_genome"].unique().tolist(),
        lib_name = config['library']['name']
    )

# Basic unstranded BigWig coverage files (Both strands CPM + Drosophila normalized)
basic_coverage_bigwigs = [
    a for a in samples_table_2.filter(regex = r"(?i).*both_strands_coverage_bw$").values.flatten().tolist()
    if pd.notna(a)
]

# Single strand BigWig coverage files (Watson, Crick, clipped, etc.)
single_strand_coverage_bigwigs = [
    a for a in samples_table_2.filter(regex = ".*coverage_bw$").values.flatten().tolist()
    if a == a and a not in basic_coverage_bigwigs
]

# All coverage BigWigs combined
all_coverage_bigwigs = basic_coverage_bigwigs + single_strand_coverage_bigwigs

# FRIP score files (calculated on broad peaks for samples with controls)
broad_frip_scores = [
    f for f in samples_table_2['broad_blk_gr_flt_FRIP'].values.tolist()
    if pd.notna(f)
] if 'broad_blk_gr_flt_FRIP' in samples_table_2 else []

# Basic blacklist-filtered peaks (narrow and broad) + broad FRIP scores
basic_peaks = [
    p for p in samples_table_2.filter(regex = r"^(narrow|broad)_peak_bl_gr_flt$").values.flatten().tolist()
    if pd.notna(p)
] + broad_frip_scores

# Peak count summary files for genomes that have peak controls
peaks_summary_files = [
    f"Results/{genome}/Analysis/Peaks_summary.tsv"
    for genome in samples_table["reference_genome"].unique()
    if (
        pd.notna(samples_table.loc[samples_table["reference_genome"] == genome, 'peak_ctrl_file_alias']) &
        (samples_table.loc[samples_table["reference_genome"] == genome, 'peak_ctrl_file_alias'] != "-")
    ).any()
]

# Advanced targets
heatmaps = samples_table_2.filter(regex="^heatmap.*").values.flatten().tolist()
heatmaps = [x for x in heatmaps if x == x]

html_reports = expand(
    "Results/{genomes}/Analysis/{lib_name}.{genomes}.{smooth}.html",
    genomes = samples_table["reference_genome"].unique().tolist(),
    lib_name = config['library']['name'],
    smooth = ("smoothed" if config['aggregate_profiles']['process_outfile']['smooth'] else "not_smoothed"),
)

#####################################################
#####           Wildcard constraints            #####
#####################################################
wildcard_constraints:
    bam_type = "Raw_bam|Processed_bam",
    bin_size = "[0-9]+",
    clipped_1bp_PE=r"((inc|exc)_16\.clipped_1_bp\.)?",
    clip_condit = r"(((inc|exc)_16\.)?clipped_1_bp\.)?",
    clip_strand = r"((inc_|exc)_16\.)?",
    cov_params ="((CPM|RPKM|none|drosNormalized)_)?bs[0-9]+_sm[0-9]+_ex[0-9]+",
    extension = r"(\.q_filt\.srt\.nodup\.mit_filt)?",
    extension_size = "[0-9]+",
    flag = "(83|163|99|147|inc_16|exc_16)",
    fused_genome="(mm10_f_d6|mm39_f_d6|hg19_f_d6|hg38_f_d6|mm10_x_CAST_EiJ_f_d6|mm39_x_CAST_EiJ)",
    genomes_not_fused="mm10|mm39|d6|hg19|hg38|mm10_x_CAST_EiJ",
    genomes_final = (
        r"((mm10|mm39|hg19|hg38|mm10_x_CAST_EiJ|mm39_x_CAST_EiJ)"
        r"(_f_d6\.(mm10|mm39|hg19|hg38|mm10_x_CAST_EiJ|mm39_x_CAST_EiJ|d6))?)"
    ),
    genomes_all=(
        "mm10|mm39|d6|hg19|hg38|mm10_f_d6|mm39_f_d6|hg19_f_d6|hg38_f_d6|mm10_x_CAST_EiJ_f_d6|"
        "mm10_x_CAST_EiJ"),
    hs_region = (
        "B6xCAST_(top_5000_pm_2000bp|"
        "PRDM9_assymetric_hs_("
                             "(invading|receiving)_strand|"
                             "mm10_aligned))|"
        "top_5000_plus_minus_2000|asymetric_(watson|crick)_strong|x_non_par|"
        "autosomal_x_non_par_ctrl"),
    int_genome = "mm10|mm39|hg19|hg38|mm10_x_CAST_EiJ",
    library_name = config['library']['name'],
    norm = "[^_/.]+",
    n_ends = "(SE|PE)",
    peak_params = "(bco_[0-9]+_)?qv_[0-9]+__[^/]*",
    peak_type = "narrow|broad",
    read = "[12]{1}",
    read_length = "Full_length_reads|1bp_clipped_reads",
    sample = "[^./ ]+",
    sample_with_extensions = "[^/]+",
    smooth = "(smoothed|not_smoothed)+",
    smooth_size = "[0-9]+",
    spike_genome="d6",
    ss_condit = r"((83-163|99-147|inc_16|exc_16)\.)?",
    ss_PE = "((83-163|99-147).)?",
    ss_SR = "(inc|exc)_16",
    strand = r"(\.(83-163|99-147|inc_16|exc_16))?",
    strands = "Both_strands|Single_strand/(1bp_clipped_reads|Full_length_reads)",
