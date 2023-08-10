import pandas as pd
import numpy as np
import os
import re

# Functions
def check_sample_table_format(samples_table):
    exit_message = "\n\n One or more errors have been detected on your 'samples_table.xlsx' file.\nErrors:\n"
    exit_script = False
    # Check names of samples don't contain / or . or finish in "_MERGED"
    if (
    (samples_table['sample_name'].str.contains("\.")) |
    (samples_table['sample_name'].str.contains("/")) |
    (samples_table['sample_name'].str.contains(" ")) |
    (samples_table['sample_name'].str.contains("_MERGED$"))).any():
        exit_message += "* One or more sample names in samples table contain a dot " +\
        "('.'), a slash ('/'), a space (' ') or ends with the word '_MERGED'. These are not allowed " +\
        "in sample names. Modify names and try again\n"
        exit_script = True
    
    # Check that there is no sample_name:genome combination repeated
    if samples_table.loc[:,['sample_name', 'reference_genome']].duplicated().any():
        exit_message += "* One or more sample names in samples table is repeated.\n" +\
        "Choose different names for all your samples.\n"
        exit_script = True
    
    # Check that "merge_with" reference is not a deduplicated/filtered file for spiked-in samples (it has to be raw bam)
    if pd.notnull(samples_table['merge_with']).any():
        if (
            (samples_table['merge_with'].str.contains(".nodup.")) &
            (samples_table['dros_spike_in'])
        ).any():
            exit_message += "* One or more samples' 'merge_with' parameter is a deduplicated/filtered " +\
                "bam file AND has drosophila spike in. When sample is spiked, the 'merge_with' " +\
                "file has to be the raw bam.\n"
            exit_script = True

    # Check that there are 2 FASTQs when sample is PE and 1 when is not
    if (
    (samples_table['PE']) &
    ( (samples_table['fastq1'] == "") | (samples_table['fastq2'] == "") )).any():
        exit_message += "* At least one sample is set as PE but only contains one " +\
        "FASTQ path. Interleaved FASTQs are not supported yet.\n"
        exit_script = True
    if ((~samples_table['PE']) & (~samples_table['fastq2'].isna())).any():
        exit_message += "* At least one sample is set as SR (PE == False) but contains " +\
        "a FASTQ path at column 'fastq2'. Please put it at column 'fastq1'\n"
        exit_script = True
    if "" in samples_table.drop(["fastq1","fastq2"], axis=1).values:
        exit_message += "* There is a cell that  is empty. There cannot be any empty cell. Please fill empty cells with '-'\n"
        exit_script = True
    
    # Check that the samples with a "dros_equalization_group" have "dros_spike_in" == T 
    #(I don't do the complementary becuase you might have the sample to equalize on another library)
    if (~samples_table['dros_equalization_group'].isnull() &
    ~samples_table['dros_spike_in']).any():
        exit_message += "* At least one of your samples has a 'dros_equalization_group' " +\
        "assigned to it but has the field 'dros_spike_in' set as False.\n"
        exit_script = True
    
    # Check that all columns with false-true are actually false true (check if col == bool)
    booleans_check=samples_table[["PE","dros_spike_in",
                     "get_single_strand",
                     "Clip_reads_to_1bp_on_5_prime",
                     "top5000_HS_heatmap"]].dtypes
    if (booleans_check != "bool").any():
        exit_message += "* One or more of the following columns has at lesat one row" +\
        "filled with something different than 'T', 'F', 'True', 'False', 'TRUE' or 'FALSE'.\n" +\
        "columns: 'PE','dros_spike_in', 'get_single_strand', 'Clip_reads_to_1bp_on_5_prime'," +\
        "'top5000_HS_heatmap'\n"
        exit_script = True
    
    # Exit if any previous condition is met
    if exit_script:
        print(exit_message)
        quit()

##############################################################
############## Functions to get inputs/params ################
##############################################################
# Function to get input
    # Cannot return fastq2="" because then snakemake looks for "" file.
    # And if I return only fastq1, then in the shell I am searching for input.fastq2, which is not there
def align_fastq_input(w):
    inpt = {"genome_path":config['genomes'][w.genomes_all]}
    
    if samples_table.loc[w.sample, "PE"]:
        inpt|= {
        "fastq1":samples_table.loc[w.sample, "fastq1"],
        "fastq2":samples_table.loc[w.sample, "fastq2"],
        }
    else:
        inpt|= {
        "fastq1":samples_table.loc[w.sample, "fastq1"],
        "fastq2":[]} # Cannot use "" here because snakemake will look for "" file
        # and give a "missing input" error. Instead, I give it an empty list
    return inpt


# Function for input/ctrl for peak calling rule
def call_peaks_macs2_input(w):
    ctrl_alias = w.peak_params.split("__") # "bco_1_qv_05__alias"
    ctrl_alias = ctrl_alias[1] # "alias"
    if ctrl_alias == "no_input":
        dict={
        "treat_bam": samples_table_2.loc[w.sample, 'dedup_flt_both_strds_bam'],
        "treat_bai": samples_table_2.loc[w.sample, 'dedup_flt_both_strds_bam'] + ".bai"}
    else:
        dict={
        "treat_bam": samples_table_2.loc[w.sample, 'dedup_flt_both_strds_bam'],
        "treat_bai": samples_table_2.loc[w.sample, 'dedup_flt_both_strds_bam'] + ".bai",
        "ctrl_bam": config['MACS2']['control'][ctrl_alias],
        "ctrl_bai": config['MACS2']['control'][ctrl_alias] + ".bai"
        }
    return dict

# Get folder name of input get_rv_fw_strand_input
def get_rv_fw_strand_input(w):
    if re.match(".*(83-163|99-147)", w.sample_with_extensions):
        strand="Single_strand/Full_length_reads"
    else:
        strand="Both_strands"
    return f"Results/{w.genomes_not_fused}/Bams/{strand}/{w.sample_with_extensions}.bam"

# Get peak parameters
def get_MACS2_params(w):
    if re.match (".*mm.*", w.genomes_not_fused):
        genome = "mm"
    elif re.match (".*hg.*", w.genomes_not_fused):
        genome = "hs" # I could change this to an regex that extracts any of the genomes

    PE = ""
    ext = ""
    parameters = w['peak_params'].split("__")
    bco_qv = parameters[0].split("_")
    
    ctrl_alias = parameters[-1]
    
    ctrl_bam = config['MACS2']['control'][ctrl_alias]
    if parameters[-1] == "no_input":
        ctrl = ""
    else:
        ctrl = f"-c {ctrl_bam}"  
    
    if samples_table.loc[w.sample, "PE"]:
        PE = "--format BAMPE"
    else:
        ext = f"--extsize {config['MACS2']['extension']}"

    if w.peak_type == "narrow":
        peak_type = "--call-summits"
        q = f"-q 0.{bco_qv[1]}"
    else:
        peak_type = "--broad"
        q = f"-q 0.{bco_qv[3]} --broad-cutoff 0.{bco_qv[1]}"

    return {
        "genome": genome,
        "PE": PE,
        "extension": ext,
        "peak_type_options": peak_type,
        "qv_bco": q,
        "ctrl": ctrl}

def process_aggregate_profiles_inputs(w):
    # Available wildcards: 'genomes_not_fused' and 'hs_region'.
    # w.hs_region can be one of the following:
        # 'top_5000_plus_minus_2000', 'x_non_par', 'autosomal_x_non_par_ctrl', 
        # 'asymetric_watson_strong', 'asymetric_crick_strong', 'top_5000_plus_minus_2000'
        # 'B6xCAST_top_5000_pm_1000bp', 'B6xCAST_PRDM9_assymetric_hs_invading_strand'
        # or 'B6xCAST_PRDM9_assymetric_hs_receiving_strand'
    
    # This rule will process all matrixes that come from the same list of HS
    matrixes_df = samples_table_2.filter(
        regex=f"^{w.hs_region}.*matrix$"
        )
    matrixes_array = matrixes_df.to_numpy().ravel()
    matrixes_list = matrixes_array[~pd.isnull(matrixes_array)].tolist()
    return matrixes_list

def sumarize_peak_count_input(w):
    genome_filtered = samples_table_2[
        samples_table_2['reference_genome'] == w.genomes_not_fused
    ]
    peak_types=["narrow", "broad"]
    peaks = {}
    for peak_type in peak_types:
        peaks |= {
        f"{peak_type}_all": genome_filtered.loc[
            genome_filtered[f'{peak_type}_peak_bl_gr_flt'].notnull(),
            f'{peak_type}_peak_bl_gr_flt'].values.tolist(),
        f"{peak_type}_hs": genome_filtered.loc[
            genome_filtered[f'{peak_type}_peak_bl_gr_flt_hs_int'].notnull(),
            f'{peak_type}_peak_bl_gr_flt_hs_int'].values.tolist(),
        }
    return peaks

def get_qctrl_bams_bais(w):
    # This function gets input for rules 'samstats' and 'samtools_flagstat'
    sample = w.sample # If I don't convert the wildcards to set, it returns a warning when I use it within loc (next line)
    genome = samples_table.loc[sample, 'reference_genome']
    if f"{w.bam_type}" == 'Raw_bam':
        if w.genomes_not_fused != "d6":
            bam = samples_table_2.loc[sample, 'raw_bam']
        else:
            bam = f"Results/{w.sample}.{genome}_f_d6.d6.bam"
    else:
        if w.genomes_not_fused != "d6":
            bam = samples_table_2.loc[sample, 'dedup_flt_both_strds_bam']
        else:
            bam = f"Results/d6/Bams/Both_strands/{w.sample}.{genome}_f_d6.d6.q_filt.srt.nodup.mit_filt.bam"
    bai = f"{bam}.bai"
    return {
    "bam": bam,
    "bai": bai
    }

def markdown_report_aggregate_profiles_input(w): # I need to add here the fastp files (now I am looking for them within the markdown file)
    if pd.notna(samples_table['peak_ctrl_file_alias']).any():
        peak_summary = f"Results/{w.genomes_not_fused}/Analysis/Peaks_summary.tsv",
    else:
        peak_summary = []
    
    mm10_agg_profiles = {}
    B6xCAST_agg_profiles = {} # What happens when the input in snakemake is an empty dictionary?

    if ( (samples_table['top5000_HS_heatmap'] &
     samples_table['get_single_strand'] & ~samples_table['B6xCAST']).any() ):
        # 'Lib_name_aggregate_profiles_data.RData' is the result of processing all top5000, 
        # asymmetrix and XnonPAR data (output of 'wrangle_X_nonPAR_asymmetric_HS' rule.
        mm10_agg_profiles={
            "mm10_top5000_asmtric_auto_XnonPAR_ag_profs": f"Results/{w.genomes_not_fused}"
            "/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/"
            f"{config['library']['name']}_aggregate_profiles_data.RData"
                            }
    elif ( (samples_table['top5000_HS_heatmap'] & ~samples_table['B6xCAST']).any() ):
        mm10_agg_profiles={
            "mm10_top_5000_ag_profs": f"Results/{w.genomes_not_fused}/Analysis"
            "/Heatmaps_and_aggregate_profiles/Hotspots/"
            f"{config['library']['name']}_top_5000_plus_minus_2000.RData"
                            }
    
    if ( (samples_table['top5000_HS_heatmap'] & samples_table['B6xCAST']).any() ):
            B6xCAST_agg_profiles={
                "B6xCAST_top_5000_pm_1000bp": f"Results/{w.genomes_not_fused}/Analysis"
                "/Heatmaps_and_aggregate_profiles/Hotspots/"
                f"{config['library']['name']}_B6xCAST_top_5000_pm_1000bp.RData",
                "B6xCAST_PRDM9_assymetric_hs_invading_strand": "Results/"
                "mm10_x_CAST_EiJ/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/"
                f"{config['library']['name']}_B6xCAST_PRDM9_assymetric_hs_invading_strand.RData",
                "B6xCAST_PRDM9_assymetric_hs_receiving_strand": "Results/"
                "mm10_x_CAST_EiJ/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/"
                f"{config['library']['name']}_B6xCAST_PRDM9_assymetric_hs_receiving_strand.RData"
                                }
        # input = [f"Results/{w.genomes_not_fused}/Analysis/Peaks_summary.tsv"] + \
        # samples_table_2['processed_flagstat'].values.tolist()
        
    return {
        "peaks_summary": peak_summary,
        "bamfiles_reads": samples_table_2['processed_flagstat'].values.tolist(),
        } | mm10_agg_profiles | B6xCAST_agg_profiles # The '|' is to concat. dictionaries

def merge_bams_input(w):
    # Information on the merge is on the 'merge_with' column from samples table,
    # and it should be a space-separated list of: the name of another sample on 
    # the samples table, or a path to a file.
    samples_orig = samples_table.loc[w.sample, 'merge_with'].split()
    samples = [f"Results/{w.sample}.{w.genomes_all}.bam"]
    # If sample to merge with is a path, use it as is, otherwise look for the
    # raw bam file on samples_table_2
    for sample in samples_orig:
        if re.search("/", sample):
            samples+=[sample]
        else:
            second_bam=f"Results/{sample}.{w.genomes_all}.bam"
            samples+=[second_bam]
    return samples

def dros_normalization_report_input(w):
    df=samples_table_2.loc[samples_table_2['dros_equalization_group'].str.match(w.dros_eq_group, na=False), :]
    return df['processed_flagstat_dros'].tolist()

def dros_normalization_input(w):
    if w.strand == "":
        bam = samples_table_2.loc[w.sample, 'dedup_flt_both_strds_bam']
    elif ("83-163" in w.strand) | ("inc_16" in w.strand):
        bam = samples_table_2.loc[w.sample, 'ss_83_or_i16_bam']
    elif ("99-147" in w.strand) | ("exc_16" in w.strand):
        bam = samples_table_2.loc[w.sample, 'ss_99_or_e16_bam']
    
    bai = bam + ".bai"
    return {
    "report" : f"Results/d6/Analysis/drosophila_normalization/{w.dros_eq_group}/drosophila_equalization_report.tsv",
    "bam" : bam,
    "bai" : bai
    }

def intersect_peaks_HSs_list_input(w):
    if samples_table.loc[w.sample, 'B6xCAST']:
        hotspots = config['references']['mm10']['dmc1']['B6xCAST_top_5000_pm_1000bp']
    else:
        hotspots = config['references']['mm10']['spo11']['top_5000_plus_minus_2000']
    return{
        "hotspots" : hotspots,
        "sample_peaks" : f"Results/{w.genomes_not_fused}/Peaks/MACS2/{w.peak_type}/{w.peak_params}/Black-grey_filtered/{w.sample}.{w.genomes_final}{w.extension}.{w.peak_type}Peak"
        }

def compute_matrix_outfiles_hs_input(w):
    hotspot_protein="spo11"
    if samples_table.loc[w.sample, 'B6xCAST']:
        hotspot_protein="dmc1"
    
    hotspots = config['references']['mm10'][hotspot_protein][w.hs_region]
    bigwig = f"Results/{w.genomes_not_fused}/Bigwigs/Coverage/{w.strands}/{w.cov_params}/{w.sample}.{w.genomes_final}{w.extension}{w.strand}.bw"
    
    return {"region": hotspots,
            "bigwig": bigwig}