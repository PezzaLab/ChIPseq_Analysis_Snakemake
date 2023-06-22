import pandas as pd
import numpy as np
import os
import re

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
    
    # Check that none of the sample names is repeated
    if samples_table['sample_name'].duplicated().any():
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
def get_input_align(w):
    if samples_table.loc[w.sample, "PE"]:
        return {
        "fastq1":samples_table.loc[w.sample, "fastq1"],
        "fastq2":samples_table.loc[w.sample, "fastq2"]
        }
    else:
        return {
        "fastq1":samples_table.loc[w.sample, "fastq1"],
        "fastq2":[]} # Cannot use "" here because snakemake will look for "" file
        # and give a "missing input" error. Instead, I give it an empty list


# Function for input/ctrl for peak calling rule
def MACS2_peak_calling_files(w):
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

def get_hs_agg_profs_inputs(w):
    samples = {}
    samples['both_strands'] = samples_table_2.loc[samples_table['top5000_HS_heatmap'],
                                'dedup_flt_both_strds_bw'].values.tolist()
    samples['both_strands_with_ss'] = samples_table_2.loc[samples_table['top5000_HS_heatmap'] &
                                samples_table['get_single_strand'],
                            'dedup_flt_both_strds_bw'].values.tolist() 
        # Here we get the "both strand" big wig, but only from samples that have
        # a single strand profile as well. This is because for all regions other
        # than top 5000 HS, we only do them when they have a single strand profile as well
    samples['83-163_or_i16'] = samples_table_2.loc[samples_table['top5000_HS_heatmap'] &
                                   samples_table['get_single_strand'],
                                   'ss_83_or_i16_bw'].values.tolist()
    samples['99-147_or_e16'] = samples_table_2.loc[samples_table['top5000_HS_heatmap'] &
                                   samples_table['get_single_strand'],
                                   'ss_99_or_e16_bw'].values.tolist()
    # Get list names
    all_files = {} # It is adviced to use dicts to make dinamically names "lists" 
    # (https://stackoverflow.com/questions/14819849/create-lists-of-unique-names-in-a-for-loop-in-python)
    
    # Get sample name from the lists above, assemble the path for the matrix file and add it
    # To the dictionary 'all_files', on the entry named as {strnd}_{hs_region}.  
    # There is 2 both strand lists because depending on the hs_region, we use one or the other.
    if len(samples['99-147_or_e16']) != 0:
        for strnd in ["83-163_or_i16", "99-147_or_e16"]:
            all_files[f"{strnd}_{w.hs_region}"] = []
            for sample in samples[strnd]:
                sample_full_name=os.path.splitext(os.path.basename(sample))[0] 
                    # basename gets filename with ext, splittext[0] takes away extension
                all_files[f"{strnd}_{w.hs_region}"] += (
                    ["Results/mm10/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/"
                    f"{w.hs_region}/Single_strand/Full_length_reads/{cov_config_params_string}/"
                    f"matrixes/{sample_full_name}.matrix"]
                    )

    all_files[f"both_strands_{w.hs_region}"] = []
    if w.hs_region != "top_5000_plus_minus_2000":
        for sample in samples["both_strands_with_ss"]:
            sample_full_name=os.path.splitext(os.path.basename(sample))[0]
            all_files[f"both_strands_{w.hs_region}"] += (
                ["Results/mm10/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/"
                f"{w.hs_region}/Both_strands/{cov_config_params_string}/"
                f"matrixes/{sample_full_name}.matrix"]
                )
    else:
        for sample in samples["both_strands"]:
            sample_full_name=os.path.splitext(os.path.basename(sample))[0] 
            all_files[f"both_strands_{w.hs_region}"] += (
                ["Results/mm10/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/"
                f"{w.hs_region}/Both_strands/{cov_config_params_string}/"
                f"matrixes/{sample_full_name}.matrix"]
                )
    if (all_files == None) | (all_files == []) :
        Print("'All files' empty")
        exit(1)
        
    return all_files

def all_peaks(w):
    genome_filtered = samples_table_2[
    samples_table_2['reference_genome'] == w.genomes_not_fused
    ]
    peaks = {}
    if pd.notna(samples_table['peak_ctrl_file_alias']).any():
        peaks = {
        "narrow_all": genome_filtered.loc[
            genome_filtered['peak_bl_gr_flt_nrw'].notnull(),
            "peak_bl_gr_flt_nrw"].values.tolist(),
        "broad_all": genome_filtered.loc[
            genome_filtered['peak_bl_gr_flt_brd'].notnull(),
            "peak_bl_gr_flt_brd"].values.tolist(),
        "narrow_hs": genome_filtered.loc[
            genome_filtered['peak_bl_gr_flt_hs_int_nrw'].notnull(),
            "peak_bl_gr_flt_hs_int_nrw"].values.tolist(),
        "broad_hs": genome_filtered.loc[
            genome_filtered['peak_bl_gr_flt_hs_int_brd'].notnull(),
            "peak_bl_gr_flt_hs_int_brd"].values.tolist()
            }
    return peaks

def get_qctrl_bams_bais(w):
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

def hmlt_report_input(w): # I need to add here the fastp files (now I am looking for them within the markdown file)
    if pd.notna(samples_table['peak_ctrl_file_alias']).any():
        peak_summary = f"Results/{w.genomes_not_fused}/Analysis/Peaks_summary.tsv",
    else:
        peak_summary = []

    if ( ((samples_table['top5000_HS_heatmap']) &
     (samples_table['get_single_strand'])).any() ):
        return {
        "agg_profiles": f"Results/{w.genomes_not_fused}/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/{config['library']['name']}_aggregate_profiles_data.RData",
        "peaks_summary": peak_summary,
        "bamfiles_reads": samples_table_2['processed_flagstat'].values.tolist()
        }
        # input = [f"Results/{w.genomes_not_fused}/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/top_5000_plus_minus_2000.RData"] + \
        # [f"Results/{w.genomes_not_fused}/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/wrangled_autosomal_x_non_par_and_asymetric.RData"] + \
        # [f"Results/{w.genomes_not_fused}/Analysis/Peaks_summary.tsv"] + \
        # samples_table_2['processed_flagstat'].values.tolist()
    elif (samples_table['top5000_HS_heatmap'].any()):
        return{
        "top_5000_ag_profs": f"Results/{w.genomes_not_fused}/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/top_5000_plus_minus_2000.RData",
        "peaks_summary": peak_summary,
        "bamfiles_reads": samples_table_2['processed_flagstat'].values.tolist()
        }
        # input = [f"Results/{w.genomes_not_fused}/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/top_5000_plus_minus_2000.RData"] + \
        # [f"Results/{w.genomes_not_fused}/Analysis/Peaks_summary.tsv"] + \
        # samples_table_2['processed_flagstat'].values.tolist()
    else:
        return {
        "peaks_summary": peak_summary,
        "bamfiles_reads": samples_table_2['processed_flagstat'].values.tolist()
        }
        # input = [f"Results/{w.genomes_not_fused}/Analysis/Peaks_summary.tsv"] + \
        # samples_table_2['processed_flagstat'].values.tolist()
    # return input

def input_merge_bams(w):
    if re.search("/",samples_table.loc[w.sample,'merge_with']):
        samples=samples_table_2.loc[w.sample, ['raw_bam', 'merge_with']].values.tolist()
    else:
        second_bam_name=samples_table_2.loc[w.sample, 'merge_with']
        second_bam=samples_table_2.loc[second_bam_name, 'raw_bam']
        first_bam=samples_table_2.loc[w.sample, 'raw_bam']
        samples=[second_bam] + [first_bam]
    return samples

def dros_norm_report_input(w):
    df=samples_table_2.loc[samples_table_2['dros_equalization_group'].str.match(w.dros_eq_group, na=False), :]
    return df['processed_flagstat_dros'].tolist()

def dros_norm_input(w):
    if w.strand == "":
        unscaled_bw = samples_table_2.loc[w.sample, 'dedup_flt_both_strds_bw']
        bam = samples_table_2.loc[w.sample, 'dedup_flt_both_strds_bam']
    elif ("83-163" in w.strand) | ("inc_16" in w.strand):
        unscaled_bw = samples_table_2.loc[w.sample, 'ss_83_or_i16_bw']
        bam = samples_table_2.loc[w.sample, 'ss_83_or_i16_bam']
    elif ("99-147" in w.strand) | ("exc_16" in w.strand):
        unscaled_bw = samples_table_2.loc[w.sample, 'ss_99_or_e16_bw']
        bam = samples_table_2.loc[w.sample, 'ss_99_or_e16_bam']

    return {
    "report" : f"Results/d6/Analysis/drosophila_normalization/{w.dros_eq_group}/drosophila_equalization_report.tsv",
    "bam" : bam,
    "unscaled_bw" : unscaled_bw
    }
    # df=samples_table_2.loc[samples_table_2['dros_equalization_group'].str.match(w.dros_eq_group, na=False), :]
    # both_strands_bams = df.loc[: , 'dedup_flt_both_strds_bam'].tolist()
    # both_strands_bais = (df.loc[: , 'dedup_flt_both_strds_bam'] + ".bai").tolist()
    # ss_bams =  df.loc[: , 'ss_83_or_i16_bam'].tolist() + df.loc[:,'ss_99_or_e16_bam'].tolist()
    # ss_bais = (df.loc[:,'ss_83_or_i16_bam'] + ".bai").tolist() + (df.loc[:,'ss_99_or_e16_bam'] + ".bai").tolist()
    #
    # norm_report=[f"Results/{w.genomes_not_fused}/d6/Analysis/{w.dros_eq_group}/drosophila_equalization_report.tsv"]
    # # Files to be returned in any case:
    #
    # default_files = both_strands_bams + both_strands_bais + norm_report
    #
    # if df['get_single_strand'].all():
    #     return default_files + ss_bams + ss_bais
    # else:
    #     return default_files
