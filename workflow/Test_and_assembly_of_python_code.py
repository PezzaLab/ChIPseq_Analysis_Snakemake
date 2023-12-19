#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 20:22:21 2022

@author: quio
"""

# %% Loading and processing
import pandas as pd
import numpy as np
import re
import yaml
import os
# Set current wd (I cannot get Spyder IDE to set wd as file's path...)
os.chdir("/Volumes/Pezza/hpc-nobackup/Agustin/test_folder/"
         "ChIPseq_Analysis_Snakemake/workflow")
import smk_functions as smkf
import suffixes as sfxs

# Read sample table
samples_table = pd.read_csv("../Config/samples.csv",
                            true_values=["True", "TRUE", "T"],
                            false_values=["False", "FALSE", "F"],
                            comment='#',
                            na_values={"dros_equalization_group": "-",
                                       "merge_with": "-",
                                       "peak_ctrl_file_alias": "-",
                                       "fastq2": "-"}).set_index("sample_name",
                                                                 drop=False)

# Add merged samples to samples_table
for sample in samples_table['sample_name']:
    if pd.notnull(samples_table.loc[sample, 'merge_with']):
        new_row = samples_table.loc[sample,]
        new_row['sample_name'] = f"{new_row['sample_name']}_MERGED"
        new_row['merge_with'] = "-"
        new_row['fastq1'] = np.nan
        new_row['fastq2'] = np.nan
        samples_table = pd.concat(
            [samples_table, new_row.to_frame().T], axis=0, join='outer')
        # to concatenate a df with a series, I need to convert series to df, and to get
        # the columns right I need to transpose the tabel  (.T). the axis=0 and join='outer'
        # are not really necessary because those are the default values for concat.
        # I put them just as a remininder and for learning purpos
        samples_table = samples_table.set_index("sample_name", drop=False)

# Read config file
with open("../Config/config.yaml", 'r') as stream:
    try:
        config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)
    finally:
        stream.close()

samples_table_2 = sfxs.generate_samples_table_2(samples_table, config)
# samples_table_2.to_csv("Config/samples_table_processed.csv")

cov_config_params_string = f"{config['coverage']['normalization']}_bs{config['coverage']['bin_size']}_sm{config['coverage']['smooth']}"

############################################
#####   Add variables to smkf module   #####
############################################
# (otherwise those variables are not accesible to that module)
smkf.samples_table = samples_table
smkf.config = config
# Remove this from here when move previous to smk_functions.py
smkf.cov_config_params_string = cov_config_params_string
smkf.samples_table_2 = samples_table_2
# %% Set wildcards
# Create wildcards to test code
class Wildcard():
    def __init__(self):
        self
w = Wildcard()
# Set wilcard attributes
w.sample = ""
w.hs_region = ""
w.genomes_not_fused = "mm10"
w.genomes_final = ""
w.extension = ".q_filt.srt.nodup.mit_filt"
w.peak_type = ""
w.peak_params = ""
w.smooth = "smoothed"
w.hs_region = "B6xCAST_efojf"
# Filter columns by regex
samples_table_2.filter(regex = ".*")
#%% Do tests
# Get input for process_aggregate_profiles_clipped_input:
# all matrix files to be processed, which there is one per each bw file
# Type of cov files to eventually deal with:
    # NAME.mm10.q_filt.srt.nodup.mit_filt.83-163.inc_16.clipped_1_bp.bw
    # NAME.mm10.q_filt.srt.nodup.mit_filt.99-147.exc_16.clipped_1_bp.bw
    # NAME.mm10.q_filt.srt.nodup.mit_filt.exc_16.clipped_1_bp.bw
    # NAME.mm10.q_filt.srt.nodup.mit_filt.inc_16.clipped_1_bp.bw
# Available wildcards in their context
Results/{genomes_not_fused}/Analysis/Heatmaps_and_aggregate_profiles/"
    "Hotspots/{libary}_{hs_region}_clipped.{smooth}.RData"
# Target path
"Results/{genomes_not_fused}/Analysis/"
            "Heatmaps_and_aggregate_profiles/Hotspots/{hs_region}/{strands}/"
                "{cov_params}/matrixes/"
                "{sample}.{genomes_final}.q_filt.srt.nodup.mit_filt."
                "{ss_condit}{clip_condit}matrix"

strands = "Both_strands|Single_strand/(1bp_clipped_reads|Full_length_reads)"
clip_condit = r"(((inc|exc)_16\.)?clipped_1_bp\.)?"

# Function
cov_params = (f"{config['coverage']['normalization']}_"
              f"bs{config['coverage']['bin_size']}_"
              f"sm{config['coverage']['smooth']}_"
              f"ex{config['coverage']['extend_reads']}")

strand_se = ["inc_16", "exc_16"]

strand_pe = [ a + "." + b
             for a in ["83-163", "99-147"]
             for b in strand_se
             ]
                
list({
  "Results/"
      + f"{w.genomes_not_fused}"
      + "/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/"
      + f"{w.hs_region}/Single_strand/1bp_clipped_reads/"
      + f"{cov_params}/matrixes/{sample}."
      + f"{genome}.q_filt.srt.nodup.mit_filt."  
      + (f"{s_pe}" if samples_table_2["PE"][sample]
         else f"{s_se}")
      + ".matrix" 
      for sample, genome in zip(samples_table_2.index, 
                                samples_table_2['final_genome'])
      if samples_table_2['Clip_reads_to_1bp_on_5_prime'][sample]
      for s_pe in strand_pe
      for s_se in strand_se
})

def a():
    return ["a", "b"]
a()
