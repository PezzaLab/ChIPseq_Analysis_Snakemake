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
import shutil as shu
# Set current wd (I cannot get Spyder IDE to set wd as file's path...)
os.chdir("/Volumes/pezza/Agustin/test_folder/ChIPseq_Analysis_Snakemake/workflow")
os.chdir("/Volumes/pezza/Agustin/test_folder/ChIPseq_Analysis_Snakemake")
import smk_functions as smkf
import suffixes as sfxs

# Read sample table
samples_table_no_merged_samples = pd.read_csv("../Config/samples.csv",
                            true_values=["True", "TRUE", "T"],
                            false_values=["False", "FALSE", "F"],
                            comment='#',
                            na_values={"dros_equalization_group": "-",
                                       "merge_with": "-",
                                       "peak_ctrl_file_alias": "-",
                                       "fastq2": "-"}).set_index("sample_name",
                                                                 drop=False)


# Add merged samples to samples_table
new_rows = samples_table_no_merged_samples.loc[
    pd.notna(samples_table_no_merged_samples['merge_with'])
]
# Generate the new names
new_rows_names = new_rows['sample_name'] + "_MERGED"
# Replace name with new name
new_rows= new_rows.assign(sample_name = new_rows_names)
# Replace `merge_with` field with '-'
new_rows= new_rows.assign(merge_with = np.nan)
# new_rows= new_rows.assign(fastq1 = np.nan)
# new_rows= new_rows.assign(fastq2 = np.nan)
# Re-index using new names
new_rows.set_index('sample_name', drop=False, inplace=True,
                    verify_integrity=True)
# Add new rows to samples table
samples_table = pd.concat([samples_table_no_merged_samples, new_rows],
                          verify_integrity=True)

####################################################
##### Add samples names with suffixes to table #####
####################################################
# All potentially generated files with their proper extensions will be
# added to  a dataframe called "samples_table_2".
# Files will be called upon using the samples_table_2 dataframe, which
# will be exported as an csv table on "Config/samples_table_processed.csv"

# Save table for debugging and exploring file names and paths
# samples_table_2.to_csv("Config/samples_table_processed.csv", index=False)

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
w.sample = "f"
w.hs_region = ""
w.genomes_not_fused = "mm10"
w.genomes_final = ""
w.extension = ".q_filt.srt.nodup.mit_filt"
w.peak_type = ""
w.peak_params = ""
w.smooth = "smoothed"
w.hs_region = "B6xCAST_efojf"
# Filter columns by regex
samples_table_2.filter(regex=".*")

#%% Do tests

    