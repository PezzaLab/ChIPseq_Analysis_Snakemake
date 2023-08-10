#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 20:22:21 2022

@author: quio
"""

#%% Loading and processing
import os
# Set current wd (I cannot get Spyder IDE to set wd as file's path...)
os.chdir("/Volumes/Pezza/hpc-nobackup/Agustin/test_folder/ChIPseq_Analysis_Snakemake/workflow")
# from snakemake.utils import Paramspace

import pandas as pd
import numpy as np
import re
import yaml
import functools as ft
import smk_functions as smkf
import suffixes as sfxs


## TEST FOLDER sample table
# samples_table = pd.read_csv(
#     "/Volumes/Pezza/hpc-nobackup/Agustin/test_folder/Snake_make/Config/samples.csv",
#     true_values=["True", "TRUE", "T"],
#     false_values=["False", "FALSE", "F"],
#     na_values={"dros_equalization_group": "-",
#                "fastq2": "-"}).set_index("sample_name", drop=False)

## Read sample table
samples_table = pd.read_csv("../Config/samples.csv",
                           true_values=["True", "TRUE", "T"],
                           false_values=["False", "FALSE", "F"],
                           na_values={"dros_equalization_group": "-",
                                      "merge_with":"-",
                                      "peak_ctrl_file_alias":"-",
                                      "fastq2": "-"}).set_index("sample_name", drop=False)

## Add merged samples to samples_table
for sample in samples_table['sample_name']:
    if pd.notnull(samples_table.loc[sample,'merge_with']):
         new_row=samples_table.loc[sample,]
         new_row['sample_name'] = f"{new_row['sample_name']}_MERGED"
         new_row['merge_with']="-"
         new_row['fastq1']= np.nan
         new_row['fastq2']= np.nan
         samples_table=pd.concat([samples_table, new_row.to_frame().T],axis=0, join='outer') 
            # to concatenate a df with a series, I need to convert series to df, and to get
            # the columns right I need to transpose the tabel  (.T). the axis=0 and join='outer'
            # are not really necessary because those are the default values for concat.
            # I put them just as a remininder and for learning purpos
         samples_table=samples_table.set_index("sample_name", drop=False)

## Read config file         
with open("../Config/config.yaml", 'r') as stream:
    try:
        config=yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)
    finally:
        stream.close()
         
samples_table_2 = sfxs.generate_samples_table_2(samples_table, config)
# samples_table_2.to_csv("Config/samples_table_processed.csv")

cov_config_params_string=f"{config['coverage']['normalization']}_bs{config['coverage']['bin_size']}_sm{config['coverage']['smooth']}"

############################################
#####   Add variables to smkf module   #####
############################################
# (otherwise those variables are not accesible to that module)
smkf.samples_table = samples_table
smkf.config = config
smkf.cov_config_params_string = cov_config_params_string # Remove this from here when move previous to smk_functions.py
smkf.samples_table_2 = samples_table_2
#%% Tests
# Create wildcards to test code
class Wildcard():
    def __init__(self):
        self


w = Wildcard()

setattr(w, 'hs_region', "B6xCAST_top_5000_pm_1000bp")

setattr(w, "genomes_not_fused", "mm10")

# %%
  # Peak files
x=samples_table_2.filter(regex=".*(narrow|broad).*annotated")
x.columns
samples_table_2[]