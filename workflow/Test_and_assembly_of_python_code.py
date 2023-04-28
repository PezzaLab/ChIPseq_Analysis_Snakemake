#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 22 20:22:21 2022

@author: quio
"""
# Set current wd (I cannot get Spyder IDE to set wd as file's path...)
os.chdir("/Volumes/Pezza/hpc-nobackup/Agustin/test_folder/ChIPseq_Analysis_Snakemake/workflow")

import pandas as pd
import numpy as np
import os
import re
import yaml
import functools as ft
import smk_functions as smkf
import suffixes as sfxs

# from snakemake.utils import Paramspace

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

#%% Tests

