#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 26 14:30:50 2022

@author: quio
"""
import os
import pandas as pd
import numpy as np

# pd.options.display.float_format = '{:,.2f}'.format

peak_types = ["narrow_all", "broad_all", "narrow_hs", "broad_hs"]
samples_table_2 = snakemake.params['samples_table_2']

# Get df with samples names as indexes
peak_count = samples_table_2[
    samples_table_2['reference_genome'] == snakemake.wildcards['genomes_not_fused']][[
    'sample_name']].set_index('sample_name')

# Count peaks and add to df        
for peak_type in peak_types:
    for file in snakemake.input[peak_type]:
        # Get file name
        file_name = file.split("/")[-1].split(".")[0]
        # Get number of peaks
        with open(file, 'r') as fp:
            number_peaks = len(fp.readlines())
        # 
        peak_count.loc[file_name, peak_type] = number_peaks

# Calculate % of peaks in HSs
peak_count["% nrw peaks at HS"] = np.where(
    peak_count['narrow_all'] > 0,
    100 * peak_count['narrow_hs'] / peak_count['narrow_all'],0)

peak_count["% brd peaks at HS"] = np.where(
    peak_count['broad_all'] > 0,
    100 * peak_count['broad_hs'] / peak_count['broad_all'],0)

# Format number
peak_count["% nrw peaks at HS"] = peak_count["% nrw peaks at HS"].map('{:,.2f}'.format)
peak_count["% brd peaks at HS"] = peak_count["% brd peaks at HS"].map('{:,.2f}'.format)

peak_count.to_csv(snakemake.output[0], sep="\t")

