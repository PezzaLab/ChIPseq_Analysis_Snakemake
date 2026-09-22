#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Makes table with information related to number of peaks for all samples

@author: quio

Wildcards
---------
genomes_not_fused: "mm10|mm39|d6|hg19|hg38|mm10_x_CAST_EiJ"

Output
-------
Table (csv file) with number of peaks (broad and narrow) and, if applicable
(only for mm10 and if asked for), number and % of peaks in hotspots.
"""
import pandas as pd
import numpy as np

# Get df with samples names
samples_table_2 = snakemake.params['samples_table_2']
selection_criteria = (
    samples_table_2['reference_genome']
    == snakemake.wildcards['genomes_not_fused']
)
peak_count = samples_table_2.loc[selection_criteria, ['sample_name']].copy()
peak_count.set_index('sample_name', drop=False, inplace=True)

# Count peaks and add to df
peak_types = ["narrow_all", "broad_all", "narrow_hs", "broad_hs"]
for peak_type in peak_types:
    if peak_type in snakemake.input.keys():
        for file in snakemake.input[peak_type]:
            # Get file name
            file_name = file.split("/")[-1].split(".")[0]
            # Get number of peaks
            with open(file, 'r') as fp:
                number_peaks = len(fp.readlines())
            peak_count.loc[file_name, peak_type] = number_peaks

# Calculate % of peaks in HSs (only for mouse if HS peaks were provided)
if (snakemake.wildcards['genomes_not_fused'] in ("mm10", "mm39") and
        "narrow_hs" in peak_count.columns and "broad_hs" in peak_count.columns):
    peak_count["% nrw peaks at HS"] = np.where(
        peak_count['narrow_all'] > 0,
        100 * peak_count['narrow_hs'] / peak_count['narrow_all'], 0)

    peak_count["% brd peaks at HS"] = np.where(
        peak_count['broad_all'] > 0,
        100 * peak_count['broad_hs'] / peak_count['broad_all'], 0)

    # Format number
    peak_count["% nrw peaks at HS"] = peak_count["% nrw peaks at HS"].map('{:,.2f}'.format)
    peak_count["% brd peaks at HS"] = peak_count["% brd peaks at HS"].map('{:,.2f}'.format)

peak_count.to_csv(snakemake.output[0], sep="\t", index=False)
