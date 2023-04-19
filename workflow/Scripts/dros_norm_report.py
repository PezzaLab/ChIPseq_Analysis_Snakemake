#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec 14 11:09:24 2022

@author: quio
"""
# The output of this script is a .tsv file with 3 columns:
# "sample", "number_dros_reads" and "scaleFactor"


import pandas as pd
import re
import os

samples_table_2 = snakemake.params['samples_table_2']
number_reads=[]
sample_name = []

for file in snakemake.input:
    with open(file) as f:
        for line in f:
            # get number of primary reads
            if re.search("primary$", line): 
                number_reads += [int(re.search("^[0-9]+", line).group(0))]
    sample_name += [samples_table_2.loc[samples_table_2['processed_flagstat_dros'] == file].index[0]]

df=pd.DataFrame(
    {"sample_name" : sample_name,
     "flagstat_file" : snakemake.input,
     "number_dros_reads" : number_reads})

min_reads = min(df['number_dros_reads'])

df["scaleFactor"]=min_reads/df['number_dros_reads']

df.to_csv(snakemake.output[0], sep='\t', index=False)

