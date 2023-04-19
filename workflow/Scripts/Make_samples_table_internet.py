#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 19:59:14 2022

@author: quio
"""
########################################################################
# USAGE:
#   First log in to o3:
#       o3-login -t LDAP_o3-pezza
#   Then run this script directly from IDE or at the command line by running:
#       ml python/3.10.2 pandas/1.4.2 && python <this/script/path>
########################################################################

import pandas as pd
import numpy as np
from pathlib import Path
import os,re


# Variables to be set by user
pe_sr=input("What sequencing technology was used in this library, choose a number.\n1 = PE\n2 = SR\n\nnote: all files have to be sequenced using the same technology\n")
pe_sr=int(pe_sr)
if ((pe_sr < 1) | (pe_sr > 2)):
    print("You need to write a number from 1 to 2. Exiting.")
    exit()
pe_dic={1:True, 2:False}

lib_tech=input("What technology was used to do the library? Choose a number\n1 = Adaptase\n2 = Regular\n3 = Other\n")
lib_tech=int(lib_tech)
if ((lib_tech < 1) | (lib_tech > 3)):
    print("You need to write a number from 1 to 3. Exiting.")
    exit()
lib_tech_dic={1: "adaptase", 2:"regular", 3: "'Adaptase' or 'regular'"}

genome=input("What is the reference genome for your samples? Choose a number\n1 = mm10\n2 = Other (fill up manually)\n")
genome=int(genome)
if ((genome < 1) | (genome > 2)):
    print("You need to write a number from 1 to 2. Exiting.")
    exit()
genome_dic={1: "mm10", 2:"'mm10' or 'hg19' or 'hg38'"}

# Create library folder and subfolders
parent_dir = os.getcwd()
paths = []
for dire in ["Results","Results/FASTQs", "Resources", "Config"]: paths += [parent_dir + "/" + dire]
for dire in paths:
    Path.mkdir(Path(dire), parents=True, exist_ok=True)

# Gets fastqs paths
dir_name = f'{parent_dir}/Results/FASTQs'
if os.path.isdir(dir_name):
    files = os. scandir(dir_name)
else:
    print(f"{dir_name} doesn't exist")
    exit()

fastqs_temp=os.listdir(dir_name)
fastqs_temp2=[x for x in fastqs_temp if re.search(".*\.fastq\.gz$", x)] # Eliminate files no ending in '.fastq.gz.'
fastqs=sorted(fastqs_temp2, key=str.lower) # Sort files in alphabetical order
fastqs=["Results/FASTQs/" + file for file in fastqs] # Add folder prefix

# Define values of table
if pe_dic:
    fastq1_paths = fastqs[0::2]
    fastq2_paths = fastqs[1::2]
else:
    fastq1_paths=fastqs[0:]
    fastq2_paths="-"

exp_names = []
for a in fastq1_paths: 
    b = os.path.basename(a).strip("\n")
    exp_names += [re.sub("(_S[0-9]+)*(_R[0-9])*(_[0-9]+)*.fastq.gz", '', b)]

# Assemble dataframe
sample_table = pd.DataFrame(
    {"sample_name" : exp_names,
    "fastq1" : fastq1_paths,
    "fastq2" : fastq2_paths,
    "PE" : pe_dic[pe_sr],
    "library_technology" : lib_tech_dic[lib_tech],
    "reference_genome" : genome_dic[genome],
    "peak_ctrl_file_alias" : "'chip_mm10_sonic_testes_adaptase_1' or 'chip_mm10_sonic_testes_regular_1' or 'chip_mm10_MNAse_testes_adaptase_7dpp' or 'chip_mm10_sonic_testes_accel_7dpp' or 'cyr_mm10_testes_adaptase_1' or 'no_input' or '-' for no peak calling",
    "dros_spike_in" : "'True' or 'False'",
    "get_single_strand" : "'True' or 'False'",
    "Clip_reads_to_1bp_on_5_prime" : "'True' or 'False'",
    "top5000_HS_heatmap" : "'True' or 'False'",
    "Size_DNA_top_5000_HS" : "'True' or 'False'",
    "merge_with" : "-",
    "dros_equalization_group" : "-"
    })

# Save table
sample_table.to_csv(f"{parent_dir}/Config/samples.csv", index=False, na_rep="-")
# sample_table.to_csv("/Users/quio/Test_folder_local/samples.csv", index=False, na_rep="-")  

# Final message
print("Done! Please review your table ('samples_table.csv') using excel or any other \
      program of the kind.\n\
      Many other fields need to be filled manually. For this, choose one of the options \
      writen in each cell.\n\
      For input files, put '-' (whithout the quotes) on the 'peak_ctrl_file_alias' \
      column so that it won't call for peaks.\n\
      If you want to merge to files you have 2 options:\n\
      if both files are on the same library, just put the sample name of the 2nd \
      file to merge on the coluum 'merge_with' of the first file to merge.\n\
      if both files are on different libraries, put the FULL PATH of the file from \
      the other library on the 'merge_with' column.\n\
      Note that if experiment has drosophila spike-in, you need to use the raw bam \
      file, not the deduplicated/filtered one.\n\
      To use a file within the library to be processed as input for peak calling \
      of the other files within that libary, then you need to put the path of the \
      deduplictated and filtered input file in the column 'peak_ctrl_file_alias'. \
      The path should be input as it will appear on snakemake, eg: 'Results/mm10/Bams/Both_strands/My_Input_zyg.mm10.q_filt.srt.nodup.mit_filt.bam'\n\
      Thanks and have a nice day!...biaatch!")

