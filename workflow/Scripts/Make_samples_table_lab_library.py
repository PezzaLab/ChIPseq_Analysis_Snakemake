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
import os
import re
import shutil as shu

# Get input by user (library name, sequencing techonology, PE/SR,)


def get_lib_name():
    library_name = input(
        "What is the library name?\nKeep in mind that it has to be"
        "the same name provided to Stuart Glenn\n")
    return {"archive": f"/archive/pezza/Agustin/{library_name}",
            "scratch": f"/s/pezzar-lab/{library_name}"}


while True:
    library_paths = get_lib_name()
    if not os.path.exists(library_paths['archive']):
        print(f"The directory '{library_paths['archive']}' does not exist. "
              "Please try again")
        continue
    elif not os.path.isdir(library_paths['archive']):
        print(f"The path '{library_paths['archive']}' is not a directory."
              "Please try again")
        continue
    else:
        break

# Get library sequencing mode (PE or SR)


def get_seq_mode():
    answer = int(
        input("Are reads paired (PE) or single (SR)?, "
              "choose a number.\n1 = PE\n2 = SR\n\n"
              "If different samples have different sequencing pairing of reads"
              " you can always change this on Config/samples_table.csv\n"
              )
    )
    return answer


while True:
    pe_sr = get_seq_mode()
    if pe_sr != 1 and pe_sr != 2:
        print(f"You need to input either '1' or '2', you can do it! ")
        continue
    else:
        break

pe_dic = {1: True, 2: False}

# Get library technology (regular or adaptase)


def get_library_tech():
    answer = int(input(
        "\nWhat technology was used to do the library? "
        "Choose a number\n1 = Adaptase\n2 = Regular\n3 = Other\n")
    )
    return answer


while True:
    lib_tech = get_library_tech()
    if lib_tech < 1 or lib_tech > 3:
        print(f"You need to input either '1', '2' or '3', you can do it! ")
        continue
    else:
        break

lib_tech_dic = {1: "adaptase", 2: "regular", 3: "'Other"}

# genome=input("What is the reference genome for most of your samples? \
#              Choose a number\n1 = mm10\n2 = Other (fill up manually)\n")
# genome=int(genome)
# if ((genome < 1) | (genome > 2)):
#     print("You need to write a number from 1 to 2. Exiting.")
#     exit()
# genome_dic={1: "mm10", 2:"'mm10' or 'hg19' or 'hg38'"}

# Create library folder and copy snakemake pipeline
source_path = (
    "/Volumes/Pezza/hpc-nobackup/Agustin/"
    "test_folder/ChIPseq_Analysis_Snakemake"
)
dest_path = library_paths['scratch']

not_copy = shu.ignore_patterns(
    '.*', 'tmp*', '_*_', 'Test_and_assembly_of_python_code.py',
    'samples_table_processed.csv', 'Results*', 'logs*', '*dry_run*',
    'commands_develop.sh', 'dag*', 'Test_code*'
)
try:
    shu.copytree(source_path, dest_path, ignore=not_copy)
except FileExistsError:
    overwrite = int(
        input(f"'{dest_path}' already exists, would you like to overwrite? \n"
              "1) Yes\n"
              "2) No\n"
              )
    )
    if overwrite:
        shu.copytree(
            source_path, dest_path, dirs_exist_ok=True,
            ignore=not_copy
        )
    else:
        print("Quiting now")
        exit
else:
    print(f"Copying {dest_path}")

# Get fastqs' filepaths
fastqs_temp = os.listdir(library_paths['archive'])
fastqs_temp2 = [x for x in fastqs_temp if re.search(".*\.fastq\.gz$", x)]
fastqs = sorted(fastqs_temp2, key=str.lower)
fastqs = [library_paths['archive'] + "/" +
          file for file in fastqs]

# Define values of table
if pe_dic[pe_sr]:
    fastq1_paths = fastqs[0::2]
    fastq2_paths = fastqs[1::2]
else:
    fastq1_paths = fastqs[0:]
    fastq2_paths = "-"

exp_names = []
for a in fastq1_paths:
    b = os.path.basename(a).strip("\n")
    exp_names += [re.sub("_S[0-9]+_R[0-9]_[0-9]+.fastq.gz", '', b)]

# Assemble dataframe
sample_table = pd.DataFrame(
    {"sample_name": exp_names,
     "fastq1": fastq1_paths,
     "fastq2": fastq2_paths,
     "PE": pe_dic[pe_sr],
     "library_technology": lib_tech_dic[lib_tech],
     "reference_genome": "mm10",
     "peak_ctrl_file_alias": ("'chip_mm10_sonic_testes_adaptase_1' or "
                              "'chip_mm10_sonic_testes_regular_1' or "
                              "'chip_mm10_MNAse_testes_adaptase_7dpp' or "
                              "'chip_mm10_sonic_testes_accel_7dpp'or "
                              "'cyr_mm10_testes_adaptase_1' or "
                              "'-' for no peak calling"),
     "dros_spike_in": "'True' or 'False'",
     "get_single_strand": "'True' or 'False'",
     "Clip_reads_to_1bp_on_5_prime": "'True' or 'False'",
     "top5000_HS_heatmap": "'True' or 'False'",
     "Size_DNA_top_5000_HS": "'True' or 'False'",
     "merge_with": "-",
     "dros_equalization_group": "-",
     "B6xCAST": "False",
     })

# Save table
sample_table.to_csv(f"{library_paths['scratch']}/Config/samples.csv",
                    index=False,
                    na_rep="-",
                    )
# sample_table.to_csv("/Users/quio/Test_folder_local/samples.csv",
#                     index=False,
#                     na_rep="-",
#                     )

# Final message
print(
    "\n\n"
    "Done! The snakemake pipeline has been copied at the following path: "
    f"{library_paths['scratch']}\n"
    "Please review the following pipeline files within the 'Config' folder:\n"
    "1) samples_table.csv.\n"
    "2) config.yaml\n"
    "3) Profiles/slurm_quio/config.yaml\n\n"
    "For more information on how to fill the files please visit "
    "https://github.com/PezzaLab/ChIPseq_Analysis_Snakemake\n"
    "Thanks and have a nice day!...biaatch!"
)

# Run this script:
# python /Volumes/Pezza/hpc-nobackup/Agustin/test_folder/ChIPseq_Analysis_Snakemake/workflow/Scripts/Make_samples_table_lab_library.py
