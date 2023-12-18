#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 19:59:14 2022

@author: quio
"""
# %% Imports
import os
import re
import shutil as shu
import pandas as pd
import yaml
import git  # To get snakepipeline current's commit hash


# %% Functions

def get_lib_name():
    library_name = input(
        "What is the library name?\nKeep in mind that it has to be "
        "the same name provided to Stuart Glenn\n")
    return {"library_name": library_name,
            "archive": f"/archive/pezza/Agustin/{library_name}",
            "scratch": f"/s/pezzar-lab/{library_name}"}


def get_seq_mode():
    answer = input(
        "Are reads paired (PE) or single (SR)?, "
        "If different samples have different sequencing pairing of reads"
        " you can always change this later on the Config/samples_table.csv\n"
        "choose a number.\n1 = PE\n2 = SR\n"
    )
    return answer


def get_library_tech():
    answer = input(
        "\nWhat technology was used to do the library? "
        "Choose a number\n1 = Adaptase\n2 = Regular\n")
    return answer


# %% Get user's input
# %%% Get input by user (library name, sequencing techonology, PE/SR,)
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

# %%% Get library sequencing mode (PE or SR)
while True:
    pe_sr = get_seq_mode()
    if pe_sr not in ["1", "2"]:
        print("You need to input either '1' or '2', you can do it! ")
        continue
    else:
        break

pe_dic = {"1": True, "2": False}

# %%% Get library technology (regular or adaptase)
while True:
    lib_tech = get_library_tech()
    if lib_tech not in ["1", "2"]:
        print("You need to input either '1' or '2', you can do it! ")
        continue
    else:
        break

lib_tech_dic = {"1": "adaptase", "2": "regular"}

# %% Create library folder and copy snakemake pipeline
source_path = (
    "/Volumes/Pezza/hpc-nobackup/Agustin/"
    "test_folder/ChIPseq_Analysis_Snakemake"
)
dest_path = library_paths['scratch']

not_copy = shu.ignore_patterns(
    '.*', 'tmp*', '_*_', 'Test_and_assembly_of_python_code.py',
    'samples_table_processed.csv', 'Results*', 'logs*', '*dry_run*',
    'commands_develop.sh', 'dag*', 'Test_code*', 'test_FASTQs',
    'rstudio-server*', 'slurm-*'
)
try:
    shu.copytree(source_path, dest_path, ignore=not_copy)
except FileExistsError:
    overwrite = input(
        f"'{dest_path}' already exists, would you like to overwrite? \n"
        "1) Yes\n"
        "2) No\n"
    )
    positive = ["1", "y", "Y", "yes", "YES", "Yes"]
    if overwrite in positive:
        shu.copytree(
            source_path, dest_path, dirs_exist_ok=True,
            ignore=not_copy
        )
    else:
        print("Quiting now")
        exit
else:
    print(f"Copying {dest_path}")

# %% Modify 'commands.sh' file with library-specific info
# %%% Get git info
# Get git curent commit hash
repo = git.Repo(
    "/Volumes/Pezza/hpc-nobackup/Agustin/test_folder/"
    "ChIPseq_Analysis_Snakemake"
)
sha = repo.head.object.hexsha
mod_time = str(repo.head.object.committed_datetime)
# %%% Modify file
# Read in the file
commands_path = f"{library_paths['scratch']}/commands.sh"
with open(commands_path, 'r') as file:
    commands = file.read()
# Replace the target strings
commands = commands.replace('{library_name}', library_paths["library_name"])
commands = commands.replace('{sha}', sha)
commands = commands.replace('{commit_date}', mod_time)
# Write the file out again
commands_new_path = (
    f"{library_paths['scratch']}/"
    f'commands_{library_paths["library_name"]}.sh'
)
with open(commands_new_path, 'w') as file:
    file.write(commands)
# Erase `commands.sh`
os.remove(commands_path)

# %% Modify config file
config_path = f"{library_paths['scratch']}/Config/config.yaml"
with open(config_path, 'r') as file:
    config = file.read()
# Replace the target string
config = config.replace('Test_library', library_paths["library_name"])
# Write the file out again
with open(config_path, 'w') as file:
    file.write(config)
# %% Do samples_table
# Get fastqs' filepaths
fastqs_temp = os.listdir(library_paths['archive'])
fastqs_temp2 = [x for x in fastqs_temp if re.search(r".*\.fastq\.gz$", x)]
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

dros = ["True" if re.search("_dros|CyR", i) else "False" for i in exp_names]
        
# Read config file
with open(config_path, 'r') as stream:
    try:
        config = yaml.safe_load(stream)
    except yaml.YAMLError as exc:
        print(exc)
    finally:
        stream.close()

# Get peak control options from config file
peak_ctrl_file_alias = list(config['MACS2']['control'].keys())
peak_ctrl_file_aliases = " or ".join(peak_ctrl_file_alias)

# Assemble dataframe
sample_table = pd.DataFrame(
    {"sample_name": exp_names,
     "fastq1": fastq1_paths,
     "fastq2": fastq2_paths,
     "PE": pe_dic[pe_sr],
     "library_technology": lib_tech_dic[lib_tech],
     "reference_genome": "mm10",
     "peak_ctrl_file_alias": peak_ctrl_file_aliases,
     "dros_spike_in": dros,
     "get_single_strand": "True",
     "Clip_reads_to_1bp_on_5_prime": "False",
     "top5000_HS_heatmap": "True",
     "Size_DNA_top_5000_HS": "False",
     "merge_with": "-",
     "dros_equalization_group": "-",
     "B6xCAST": "False",
     })

# Save table
sample_table.to_csv(f"{library_paths['scratch']}/Config/samples.csv",
                    index=False,
                    na_rep="-",
                    )

# %% Final message
print(
    "\n\n"
    "Done! The snakemake pipeline has been copied at the following path: "
    f"{library_paths['scratch']}\n"
    "Please review the following pipeline files within the 'Config' folder:\n"
    "1) samples_table.csv.\n"
    "2) config.yaml\n"
    "3) Profiles/slurm_quio/config.yaml\n\n"
    "For more information on how to fill the files please visit "
    "https://github.com/PezzaLab/ChIPseq_Analysis_Snakemake\n\n"
    "Thanks and have a nice day!...biaatch!"
)


# Run this script:
# ml slurm python/3.10.2 pandas/1.4.2 && python /Volumes/Pezza/hpc-nobackup/Agustin/test_folder/ChIPseq_Analysis_Snakemake/workflow/Scripts/Make_samples_table_lab_library.py
