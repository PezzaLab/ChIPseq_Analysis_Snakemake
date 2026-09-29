#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 20 19:59:14 2022

@author: quio
"""
########################################################################
# USAGE:
#   Create a directory on scratch with the library name (or user/library),
#   and another directory within called 'FASTQs', where all .fastq.gz
#   files should be found. e.g.: /s/pezzar-lab/agustin/lib2/FASTQs
#   Run this script directly from IDE or at the command line by running:
#       ml python/3.10.2 pandas/1.4.2 && python <this/script/path>
########################################################################

# %% Imports
import os
import re
import shutil as shu
import pandas as pd
import yaml
import git  # To get snakepipeline current's commit hash
from datetime import datetime

# %% Functions and default values

base_path = "/s/pezzar-lab"
pipeline_path = "/hpc-prj/pezza/Agustin/test_folder/ChIPseq_Analysis_Snakemake"


def get_lib_name():
    raw_input = input(
        f"What is the library name?\nLibrary should be on {base_path}/ "
        "it should contain a directory called 'FASTQs' with all the fastq.gz "
        f"files in it. e.g.: {base_path}/<library_name>/FASTQs or "
        f"{base_path}/<user>/<library_name>/FASTQs\n"
    ).strip()

    norm_path = os.path.normpath(raw_input)
    if os.path.isabs(norm_path):
        if norm_path.startswith(base_path):
            rel_path = os.path.relpath(norm_path, base_path)
            scratch_path = norm_path
        elif norm_path.startswith("/s/pezza"):
            rel_path = os.path.relpath(norm_path, "/s/pezza")
            scratch_path = os.path.join(base_path, rel_path)
        else:
            rel_path = norm_path.lstrip("/")
            scratch_path = norm_path
    else:
        rel_path = norm_path.lstrip("/")
        scratch_path = os.path.join(base_path, rel_path)

    lib_name = os.path.basename(scratch_path)
    fastqs_path = os.path.join(scratch_path, "FASTQs")

    return {
        "library_name": lib_name,
        "rel_path": rel_path,
        "fastqs": fastqs_path,
        "scratch": scratch_path,
    }


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
# %%% Get library name
while True:
    library_paths = get_lib_name()
    if not os.path.exists(library_paths['fastqs']):
        print(f"The directory '{library_paths['fastqs']}' does not exist. "
              "Please try again")
        continue
    elif not os.path.isdir(library_paths['fastqs']):
        print(f"The path '{library_paths['fastqs']}' is not a directory."
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

# %% Create library folder and copy snakemake pipeline (not samples.csv)
dest_path = library_paths['scratch']

# Files and dirs not to copy (not copy samples.csv, which might be in
# destination already):
not_copy = shu.ignore_patterns(
    '.*', 'tmp*', '_*_', 'Test_and_assembly_of_python_code.py',
    'samples_table_processed.csv', 'Results*', 'logs*', '*dry_run*',
    '*Dry_run*', 'commands*', 'dag*', 'Test_code*', 'test_FASTQs*',
    'rstudio-server*', 'slurm-*', 'benchmarks*', 'samples.csv', 'log*',
    'renaming_info.tsv', 'rulegraph.svg', 'Make_samples_table*',
    'setup_vm*', 'papers*', '*env_info*', '__pycache__', '*.pyc',
    '*.Rhistory', '*.RData', 'standalone_*'
)

shu.copytree(
    pipeline_path, dest_path, dirs_exist_ok=True,
    ignore=not_copy
)

# %% Write or update 'commands_X.sh'

# Get pipeline-git info
repo = git.Repo(
    pipeline_path
)
sha = repo.head.object.hexsha
mod_time = str(repo.head.object.committed_datetime)

# If commands_X file already there just update with new pipeline-git info,
commands_source_path = f"{pipeline_path}/commands.sh"
commands_destination_path = os.path.join(
    library_paths['scratch'], f'commands_{library_paths["library_name"]}.sh'
    )

if os.path.exists(commands_destination_path):
    # Get the current date
    current_date = datetime.now()
    # Format the date as "yyyy-mm-dd"
    formatted_date = current_date.strftime("%Y-%m-%d")

    text_append = (f"\n# On {formatted_date}, all pipeline was updated with "
                   f"commit hash:\n# {sha}\n")
    with open(commands_destination_path, 'a') as file:
        file.write(text_append)
else:
    # Read in the file
    with open(commands_source_path, 'r') as file:
        commands = file.read()
    # Replace the target strings
    commands = commands.replace('/s/pezzar-lab/{library_name}',
                                library_paths['scratch'])
    commands = commands.replace('/s/pezza/{library_name}',
                                library_paths['scratch'])
    commands = commands.replace('dropboxOMRF:Bioinformatics/Libraries/{library_name}',
                                f'dropboxOMRF:Bioinformatics/Libraries/{library_paths["rel_path"]}')
    commands = commands.replace('{library_name}',
                                library_paths["library_name"])
    commands = commands.replace('{sha}', sha)
    commands = commands.replace('{commit_date}', mod_time)
    # Write the file out again

    with open(commands_destination_path, 'w') as file:
        file.write(commands)

# %% Modify config file
config_path = os.path.join(library_paths['scratch'], "Config", "config.yaml")
with open(config_path, 'r') as file:
    config = file.read()
# Replace the target string
config = config.replace('Test_library', library_paths["library_name"])
# Write the file out again
with open(config_path, 'w') as file:
    file.write(config)
# %% Do samples_table

fastqs_temp = os.listdir(library_paths['fastqs'])
fastqs_temp2 = [x for x in fastqs_temp if re.search(r".*\.fastq\.gz$", x)]
fastqs = sorted(fastqs_temp2, key=str.lower)
fastqs = [os.path.join(library_paths['fastqs'], file) for file in fastqs]

# Verification step
if pe_dic[pe_sr] and len(fastqs) % 2 != 0:
    raise ValueError(
        f"Paired-End (PE) mode selected, but an odd number of FASTQ files ({len(fastqs)}) "
        f"was found in '{library_paths['fastqs']}'. Each sample requires both mates."
    )

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
    exp_names += [re.sub(".fastq.gz", '', b)]

dros = ["True" if re.search("_dros|CyR", i) else "False" for i in exp_names]
B6xCAST = ["True" if re.search("_B6XCAST_", i, flags=re.IGNORECASE)
           else "False" for i in exp_names]

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
     "B6xCAST": B6xCAST,
     })

# Save table
try:
    sample_table.to_csv(os.path.join(library_paths['scratch'], "Config", "samples.csv"),
                        index=False,
                        na_rep="-",
                        mode="x"  # Do not overwrite
                        )
except FileExistsError:
    overwrite = input(
        "'samples.csv' already exists, would you like to overwrite? \n"
        "1) Yes\n"
        "2) No\n"
    )
    positive = ["1", "y", "Y", "yes", "YES", "Yes"]
    if overwrite in positive:
        sample_table.to_csv(f"{library_paths['scratch']}/Config/samples.csv",
                            index=False,
                            na_rep="-",
                            mode="w"  # Overwrites file
                            )
        print("samples.csv copied")
    else:
        print("Keeping original samples.csv file")
else:
    print("samples.csv copied")


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
# ml slurm python/3.10.2 pandas/1.4.2 && python /hpc-prj/pezza/Agustin/test_folder/ChIPseq_Analysis_Snakemake/workflow/Scripts/Make_samples_table_internet.py
