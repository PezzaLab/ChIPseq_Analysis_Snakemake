#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 14:59:00 2024

@author: quio

Usage:
    cd /project/path && python workflow/Scripts/rename_samples.py

Usage example in the cluster:
    ml slurm python/3.14.7 && \
    cd /project/path && \
    python workflow/Scripts/rename_samples.py
    
Requirements:
    ./Config/renaming_info.tsv: 
        1st column should contain the old sample name (full), 2nd column the 
        new sample name
    python/3.14.7:
        This script should be run with python/3.14.7 loaded (which includes pandas and numpy)

Description:
    Script will rename files and folders within the '../../Results/'
    directory using information stored at "./Config/renaming_info.tsv".
    It will also modify the 'samples.csv' file with the new names and save a
    copy of the old samples.csv file named as 'samples.original_i.csv', where
    i is the iteration index.
    How do the script changes the file names?
    It will do so by looking for a specific STRING (first column of
    renaming_info.tsv) and replacing it by another STRING (second column in
    renaming_info.tsv) in all files and subdirectories' names within ./Results
    How to fill "./Config/renaming_info.tsv"?
    First column is the word to be replaced and second the replacement word.
    What if there there is already a 'renaming_info.tsv' file because this is
    the second renaming event?
    In this case, you should save the original file with a different name and
    modify the existinf file with the new names.
    samples.csv and save it as 'samples.original_i.csv', being i the iteration
    number.

Output:
    './Config/samples.original_i.csv': 
        copy of the previous 'sample.csv' file, where i is a number that grows
        with every renaming iteration.
    
Tips:
    To make sure you don't rename something unintentionally, always use full
    experiment name in first column of "./Config/renaming_info.tsv".
"""

import os
import shutil as shu
import pandas as pd


def rename_files(folder_path, old_word, new_word):
    """
    This function recursively walks through a directory and renames all files
    and directories within according to a name-conversion table.

    Args:
    -----
        folder_path <str>:
            The path to the directory to search.
        old_word <str>:
            The word to be replaced in the filenames.
        new_word <str>:
            The word to replace the old word with.
    """
    old_files = []
    new_files = []
    old_dirs = []
    new_dirs = []

    # Change names
    for root, dirnames, filenames in os.walk(folder_path, topdown=False):
        for filename in filenames:
            if old_word in filename:
                new_filename = filename.replace(old_word, new_word)
                old_filepath = os.path.join(root, filename)
                new_filepath = os.path.join(root, new_filename)
                os.rename(old_filepath, new_filepath)
                old_files.append(old_filepath)
                new_files.append(new_filepath)
        for dirname in dirnames:
            if old_word in dirname:
                old_dirpath = os.path.join(root, dirname)
                new_dirname = dirname.replace(old_word, new_word)
                new_dirpath = os.path.join(root, new_dirname)
                os.rename(old_dirpath, new_dirpath)
                old_dirs.append(old_dirpath)
                new_dirs.append(new_dirpath)
    if (len(old_dirs) > 0) & (len(new_dirs) > 0):
        print("Directories:\n")
        for old, new in zip(old_dirs, new_dirs):
            print(f"Renamed dir: {old} to {new}\n")
        print("\n\nFiles:\n")
        for old, new in zip(old_files, new_files):
            print(f"Renamed file: {old} to {new}")
    elif len(old_dirs) > 0:
        print("No files were found for this name.\n")
        print("Directories:\n")
        for old, new in zip(old_dirs, new_dirs):
            print(f"Renamed dir: {old} to {new}\n")
    else:
        print("No dirs or files were found for this name.\n")


#%% Rename files
# rename_info_tsv_path = "/Volumes/pezza/Agustin/test_folder/ChIPseq_Analysis_Snakemake/Config/renaming_info.tsv"
rename_info_tsv_path = "./Config/renaming_info.tsv"
renaming_info = pd.read_csv(
    rename_info_tsv_path,
    sep="\t",
    names=["old_name", "new_name"],
    index_col=False,
    ).set_index("old_name", drop=False)

# folder_path = "/Volumes/pezza/Agustin/test_folder/ChIPseq_Analysis_Snakemake"
folder_path = "./Results"
for i in renaming_info.index:
    old_name, new_name = renaming_info.loc[i, ["old_name", "new_name"]]
    print("----------------------------------")
    print(old_name)
    rename_files(folder_path, old_name, new_name)
    print("\n\n")

#%% Make copy of original samples.csv
# In case more than one consecutive renaming events occur, the original 
# samples.csv will be saved as samples.original_i.csv, i being the interation
# index.
# In the first renaming event, the script will make a copy of samples table
# called 'samples.original_1.csv', and modify the 'samples.csv' with new names.

samples_csv_dest = os.path.join("Config", "samples.original_1.csv")
samples_csv_source = os.path.join("Config", "samples.csv")

i = 1
while (os.path.exists(samples_csv_dest)) & (i < 20):  # i<20 just in case
    i += 1
    samples_csv_dest = os.path.join("Config", f"samples.original_{i}.csv")
    if i > 1:
        samples_csv_source = os.path.join("Config",
                                          f"samples.original_{i-1}.csv")

print(f"Making copy of {samples_csv_source} in file {samples_csv_dest}")
shu.copy(samples_csv_source, samples_csv_dest)

#%% Modify samples.csv['sample_name']

# Change column 'sample_name'
samples_table = pd.read_csv(
    "Config/samples.csv",
    true_values=["True", "TRUE", "T"],
    false_values=["False", "FALSE", "F"],
    na_values={"merge_with": "-",
               "fastq2": "-",
               "peak_ctrl_file_alias": "-",
               },
    comment='#',
    ).set_index("sample_name", drop=False)

for i in renaming_info.index:
    old_name, new_name = renaming_info.loc[i, ["old_name", "new_name"]]
    samples_table['sample_name'] = samples_table['sample_name'].str.replace(
        old_name,
        new_name)
# Save modified samples.csv
samples_table.to_csv(path_or_buf="./Config/samples.csv", index=False)

print("Finished!")
