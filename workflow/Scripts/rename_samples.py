#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar  6 14:59:00 2024

@author: quio

Usage:
    cd /project/path && python workflow/Scripts/rename_samples.py

Usage example in the cluster:
    ml python pandas numpy && \
    cd /project/path && \
    python workflow/Scripts/rename_samples.py

Description:
    Script will rename files and folders within the '../../Results/'
    directory using information stored "./Config/renaming_info.tsv".
    It will do so by looking for a specific STRING (first column of
    renaming_info.tsv) and replacing it by another STRING (second column in
    renaming_info.tsv) in all files and subdirectories' names.
    How to fill "./Config/renaming_info.tsv"?
    First column is the word to be replaced and second the replacement word.

Tips:
    To make sure you don't rename something unintentionally, always use full
    experiment name in first column of "./Config/renaming_info.tsv".

Requirements:
    This script should be run on an environment with pandas installed
"""

import os
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
    
    # Set default values
    no_files_with_old_name = True
    no_dirs_with_old_name = True

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


# Read table with renaming info
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
print("Finished!")
