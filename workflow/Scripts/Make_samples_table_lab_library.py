#!/usr/bin/env python3

import os
import re
import shutil
import yaml
import git
import pandas as pd
from datetime import datetime

# ---------------------- CONFIG ----------------------
FASTQ_BASE = "/archive/pezza/Agustin/"
RESULTS_BASE = "/s/pezzar-lab/"
PIPELINE_SRC = "/hpc-prj/pezza/Agustin/test_folder/ChIPseq_Analysis_Snakemake"
IGNORED_FILES = shutil.ignore_patterns(
    '.*', 'tmp*', '_*_', 'Test_and_assembly_of_python_code.py',
    'samples_table_processed.csv', 'Results*', 'logs*', '*dry_run*',
    'commands*', 'dag*', 'Test_code*', 'test_FASTQs', 'rstudio-server*',
    'slurm-*', 'benchmarks*', 'samples.csv', 'log*', 'renaming_info.tsv',
    'rulegraph.svg', 'Make_samples_table*'
)

# ---------------------- INPUT ----------------------
def ask_input(prompt: str, valid: list[str] = None) -> str:
    while True:
        val = input(prompt).strip()
        if not valid or val in valid:
            return val
        print(f"Invalid input. Please enter one of: {valid}")

def get_library_paths() -> dict:
    name = input("Library name (as provided to Stuart Glenn): ").strip()
    return {
        "name": name,
        "archive": os.path.join(FASTQ_BASE, name),
        "scratch": os.path.join(RESULTS_BASE, name),
    }

def choose_subfolder(base_path: str, dirs: list[str]) -> str:
    # Pair each directory name with its creation timestamp
    dir_info = []
    for d in dirs:
        dir_path = os.path.join(base_path, d)
        creation_ts = os.path.getctime(dir_path)
        dir_info.append((d, creation_ts))
    
    # Sort by timestamp descending (newest first)
    dir_info.sort(key=lambda x: x[1], reverse=True)
    
    print("\nMultiple subdirectories found:")
    for i, (d, ts) in enumerate(dir_info, 1):
        creation_date = datetime.fromtimestamp(ts).strftime("%Y/%m/%d")
        print(f"{i}) {d}  (created: {creation_date})")
    
    # Map back to sorted order for selection
    valid_choices = [str(i) for i in range(1, len(dir_info) + 1)]
    index = int(ask_input("Choose a number: ", valid_choices))
    return os.path.join(base_path, dir_info[index - 1][0])

# ---------------------- SETUP ----------------------
def find_fastq_folder(archive: str) -> str:
    contents = os.listdir(archive)
    dirs = [d for d in contents if os.path.isdir(os.path.join(archive, d))]
    files = [f for f in contents if os.path.isfile(os.path.join(archive, f))]

    if not dirs and not files:
        raise RuntimeError(f"No files or folders found in {archive}")
    elif len(dirs) == 1:
        return os.path.join(archive, dirs[0])
    elif len(dirs) > 1:
        return choose_subfolder(archive, dirs)
    else:
        return archive

def copy_pipeline(dest: str):
    if os.path.exists(dest):
        opt = ask_input(f"{dest} exists. Overwrite? (y/n): ", ["y", "n"])
        if opt == "n":
            print("Aborting.")
            exit(1)
    shutil.copytree(PIPELINE_SRC, dest, dirs_exist_ok=True, ignore=IGNORED_FILES)

# ---------------------- METADATA ----------------------
def write_command_script(dest: str, lib_name: str):
    repo = git.Repo(PIPELINE_SRC)
    sha = repo.head.object.hexsha
    date = str(repo.head.object.committed_datetime)

    src_cmd = os.path.join(PIPELINE_SRC, "commands.sh")
    dst_cmd = os.path.join(dest, f"commands_{lib_name}.sh")

    if os.path.exists(dst_cmd):
        note = f"\n# Updated on {datetime.now().date()} with commit {sha}\n"
        with open(dst_cmd, 'a') as f:
            f.write(note)
    else:
        with open(src_cmd) as f:
            content = f.read().replace('{library_name}', lib_name)
            content = content.replace('{sha}', sha).replace('{commit_date}', date)
        with open(dst_cmd, 'w') as f:
            f.write(content)

def update_config_yaml(path: str, lib_name: str):
    with open(path, 'r') as f:
        config = f.read().replace("Test_library", lib_name)
    with open(path, 'w') as f:
        f.write(config)

# ---------------------- SAMPLES TABLE ----------------------
def guess_sample_names(fq1_paths: list[str]) -> list[str]:
    return [
        re.sub("_S[0-9]+_R[0-9]_[0-9]+\.fastq\.gz", '', os.path.basename(f))
        for f in fq1_paths
    ]

def build_sample_table(fastqs: list[str], paired: bool, tech: str, config_path: str) -> pd.DataFrame:
    fastq1 = fastqs[::2] if paired else fastqs
    fastq2 = fastqs[1::2] if paired else ['-'] * len(fastq1)
    names = guess_sample_names(fastq1)

    with open(config_path) as f:
        config = yaml.safe_load(f)
    ctrl_keys = list(config['MACS2']['control'].keys())
    ctrl_alias = " or ".join(ctrl_keys)

    table = pd.DataFrame({
        "sample_name": names,
        "fastq1": fastq1,
        "fastq2": fastq2,
        "PE": paired,
        "library_technology": tech,
        "reference_genome": "mm10",
        "peak_ctrl_file_alias": ctrl_alias,
        "dros_spike_in": ["True" if re.search("_dros|CyR", n) else "False" for n in names],
        "get_single_strand": "True",
        "Clip_reads_to_1bp_on_5_prime": "False",
        "top5000_HS_heatmap": "True",
        "Size_DNA_top_5000_HS": "False",
        "merge_with": "-",
        "dros_equalization_group": "-",
        "B6xCAST": ["True" if re.search("_B6XCAST_", n, re.IGNORECASE) else "False" for n in names],
    })
    return table

def save_sample_table(df: pd.DataFrame, path: str):
    if os.path.exists(path):
        opt = ask_input("samples.csv exists. Overwrite? (y/n) \n"
                        "If so, will create "
                        "backup of existing samples.csv : ", ["y", "n"])
        if opt == "y":
            shutil.copy(path, path.replace(".csv", "_backup.csv"))
            df.to_csv(path, index=False, na_rep="-")
    else:
        df.to_csv(path, index=False, na_rep="-")

# ---------------------- MAIN ----------------------
def main():
    paths = get_library_paths()
    if not os.path.isdir(paths['archive']):
        raise RuntimeError(f"Invalid archive path: {paths['archive']}")

    fastq_dir = find_fastq_folder(paths['archive'])
    fastq_files = sorted([
        os.path.join(fastq_dir, f) for f in os.listdir(fastq_dir)
        if re.search(r"R[12]_\d+\.fastq\.gz$", f)
    ])

    pe = ask_input("1 = PE\n2 = SR\nPaired-end?: ", ["1", "2"]) == "1"
    tech = {"1": "adaptase", "2": "regular"}[ask_input("1 = Adaptase\n2 = Regular\nLibrary tech?: ", ["1", "2"])]

    copy_pipeline(paths['scratch'])
    update_config_yaml(os.path.join(paths['scratch'], "Config", "config.yaml"), paths['name'])
    write_command_script(paths['scratch'], paths['name'])

    df = build_sample_table(fastq_files, pe, tech, os.path.join(paths['scratch'], "Config", "config.yaml"))
    save_sample_table(df, os.path.join(paths['scratch'], "Config", "samples.csv"))

    print(f"\nSetup complete at {paths['scratch']}\nReview Config/samples.csv and Config/config.yaml")

if __name__ == "__main__":
    main()
