#!/usr/bin/env python3
"""
Export standalone Basic ChIP-seq / CUT&RUN Pipeline.

This script extracts only the files required for standard ChIP-seq / CUT&RUN
processing (up through trimming, alignment, deduplication, filtering, unstranded
BigWigs, MACS2 peaks, blacklist filtering, peak summaries, and MultiQC) into a
self-contained target directory.

Usage:
    python3 workflow/Scripts/export_basic_pipeline.py [output_dir]
Default output_dir:
    ./standalone_basic_pipeline
"""

import sys
import os
import shutil
from pathlib import Path

def export_basic(dest_dir_str="standalone_basic_pipeline"):
    dest_dir = Path(dest_dir_str).resolve()
    repo_root = Path(__file__).resolve().parent.parent.parent

    print(f">>> Exporting standalone Basic pipeline from {repo_root} to {dest_dir}...")
    dest_dir.mkdir(parents=True, exist_ok=True)

    # 1. Rules directory
    rules_src = repo_root / "workflow" / "rules"
    rules_dest = dest_dir / "workflow" / "rules"
    rules_dest.mkdir(parents=True, exist_ok=True)

    basic_rules = [
        "common.smk",
        "trimming.smk",
        "alignment.smk",
        "filtering.smk",
        "qc.smk",
        "coverage.smk",
        "peaks.smk",
    ]
    for r in basic_rules:
        shutil.copy2(rules_src / r, rules_dest / r)
        print(f"  [COPIED] workflow/rules/{r}")

    # 2. Main Snakefile (copy basic.smk as Snakefile)
    shutil.copy2(repo_root / "workflow" / "basic.smk", dest_dir / "workflow" / "Snakefile")
    print("  [COPIED] workflow/basic.smk -> workflow/Snakefile")

    # 3. Python helper modules & scripts
    for mod in ["smk_functions.py", "suffixes.py"]:
        shutil.copy2(repo_root / "workflow" / mod, dest_dir / "workflow" / mod)
        print(f"  [COPIED] workflow/{mod}")

    scripts_dest = dest_dir / "workflow" / "Scripts"
    scripts_dest.mkdir(parents=True, exist_ok=True)
    basic_scripts = [
        "dros_norm_report.py",
        "peaks_summary.py",
    ]
    for s in basic_scripts:
        shutil.copy2(repo_root / "workflow" / "Scripts" / s, scripts_dest / s)
        print(f"  [COPIED] workflow/Scripts/{s}")

    # 4. Config directory
    config_dest = dest_dir / "Config"
    config_dest.mkdir(parents=True, exist_ok=True)
    for c in ["config.yaml", "samples.csv", "multiqc_config.yaml"]:
        src_file = repo_root / "Config" / c
        if src_file.exists():
            shutil.copy2(src_file, config_dest / c)
            print(f"  [COPIED] Config/{c}")

    profiles_src = repo_root / "Config" / "Profiles"
    if profiles_src.exists():
        shutil.copytree(profiles_src, config_dest / "Profiles", dirs_exist_ok=True)
        print("  [COPIED] Config/Profiles/")

    # 5. Bed files: blacklist and chromosomes only
    res_bed_src = repo_root / "Resources" / "bed_files"
    if res_bed_src.exists():
        for genome_dir in res_bed_src.iterdir():
            if genome_dir.is_dir():
                target_gdir = dest_dir / "Resources" / "bed_files" / genome_dir.name
                target_gdir.mkdir(parents=True, exist_ok=True)
                for f in genome_dir.iterdir():
                    if f.is_file() and ("blacklist" in f.name or "chromosomes" in f.name or "excluderanges" in f.name):
                        shutil.copy2(f, target_gdir / f.name)
                        print(f"  [COPIED] Resources/bed_files/{genome_dir.name}/{f.name}")

    # 6. Basic README
    readme_content = """# Basic ChIP-seq / CUT&RUN Pipeline

This is a self-contained, reproducible Snakemake workflow for standard ChIP-seq and CUT&RUN data analysis.

## Features
- FastQ adapter & quality trimming (`fastp`, `cutadapt`)
- Alignment & Sorting (`bwa mem`, `samtools`)
- Duplicate marking & MAPQ filtering (`picard`, `sambamba`, `samtools`)
- Host vs Drosophila spike-in genome read separation
- Library complexity & insert size metrics (`picard`)
- Unstranded BigWig generation (CPM-normalized and Drosophila spike-in normalized)
- Peak calling with controls (`macs2 callpeak` broad and narrow)
- Blacklist filtering (`bedtools intersect`)
- Quality control aggregation (`multiqc`)
- Peak count summary report (`Peaks_summary.tsv`)

## Execution
```bash
# Dry-run
snakemake -n basic

# Run on SLURM cluster
snakemake --profile Config/Profiles/slurm_quio_repeat_3 basic
```
"""
    with open(dest_dir / "README.md", "w") as fp:
        fp.write(readme_content)
    print("  [CREATED] README.md")

    print(f"\n>>> Successfully created self-contained basic pipeline at: {dest_dir}")

if __name__ == "__main__":
    dest = sys.argv[1] if len(sys.argv) > 1 else "standalone_basic_pipeline"
    export_basic(dest)
