# Standalone entrypoint for Basic ChIP-seq / CUT&RUN Analysis
# Self-contained for publication or standard ChIP-seq analyses without meiotic specialization

configfile: "Config/config.yaml"

# Set basic mode
config["mode"] = "basic"

# Common setup, metadata, and wildcard constraints
include: "rules/common.smk"

# Basic processing rules
include: "rules/trimming.smk"
include: "rules/alignment.smk"
include: "rules/filtering.smk"
include: "rules/qc.smk"
include: "rules/coverage.smk"
include: "rules/peaks.smk"
include: "rules/checksums.smk"

ruleorder: dedup_bam > index
ruleorder: merge_bams > align_fastq

#####################################################
#####                TARGET RULE                #####
#####################################################

rule basic:
    input:
        multiqc,
        basic_coverage_bigwigs,
        basic_peaks,
        peaks_summary_files,
        "Results/checksums.tsv",

