# Rule to generate checksums.tsv manifest for Results/ deliverables exported to Dropbox via rclone
# Conforms to specifications from dialog 32cfd667-47a8-4121-88cf-f7cc23f687fd

def get_pipeline_all_outputs(wildcards=None):
    """
    Returns all expected output files depending on pipeline mode (basic vs full).
    Ensures generate_checksums only runs once all upstream deliverables have completed.
    """
    outputs = [
        *multiqc,
        *basic_coverage_bigwigs,
        *basic_peaks,
        *peaks_summary_files,
    ]
    if not is_basic_mode:
        outputs.extend(single_strand_coverage_bigwigs)
        outputs.extend(heatmaps)
        outputs.extend(html_reports)
    return outputs

rule generate_checksums:
    input:
        get_pipeline_all_outputs,
    output:
        "Results/checksums.tsv",
    log:
        "logs/generate_checksums/generate_checksums.log",
    benchmark:
        "benchmarks/generate_checksums/generate_checksums.tsv",
    threads: 2
    envmodules:
        config.get("rclone", "rclone"),
    shell:
        """
        bash workflow/Scripts/generate_rclone_checksums.sh . {output} > {log} 2>&1
        """

# Convenience target alias
rule checksums:
    input:
        "Results/checksums.tsv",
