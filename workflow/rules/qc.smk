# Quality control rules: Picard, samtools stats/flagstat, MultiQC

rule samstats:
    input:
        unpack(smkf.samstats_samtools_flagstat_input),
    output:
        "Results/{genomes_not_fused}/Qctrl/{sample}/{bam_type}/"
        "{sample}.samstats.txt",
    envmodules:
        config["samtools"],
    benchmark:
        "benchmarks/samstats/{genomes_not_fused}_{bam_type}/"
        "{sample}.tsv",
    shell:
        """
        samtools stats {input.bam} > {output}
        """

rule samtools_flagstat:
    input:
        unpack(smkf.samstats_samtools_flagstat_input),
    output:
        "Results/{genomes_not_fused}/Qctrl/{sample}/{bam_type}/"
        "{sample}.flagstat.txt",
    envmodules:
        config["samtools"],
    benchmark:
        "benchmarks/samtools_flagstat/{genomes_not_fused}_{bam_type}/"
        "{sample}.tsv",
    shell:
        """
        samtools flagstat {input.bam} > {output}
        """

rule insert_size_picard:
    input:
        lambda w: samples_table_2.loc[w.sample, 'dedup_flt_both_strds_bam'],
    output:
        tab = "Results/{genomes_not_fused}/Qctrl/{sample}/Processed_bam/"
              "{sample}.insert_size_picard.tab",
        pdf = "Results/{genomes_not_fused}/Qctrl/{sample}/Processed_bam/"
              "{sample}.insert_size_picard.pdf",
    envmodules:
        config["picard"],
        config["R"],
    benchmark:
        "benchmarks/insert_size_picard/{genomes_not_fused}/"
        "{sample}.tsv",
    shell:
        """
        picard CollectInsertSizeMetrics M=0.25 I={input} O={output.tab} \\
          H={output.pdf} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true \\
          VERBOSITY=INFO
        """

rule library_complexity_picard:
    input:
        lambda w: samples_table_2.loc[w.sample, 'raw_bam'],
    output:
        "Results/{genomes_not_fused}/Qctrl/{sample}/Raw_bam/"
        "{sample}.picard_library_complexity.tab",
    envmodules:
        config["picard"],
    benchmark:
        "benchmarks/library_complexity_picard/{genomes_not_fused}/"
        "{sample}.tsv",
    shell:
        """
        picard EstimateLibraryComplexity I={input} O={output} \\
          VALIDATION_STRINGENCY=LENIENT VERBOSITY=INFO
        """

rule multiqc:
    input:
        smkf.multiqc_input,
    output:
        directory(
            "Results/{genomes_not_fused}/Qctrl/multiqc_report_{library_name}_data"
        ),
        html = "Results/{genomes_not_fused}/Qctrl/"
               "multiqc_report_{library_name}.html",
    envmodules:
        config["multiqc"],
    benchmark:
        "benchmarks/multiqc/{genomes_not_fused}/"
        "{library_name}.tsv",
    shell:
        """
        config_file_path="$(pwd)/Config/multiqc_config.yaml"
        cd Results/{wildcards.genomes_not_fused}/Qctrl && \\
        multiqc --dirs --dirs-depth 1 \\
          --filename "multiqc_report_{wildcards.library_name}.html" \\
          --config $config_file_path ./ ../../FASTQ_reports
        """
