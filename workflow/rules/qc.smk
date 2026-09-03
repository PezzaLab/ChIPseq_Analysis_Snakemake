# Quality control rules: Picard, samtools stats/flagstat, MultiQC, and FRIP

rule samstats:
    input:
        unpack(smkf.samstats_samtools_flagstat_input),
    output:
        "Results/{genomes_not_fused}/Qctrl/{sample}/{bam_type}/"
        "{sample}.samstats.txt",
    log:
        "logs/samstats/{genomes_not_fused}_{bam_type}/{sample}.log",
    threads: 2
    envmodules:
        config["samtools"],
    benchmark:
        "benchmarks/samstats/{genomes_not_fused}_{bam_type}/"
        "{sample}.tsv",
    shell:
        """
        samtools stats -@ {threads} {input.bam} > {output} 2> {log}
        """

rule samtools_flagstat:
    input:
        unpack(smkf.samstats_samtools_flagstat_input),
    output:
        "Results/{genomes_not_fused}/Qctrl/{sample}/{bam_type}/"
        "{sample}.flagstat.txt",
    log:
        "logs/samtools_flagstat/{genomes_not_fused}_{bam_type}/{sample}.log",
    threads: 2
    envmodules:
        config["samtools"],
    benchmark:
        "benchmarks/samtools_flagstat/{genomes_not_fused}_{bam_type}/"
        "{sample}.tsv",
    shell:
        """
        samtools flagstat -@ {threads} {input.bam} > {output} 2> {log}
        """

rule insert_size_picard:
    input:
        lambda w: samples_table_2.loc[w.sample, 'dedup_flt_both_strds_bam'],
    output:
        tab = "Results/{genomes_not_fused}/Qctrl/{sample}/Processed_bam/"
              "{sample}.insert_size_picard.tab",
        pdf = "Results/{genomes_not_fused}/Qctrl/{sample}/Processed_bam/"
              "{sample}.insert_size_picard.pdf",
    log:
        "logs/insert_size_picard/{genomes_not_fused}/{sample}.log",
    threads: 2
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
          VERBOSITY=INFO > {log} 2>&1
        """

rule library_complexity_picard:
    input:
        lambda w: samples_table_2.loc[w.sample, 'raw_bam'],
    output:
        "Results/{genomes_not_fused}/Qctrl/{sample}/Raw_bam/"
        "{sample}.picard_library_complexity.tab",
    log:
        "logs/library_complexity_picard/{genomes_not_fused}/{sample}.log",
    threads: 2
    envmodules:
        config["picard"],
    benchmark:
        "benchmarks/library_complexity_picard/{genomes_not_fused}/"
        "{sample}.tsv",
    shell:
        """
        picard EstimateLibraryComplexity I={input} O={output} \\
          VALIDATION_STRINGENCY=LENIENT VERBOSITY=INFO > {log} 2>&1
        """

rule FRIP:
    input:
        unpack(smkf.FRIP_input),
    output:
        "Results/{genomes_not_fused}/Qctrl/{sample}/Processed_bam/FRIP/"
        "MACS2_{peak_type}_{peak_params}_bl-gr_flt/"
        "{sample}.{genomes_final}{extension}.FRIP.txt",
    wildcard_constraints:
        peak_type = "broad|narrow",
        extension = r"(\.q_filt\.srt\.nodup\.mit_filt)?",
    log:
        "logs/FRIP/{genomes_not_fused}/{sample}.{genomes_final}{extension}.{peak_type}_{peak_params}.log",
    threads: 2
    envmodules:
        config["sambamba"],
        config["samtools"],
    benchmark:
        "benchmarks/FRIP/{genomes_not_fused}/{sample}.{genomes_final}{extension}.{peak_type}_{peak_params}.tsv",
    shell:
        """
        (
        total_reads=$(sambamba view -t {threads} -c {input.bam})
        peak_reads=$(samtools view -@ {threads} -L {input.peak} -c {input.bam})

        if [ "$total_reads" -gt 0 ]; then
            FRIP=$(awk -v pr="$peak_reads" -v tr="$total_reads" 'BEGIN {{ printf "%.2f", 100 * pr / tr }}')
        else
            FRIP="0.00"
        fi

        echo "Sample name: $(basename {input.bam} .bam)" > {output}
        echo "Total reads: $total_reads" >> {output}
        echo "Reads in peaks: $peak_reads" >> {output}
        echo "FRIP (%): $FRIP" >> {output}
        ) > {log} 2>&1
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
    log:
        "logs/multiqc/{genomes_not_fused}/{library_name}.log",
    threads: 2
    params:
        config_file = os.path.abspath("Config/multiqc_config.yaml"),
    envmodules:
        config["multiqc"],
    benchmark:
        "benchmarks/multiqc/{genomes_not_fused}/"
        "{library_name}.tsv",
    shell:
        """
        (cd Results/{wildcards.genomes_not_fused}/Qctrl && \\
        multiqc --dirs --dirs-depth 1 \\
          --filename "multiqc_report_{wildcards.library_name}.html" \\
          --config {params.config_file} ./ ../../FASTQ_reports) > {log} 2>&1
        """
