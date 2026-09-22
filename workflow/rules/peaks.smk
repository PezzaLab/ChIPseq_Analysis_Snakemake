# Peak calling and filtering rules (MACS2, bedtools blacklist filtering, and summary)

localrules: filter_peaks_blk_grey_list, summarize_peak_count

rule call_peaks_macs2:
    input:
        unpack(smkf.call_peaks_macs2_input),
    output:
        "Results/{genomes_not_fused}/Peaks/MACS2/{peak_type}/{peak_params}/"
        "{sample}.{genomes_final}{extension}_peaks.{peak_type}Peak",
    log:
        "logs/call_peaks_macs2/{genomes_not_fused}_{peak_type}_{peak_params}/"
        "{sample}.{genomes_final}{extension}.log",
    threads: 2
    params:
        p = smkf.call_peaks_macs2_params,
    wildcard_constraints:
        extension = r"(\.q_filt\.srt\.nodup\.mit_filt)?",
    envmodules:
        config["macs2"],
    benchmark:
        "benchmarks/call_peaks_macs2/{genomes_not_fused}_{peak_type}_{peak_params}/"
        "{sample}.{genomes_final}{extension}_peaks.{peak_type}Peak.tsv",
    shell:
        """
        (directory=$(dirname {output})
        macs2 callpeak --outdir $directory \\
         --name {wildcards.sample}.{wildcards.genomes_final}{wildcards.extension} \\
         -t {input.treat_bam} {params.p[ctrl]} \\
         -g {params.p[genome]} \\
         {params.p[qv_bco]} {params.p[PE]} {params.p[peak_type_options]} {params.p[extension]}) > {log} 2>&1
        """

rule filter_peaks_blk_grey_list:
    input:
        unpack(smkf.filter_peaks_blk_grey_list_input),
    output:
        blacklist_filtered = "Results/{genomes_not_fused}/Peaks/MACS2/{peak_type}/"
                             "{peak_params}/blacklist_filtered/{sample}."
                             "{genomes_final}{extension}.{peak_type}Peak",
        blacklist_overlap = "Results/{genomes_not_fused}/Peaks/MACS2/{peak_type}/"
                            "{peak_params}/blacklist_overlap/{sample}."
                            "{genomes_final}{extension}.{peak_type}Peak",
    log:
        "logs/filter_peaks_blk_grey_list/{genomes_not_fused}_{peak_type}_{peak_params}/"
        "{sample}.{genomes_final}{extension}.log",
    threads: 1
    wildcard_constraints:
        extension = r"(\.q_filt\.srt\.nodup\.mit_filt)?",
    envmodules:
        config["bedtools"],
    shell:
        """
        (bedtools intersect -v -a {input.peaks} -b {input.black_list} | \\
        awk 'BEGIN {{ OFS = "\\t"}} {{print $0, ($1 ":" $2 "-" $3)}}' - > \\
        {output.blacklist_filtered}

        bedtools intersect -wa -u -a {input.peaks} -b {input.black_list} | \\
        awk 'BEGIN {{ OFS = "\\t"}} {{print $0, ($1 ":" $2 "-" $3)}}' - > \\
        {output.blacklist_overlap}) 2> {log}
        """

rule summarize_peak_count:
    input:
        unpack(smkf.summarize_peak_count_input),
    output:
        "Results/{genomes_not_fused}/Analysis/Peaks_summary.tsv",
    log:
        "logs/summarize_peak_count/{genomes_not_fused}.log",
    threads: 1
    params:
        samples_table_2 = samples_table_2,
    script:
        "../Scripts/peaks_summary.py"
