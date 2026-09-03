# Peak calling and filtering rules (MACS2, bedtools blacklist filtering, and summary)

localrules: filter_peaks_blk_grey_list, sumarize_peak_count

rule call_peaks_macs2:
    input:
        unpack(smkf.call_peaks_macs2_input)
    output:
        "Results/{genomes_not_fused}/Peaks/MACS2/{peak_type}/{peak_params}/"
        "{sample}.{genomes_final}{extension}_peaks.{peak_type}Peak",
    params:
        smkf.call_peaks_macs2_params,
    wildcard_constraints:
        extension = r"(\.q_filt\.srt\.nodup\.mit_filt)?",
    envmodules:
        config["macs2"],
    benchmark:
        "benchmarks/call_peaks_macs2/{genomes_not_fused}_{peak_type}_{peak_params}/"
        "{sample}.{genomes_final}{extension}_peaks.{peak_type}Peak.tsv",
    shell:
        """
        directory=$(dirname {output})
        macs2 callpeak --outdir $directory \\
         --name {wildcards.sample}.{wildcards.genomes_final}{wildcards.extension}\\
         -t {input.treat_bam} {params[0][ctrl]} \\
         -g {params[0][genome]} \\
         {params[0][qv_bco]} {params[0][PE]} {params[0][peak_type_options]} {params[0][extension]}
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
    wildcard_constraints:
        extension = r"(\.q_filt\.srt\.nodup\.mit_filt)?",
    envmodules:
        config["bedtools"],
    shell:
        """
        bedtools intersect -v -a {input.peaks} -b {input.black_list} | \\
        awk 'BEGIN {{ OFS = "\\t"}} {{print $0, ($1 ":" $2 "-" $3)}}' - > \\
        {output.blacklist_filtered}

        bedtools intersect -wa -u -a {input.peaks} -b {input.black_list} | \\
        awk 'BEGIN {{ OFS = "\\t"}} {{print $0, ($1 ":" $2 "-" $3)}}' - > \\
        {output.blacklist_overlap}
        """

rule sumarize_peak_count:
    input:
        unpack(smkf.sumarize_peak_count_input),
    output:
        "Results/{genomes_not_fused}/Analysis/Peaks_summary.tsv",
    params:
        samples_table_2 = samples_table_2,
    script:
        "../Scripts/peaks_summary.py"
