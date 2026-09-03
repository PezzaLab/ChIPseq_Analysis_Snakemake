# Recombination hotspots, deepTools matrix, heatmaps, R aggregate profile wrangling, and HTML report

localrules: intersect_peaks_HSs_list

rule intersect_peaks_HSs_list:
    input:
        unpack(smkf.intersect_peaks_HSs_list_input),
    output:
        "Results/{genomes_not_fused}/Peaks/MACS2/{peak_type}/"
        "{peak_params}/blacklist_filtered/Intersect_HSs_plus_minus_2000_bp/"
        "{sample}.{genomes_final}{extension}.{peak_type}Peak",
    log:
        "logs/intersect_peaks_HSs_list/{genomes_not_fused}_{peak_type}_{peak_params}/"
        "{sample}.{genomes_final}{extension}.log",
    threads: 1
    envmodules:
        config["bedtools"],
    shell:
        """
        bedtools intersect -wa -a {input.sample_peaks} -b {input.hotspots} > {output} 2> {log}
        """

rule compute_matrix_outfiles_hs:
    input:
        unpack(smkf.compute_matrix_outfiles_hs_input),
    output:
        filename = temp(
            "Results/{genomes_not_fused}/Analysis/"
            "Heatmaps_and_aggregate_profiles/Hotspots/{hs_region}/{strands}/"
            "{cov_params}/filenames/"
            "{sample}.{genomes_final}.q_filt.srt.nodup.mit_filt."
            "{ss_condit}{clip_condit}filename"
        ),
        matrix = "Results/{genomes_not_fused}/Analysis/"
                 "Heatmaps_and_aggregate_profiles/Hotspots/{hs_region}/{strands}/"
                 "{cov_params}/matrixes/"
                 "{sample}.{genomes_final}.q_filt.srt.nodup.mit_filt."
                 "{ss_condit}{clip_condit}matrix",
    log:
        "logs/compute_matrix_outfiles_hs/{genomes_not_fused}_{hs_region}_{strands}_{cov_params}/"
        "{sample}.{genomes_final}.{ss_condit}{clip_condit}.log",
    threads: 4
    params:
        before_region = config['heatmaps']['before_region'],
        after_region = config['heatmaps']['after_region'],
        bin_size = config['heatmaps']['bin_size'],
    envmodules:
        config["deeptools"],
    benchmark:
        "benchmarks/compute_matrix_outfiles_hs/{genomes_not_fused}_"
        "{hs_region}_{strands}_{cov_params}/"
        "{sample}.{genomes_final}.{ss_condit}{clip_condit}.tsv",
    shell:
        """
        computeMatrix reference-point -p {threads} -R {input.region} \\
          -S {input.bigwig} \\
          --sortRegions keep \\
          -b {params.before_region} \\
          -a {params.after_region} \\
          --outFileName {output.filename} \\
          --outFileNameMatrix {output.matrix} \\
          --referencePoint center \\
          --binSize {params.bin_size} > {log} 2>&1
        """

rule process_aggregate_profiles:
    input:
        unpack(smkf.process_aggregate_profiles_inputs),
    output:
        "Results/{genomes_not_fused}/Analysis/Heatmaps_and_aggregate_profiles/"
        "Hotspots/{library_name}_{hs_region}.{smooth}.RData",
    log:
        "logs/process_aggregate_profiles/{genomes_not_fused}/"
        "{library_name}_{hs_region}.{smooth}.log",
    threads: 1
    envmodules:
        config["R"],
        config["bioconductor"],
    params:
        genomes = config['genomes'].keys(),
    benchmark:
        "benchmarks/process_aggregate_profiles/{genomes_not_fused}/"
        "{library_name}_{hs_region}.{smooth}.tsv",
    script:
        "../Scripts/process_outfile_matrix.R"

rule wrangle_X_nonPAR_asymmetric_HS:
    input:
        "Results/{genomes_not_fused}/Analysis/Heatmaps_and_aggregate_profiles/"
        "Hotspots/{library_name}_x_non_par.{smooth}.RData",
        "Results/{genomes_not_fused}/Analysis/Heatmaps_and_aggregate_profiles/"
        "Hotspots/{library_name}_autosomal_x_non_par_ctrl.{smooth}.RData",
        "Results/{genomes_not_fused}/Analysis/Heatmaps_and_aggregate_profiles/"
        "Hotspots/{library_name}_asymetric_watson_strong.{smooth}.RData",
        "Results/{genomes_not_fused}/Analysis/Heatmaps_and_aggregate_profiles/"
        "Hotspots/{library_name}_asymetric_crick_strong.{smooth}.RData",
        "Results/{genomes_not_fused}/Analysis/Heatmaps_and_aggregate_profiles/"
        "Hotspots/{library_name}_top_5000_plus_minus_2000.{smooth}.RData",
    output:
        "Results/{genomes_not_fused}/Analysis/Heatmaps_and_aggregate_profiles/"
        "Hotspots/{library_name}_aggregate_profiles_data.{smooth}.RData",
    log:
        "logs/wrangle_X_nonPAR_asymmetric_HS/{genomes_not_fused}/"
        "{library_name}_{smooth}.log",
    threads: 1
    envmodules:
        config["R"],
        config["bioconductor"],
    benchmark:
        "benchmarks/wrangle_X_nonPAR_asymmetric_HS/{genomes_not_fused}/"
        "{library_name}_{smooth}.tsv",
    script:
        "../Scripts/wrangle_X_nonPAR_asymmetric_HS.R"

rule process_aggregate_profiles_clipped:
    input:
        smkf.process_aggregate_profiles_clipped_input,
    output:
        "Results/{genomes_not_fused}/Analysis/Heatmaps_and_aggregate_profiles/"
        "Hotspots/{library_name}_{hs_region}_clipped.{smooth}.RData",
    log:
        "logs/process_aggregate_profiles_clipped/{genomes_not_fused}/"
        "{library_name}_{hs_region}_{smooth}.log",
    threads: 1
    envmodules:
        config["R"],
        config["bioconductor"],
    benchmark:
        "benchmarks/process_aggregate_profiles_clipped/{genomes_not_fused}/"
        "{library_name}_{hs_region}_{smooth}.tsv",
    script:
        "../Scripts/process_outfile_matrix_clipped.R"

rule plot_heatmaps:
    input:
        "Results/{genomes_not_fused}/Analysis/Heatmaps_and_aggregate_profiles/"
        "Hotspots/{hs_region}/Both_strands/{cov_params}/filenames/"
        "{sample}.{genomes_final}{extension}{strand}.filename",
    output:
        "Results/{genomes_not_fused}/Analysis/Heatmaps_and_aggregate_profiles/"
        "Hotspots/{hs_region}/Both_strands/{cov_params}/Heatmaps/"
        "{sample}.{genomes_final}{extension}{strand}.png",
    log:
        "logs/plot_heatmaps/{genomes_not_fused}_{hs_region}_{cov_params}/"
        "{sample}.{genomes_final}{extension}{strand}.log",
    threads: 1
    wildcard_constraints:
        extension = "(\\.q_filt\\.srt\\.nodup\\.mit_filt)?",
    envmodules:
        config["deeptools"],
    benchmark:
        "benchmarks/plot_heatmaps/{genomes_not_fused}_{hs_region}_"
        "{cov_params}/{sample}.{genomes_final}{extension}{strand}.tsv",
    shell:
        """
        plotHeatmap --matrixFile {input} --sortRegions keep \\
          --plotTitle "{wildcards.sample}, {wildcards.strand} at {wildcards.hs_region}" \\
          -out {output} \\
          --legendLocation none > {log} 2>&1
        """

rule markdown_report_aggregate_profiles:
    input:
        unpack(smkf.markdown_report_aggregate_profiles_input),
    output:
        html_report = "Results/{genomes_not_fused}/Analysis/{lib_name}.{genomes_not_fused}.{smooth}.html",
        xlsx_table = "Results/{genomes_not_fused}/Analysis/"
                     "{lib_name}_for_NGS_files.{smooth}.xlsx",
    log:
        "logs/markdown_report_aggregate_profiles/{genomes_not_fused}/{lib_name}.{smooth}.log",
    threads: 1
    envmodules:
        config["R"],
        config["bioconductor"],
    benchmark:
        "benchmarks/markdown_report_aggregate_profiles/{genomes_not_fused}/"
        "{lib_name}.{genomes_not_fused}.{smooth}.tsv",
    script:
        "../Scripts/html_report.R"
