# Coverage BigWig generation: unstranded CPM & Drosophila spike-in normalization

rule create_coverage_bw:
    input:
        bam="Results/{genomes_not_fused}/Bams/{strands}/"
            "{sample}.{genomes_final}.q_filt.srt.nodup.mit_filt."
                "{ss_condit}{clip_condit}bam",
        bai="Results/{genomes_not_fused}/Bams/{strands}/"
            "{sample}.{genomes_final}.q_filt.srt.nodup.mit_filt."
                "{ss_condit}{clip_condit}bam.bai",
    output:
        "Results/{genomes_not_fused}/Bigwigs/Coverage/{strands}/"
        "{norm}_bs{bin_size}_sm{smooth_size}_ex{extension_size}/"
        "{sample}.{genomes_final}.q_filt.srt.nodup.mit_filt."
        "{ss_condit}{clip_condit}bw",
    envmodules:
        config["deeptools"],
    benchmark:
        "benchmarks/create_coverage_bw/{genomes_not_fused}_{strands}_{norm}_bs"
        "{bin_size}_sm{smooth_size}_ex{extension_size}/"
        "{sample}.{genomes_final}.q_filt.srt.nodup.mit_filt."
        "{ss_condit}{clip_condit}.tsv",
    shell:
        """
        # The coverage will be different depending on whether the bam has
        # been clipped or not. If clipped, then cov smooth and bin == 1
        # and no read extension is done, otherwise use config file parameters

        if [ -z "{wildcards.clip_condit}" ]
        then
            echo "Doing smooth and bin according to config file (not clipped profile)"
            bamCoverage --numberOfProcessors max -b {input.bam} \\
             -bs {wildcards.bin_size} \\
             --smoothLength {wildcards.smooth_size} \\
             --extendReads {wildcards.extension_size} \\
             --normalizeUsing {wildcards.norm} \\
             -o {output}
        else
            echo "Doing smooth 1 and bin 1 (clipped profile)"
            bamCoverage \\
             --numberOfProcessors max \\
             -b {input.bam} -bs 1 \\
             --smoothLength 1 \\
             --normalizeUsing {config[coverage][normalization]} \\
             -o {output}
        fi
        """

rule dros_normalization_report:
    input:
        smkf.dros_normalization_report_input,
    output:
        "Results/d6/Analysis/drosophila_normalization/drosophila_100K_reads/"
        "drosophila_equalization_report.tsv",
    params:
        samples_table_2 = samples_table_2,
    script:
        "../Scripts/dros_norm_report.py"

rule dros_normalization:
    input:
        unpack(smkf.dros_normalization_input),
    output:
        "Results/{genomes_not_fused}/Bigwigs/Coverage/{strands}/"
        "drosNormalized_bs{bin_size}_sm{smooth_size}_ex{extension_size}/"
        "drosophila_100K_reads/{sample}.{genomes_final}{extension}{strand}."
        "dros_norm.bw",
    envmodules:
        config["deeptools"],
    benchmark:
        "benchmarks/dros_normalization/{genomes_not_fused}_{strands}_"
        "{bin_size}_sm{smooth_size}_ex{extension_size}/"
        "{sample}.{genomes_final}{extension}{strand}.tsv",
    shell:
        """
        scaleFactor=$(awk 'BEGIN {{FS="\\t"}} $1 == "{wildcards.sample}" \\
         {{print $4}}' {input.report})

        if [ -z "${{scaleFactor}}" ]; then
            echo "Empty/non-existent scale factor."
            echo "Most likely there were no drosophila reads on your sample"
            exit 1
        fi

        echo -e "scaleFactor = $scaleFactor \\n"
        echo "Sample will be scaled down with a scale factor of $scaleFactor"

        bamCoverage --numberOfProcessors max -b {input.bam} \\
         -bs {wildcards.bin_size} \\
         --smoothLength {wildcards.smooth_size} \\
         --extendReads {wildcards.extension_size} \\
         -o {output} \\
         --scaleFactor $scaleFactor
        """
