# Spike-in separation, deduplication, and quality filtering rules

rule filter_spiked_in_bam_int_genome:
    input:
        bam="Results/{sample}.{fused_genome}.bam",
        bai="Results/{sample}.{fused_genome}.bam.bai",
    output:
        "Results/{sample}.{fused_genome}.{int_genome}.bam",
    log:
        "logs/filter_spiked_in_bam_int_genome/{sample}.{fused_genome}.{int_genome}.log",
    threads: 8
    params:
        spike_suffix= lambda w: config['spike_in_suffixes'][w.fused_genome],
    envmodules:
        config["samtools"],
        config["sambamba"],
    benchmark:
        "benchmarks/filter_spiked_in_bam_int_genome/{sample}.{fused_genome}.{int_genome}.tsv",
    shell:
        """
        (int_chromosomes=$(samtools idxstats {input.bam} | \\
         cut -f 1 | grep -v .*{params.spike_suffix} | \\
         grep -v \\* | paste -sd ' ')
        sambamba view --nthreads {threads} --format bam \\
         -o {output} {input.bam} $int_chromosomes) > {log} 2>&1
        """

rule filter_spiked_in_bam_spike_in_genome:
    input:
        bam="Results/{sample}.{fused_genome}.bam",
        bai="Results/{sample}.{fused_genome}.bam.bai",
    output:
        "Results/{sample}.{fused_genome}.{spike_genome}.bam",
    log:
        "logs/filter_spiked_in_bam_spike_in_genome/{sample}.{fused_genome}.{spike_genome}.log",
    threads: 8
    params:
        spike_suffix= lambda w: config['spike_in_suffixes'][w.fused_genome],
    envmodules:
        config["sambamba"],
        config["samtools"],
    benchmark:
        "benchmarks/filter_spiked_in_bam_spike_in_genome/{sample}.{fused_genome}.{spike_genome}.tsv",
    shell:
        """
        (spike_chromosomes=$(samtools idxstats {input.bam} | cut -f 1 | \\
         grep .*{params.spike_suffix} | grep -v \\* | paste -sd ' ')

        sambamba view \\
          --nthreads {threads} \\
          --format bam \\
          -o {output} \\
          {input.bam} \\
          $spike_chromosomes) > {log} 2>&1
        """

rule dedup_bam:
    input:
        bam = "Results/{sample}.{genomes_final}.bam",
        bai = "Results/{sample}.{genomes_final}.bam.bai",
    output:
        bam = temp(
            "Results/{genomes_not_fused}/{sample}.{genomes_final}.nodup.bam"
        ),
        qctrl = "Results/{genomes_not_fused}/Qctrl/{sample}/Raw_bam/"
                "{sample}.{genomes_final}.markDup_raw.txt",
        bai = temp(
            "Results/{genomes_not_fused}/{sample}.{genomes_final}.nodup.bai"
        ),
    log:
        "logs/dedup_bam/{genomes_not_fused}/{sample}.{genomes_final}.log",
    threads: 2
    envmodules:
        config["picard"],
    benchmark:
        "benchmarks/dedup_bam/{genomes_not_fused}/{genomes_final}/{sample}.{genomes_final}.tsv",
    shell:
        """
        picard MarkDuplicates I={input.bam} O={output.bam} M={output.qctrl} \\
          VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true \\
          REMOVE_DUPLICATES=TRUE VERBOSITY=INFO CREATE_INDEX=TRUE \\
          > {log} 2>&1
        """

rule filter_bam:
    input:
        bam = "Results/{genomes_not_fused}/{sample}.{genomes_final}.nodup.bam",
        bai = "Results/{genomes_not_fused}/{sample}.{genomes_final}.nodup.bam.bai",
        chromosomes = "Resources/bed_files/{genomes_not_fused}/{genomes_not_fused}_chromosomes.txt",
    output:
        "Results/{genomes_not_fused}/Bams/Both_strands/"
        "{sample}.{genomes_final}.q_filt.srt.nodup.mit_filt.bam",
    log:
        "logs/filter_bam/{genomes_not_fused}/{sample}.{genomes_final}.log",
    threads: 8
    params:
        min_mapq = config['min_mapping_qual'],
        flags = smkf.filter_bam_params,
    envmodules:
        config["samtools"],
        config["sambamba"],
    benchmark:
        "benchmarks/filter_bam/{genomes_not_fused}/{sample}.{genomes_final}.tsv",
    shell:
        r"""
        # Filter out mitochondrial and unplaced/unlocalized scaffolds
        # Keep only unique properly mapped reads (see flags)
        samtools view \
          -q {params.min_mapq} \
          -b \
          -@ {threads} \
          {params.flags} \
          -o {output} \
          {input.bam} \
          $(cat {input.chromosomes}) \
          > {log} 2>&1
        """
