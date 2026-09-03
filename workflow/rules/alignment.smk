# Alignment, indexing, and BAM merging rules

rule align_fastq:
    input:
        unpack(smkf.align_fastq_input),
    params:
        # Take one of the index files (e.g., .../genome.fa.bwt) and drop only the final ext → .../genome.fa
        genome_prefix = lambda w, input: re.sub(
            r'\.(sa|pac|bwt|ann|amb)$', '', input.reference_genome_indexed_files[0]
        ),
    output:
        bam = temp("Results/{sample}.{genomes_all}.bam"),
    envmodules:
        config["bwa"],
        config["samtools"],
    benchmark:
        "benchmarks/align_fastq/{sample}.{genomes_all}.tsv",
    shell:
        """
        bwa mem -v 1 -M \\
          -t $SLURM_CPUS_ON_NODE {params.genome_prefix} \\
          {input.fq1} {input.fq2} |\\
        samtools sort -@ $SLURM_CPUS_ON_NODE -O bam \\
          -o {output.bam}

        # bwa
            # -M: Mark shorter split hits as secondary (for Picard
            #     compatibility)
        """

rule index:
    input:
        "{sample}.bam",
    output:
        "{sample,.+}.bam.bai",
    envmodules:
        config['sambamba'],
    shell:
        """
        sambamba index {input}
        """

rule merge_bams:
    input:
        smkf.merge_bams_input,
    output:
        temp("Results/{sample}_MERGED.{genomes_all}.bam"),
    envmodules:
        config["samtools"],
    benchmark:
        "benchmarks/merge_bams/{sample}"
        "{sample}.{genomes_all}.tsv"
    shell:
        """
        samtools merge -o {output} {input}
        """
