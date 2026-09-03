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
    log:
        "logs/align_fastq/{sample}.{genomes_all}.log",
    threads: 12
    envmodules:
        config["bwa"],
        config["samtools"],
    benchmark:
        "benchmarks/align_fastq/{sample}.{genomes_all}.tsv",
    shell:
        """
        (bwa mem -v 1 -M \\
          -t {threads} {params.genome_prefix} \\
          {input.fq1} {input.fq2} |\\
        samtools sort -@ {threads} -O bam \\
          -o {output.bam}) > {log} 2>&1
        """

rule index:
    input:
        "{sample}.bam",
    output:
        "{sample,.+}.bam.bai",
    log:
        "logs/index/{sample}.log",
    threads: 2
    envmodules:
        config['sambamba'],
    shell:
        """
        sambamba index -t {threads} {input} > {log} 2>&1
        """

rule merge_bams:
    input:
        smkf.merge_bams_input,
    output:
        temp("Results/{sample}_MERGED.{genomes_all}.bam"),
    log:
        "logs/merge_bams/{sample}.{genomes_all}.log",
    threads: 2
    envmodules:
        config["samtools"],
    benchmark:
        "benchmarks/merge_bams/{sample}.{genomes_all}.tsv",
    shell:
        """
        samtools merge -@ {threads} -o {output} {input} > {log} 2>&1
        """
