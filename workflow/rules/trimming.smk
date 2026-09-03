# Trimming and adapter removal rules (fastp and cutadapt)

rule trim_adapters_SE:
    input:
        lambda w: samples_table.loc[w.sample, "fastq1"],
    output:
        fq1 = temp("Results/{sample}.adap_trimmed.R1.SE.fq.gz"),
        html = "Results/FASTQ_reports/{sample}.SE.fastp.html",
        json = "Results/FASTQ_reports/{sample}.SE.fastp.json",
    envmodules:
        config["fastp"],
    benchmark:
        "benchmarks/trim_adapters_SE/{sample}.tsv",
    shell:
        """
        fastp \\
          --thread $SLURM_CPUS_ON_NODE \\
          -i {input} \\
          --html {output.html} \\
          --json {output.json} \\
          -o {output.fq1} \\

          # For SE data, the adapters are evaluated by analyzing the
          # tails of first ~1M reads. FASTP has a built-in list of the
          # most common adapters used. If adapters are 'weird', you can
          # specify them using --adapter_sequence option. For this the
          # snakemake rule should be changed.
        """

rule trim_adapters_PE:
    input:
        unpack(smkf.trim_adapters_PE_input),
    output:
        fq1 = temp("Results/{sample}.adap_trimmed.R1.PE.fq.gz"),
        fq2 = temp("Results/{sample}.adap_trimmed.R2.PE.fq.gz"),
        html = "Results/FASTQ_reports/{sample}.PE.fastp.html",
        json = "Results/FASTQ_reports/{sample}.PE.fastp.json",
    envmodules:
        config["fastp"],
    benchmark:
        "benchmarks/trim_adapters_PE/{sample}.tsv",
    shell:
        """
        fastp \\
          --thread $SLURM_CPUS_ON_NODE \\
          -i {input.fastq1} \\
          -I {input.fastq2} \\
          --html {output.html} \\
          --json {output.json} \\
          -o {output.fq1} \\
          -O {output.fq2}
        """

rule trim_adaptase:
    input:
        "Results/{sample}.adap_trimmed.R{read}.{n_ends}.fq.gz",
    output:
        fq1 = temp("Results/{sample}.adap_trimmed.adaptase_trimmed.R{read}.{n_ends}.fq.gz"),
    envmodules:
        config["cutadapt"],
    benchmark:
        "benchmarks/trim_adaptase/{sample}.{read}.{n_ends}.tsv",
    shell:
        """
        cutadapt \\
          --cores 0 \\
          -u 10  \\
          -o {output} \\
          {input}

        # -u : remove bases from the beginning or end of each read. If
        #      the given length is positive, the bases are removed from
        #      the beginning of each read.
        # --cores: Use 0 to auto-detect the number of available cores
        """
