# Trimming and adapter removal rules (fastp and cutadapt)

rule trim_adapters_SE:
    input:
        lambda w: samples_table.loc[w.sample, "fastq1"],
    output:
        fq1 = temp("Results/{sample}.adap_trimmed.R1.SE.fq.gz"),
        html = "Results/FASTQ_reports/{sample}.SE.fastp.html",
        json = "Results/FASTQ_reports/{sample}.SE.fastp.json",
    log:
        "logs/trim_adapters_SE/{sample}.log",
    threads: 6
    envmodules:
        config["fastp"],
    benchmark:
        "benchmarks/trim_adapters_SE/{sample}.tsv",
    shell:
        """
        fastp \\
          --thread {threads} \\
          -i {input} \\
          --html {output.html} \\
          --json {output.json} \\
          -o {output.fq1} \\
          > {log} 2>&1
        """

rule trim_adapters_PE:
    input:
        unpack(smkf.trim_adapters_PE_input),
    output:
        fq1 = temp("Results/{sample}.adap_trimmed.R1.PE.fq.gz"),
        fq2 = temp("Results/{sample}.adap_trimmed.R2.PE.fq.gz"),
        html = "Results/FASTQ_reports/{sample}.PE.fastp.html",
        json = "Results/FASTQ_reports/{sample}.PE.fastp.json",
    log:
        "logs/trim_adapters_PE/{sample}.log",
    threads: 6
    envmodules:
        config["fastp"],
    benchmark:
        "benchmarks/trim_adapters_PE/{sample}.tsv",
    shell:
        """
        fastp \\
          --thread {threads} \\
          -i {input.fastq1} \\
          -I {input.fastq2} \\
          --html {output.html} \\
          --json {output.json} \\
          -o {output.fq1} \\
          -O {output.fq2} \\
          > {log} 2>&1
        """

rule trim_adaptase:
    input:
        "Results/{sample}.adap_trimmed.R{read}.{n_ends}.fq.gz",
    output:
        fq1 = temp("Results/{sample}.adap_trimmed.adaptase_trimmed.R{read}.{n_ends}.fq.gz"),
    log:
        "logs/trim_adaptase/{sample}.{read}.{n_ends}.log",
    threads: 4
    envmodules:
        config["cutadapt"],
    benchmark:
        "benchmarks/trim_adaptase/{sample}.{read}.{n_ends}.tsv",
    shell:
        """
        cutadapt \\
          --cores {threads} \\
          -u 10 \\
          -o {output} \\
          {input} \\
          > {log} 2>&1
        """
