# Single-strand read separation, Watson/Crick merging, and 1bp-clipping rules

rule get_strand_sep_bams:
    input:
        bam="Results/{genomes_not_fused}/Bams/Both_strands/"
            "{sample}.{genomes_final}.q_filt.srt.nodup.mit_filt.bam",
    output:
        bam=temp(
            "Results/{genomes_not_fused}/Bams/Single_strand/"
            "{read_length}/{sample}.{genomes_final}."
            "q_filt.srt.nodup.mit_filt.{flag}.bam"
        ),
        bai=temp(
            "Results/{genomes_not_fused}/Bams/Single_strand/"
            "{read_length}/{sample}.{genomes_final}."
            "q_filt.srt.nodup.mit_filt.{flag}.bam.bai"
        ),
    params:
        smkf.get_strand_sep_bams_params
    envmodules:
        config["samtools"],
    benchmark:
        "benchmarks/get_strand_sep_bams/{genomes_not_fused}_{read_length}/"
        "{sample}.{genomes_final}.{flag}.tsv",
    shell:
        """
        samtools view -b {params} {input.bam} > {output.bam}
        samtools index {output.bam}
        """

rule merge_watson_or_crick:
    input:
        bam_1="Results/{genomes_not_fused}/Bams/Single_strand/"
              "Full_length_reads/{sample}.{genomes_final}."
              "q_filt.srt.nodup.mit_filt.{flag1}.bam",
        index_1="Results/{genomes_not_fused}/Bams/Single_strand/"
                "Full_length_reads/{sample}.{genomes_final}."
                "q_filt.srt.nodup.mit_filt.{flag1}.bam.bai",
        bam_2="Results/{genomes_not_fused}/Bams/Single_strand/"
              "Full_length_reads/{sample}.{genomes_final}."
              "q_filt.srt.nodup.mit_filt.{flag2}.bam",
        index_2="Results/{genomes_not_fused}/Bams/Single_strand/"
                "Full_length_reads/{sample}.{genomes_final}."
                "q_filt.srt.nodup.mit_filt.{flag2}.bam.bai",
    output:
        "Results/{genomes_not_fused}/Bams/Single_strand/Full_length_reads/"
        "{sample}.{genomes_final}.q_filt.srt.nodup.mit_filt."
        "{flag1}-{flag2}.bam",
    wildcard_constraints:
        flag1="(83|99)",
        flag2="(163|147)",
    envmodules:
        config["samtools"],
    benchmark:
        "benchmarks/merge_watson_or_crick/{genomes_not_fused}/"
        "{sample}.{genomes_final}.{flag1}-{flag2}.tsv",
    shell:
        """
        samtools merge {output} {input.bam_1} {input.bam_2}
        """

rule clip_1bp:
    input:
        smkf.clip_1bp_input,
    output:
        bam="Results/{genomes_not_fused}/Bams/Single_strand/1bp_clipped_reads/"
            "{sample}.{genomes_final}.q_filt.srt.nodup.mit_filt.{ss_PE}{ss_SR}"
            ".clipped_1_bp.bam",
        tmp_dir=temp(
            directory(
                "Results/{genomes_not_fused}/Bams/Single_strand/1bp_clipped_reads/"
                "{sample}.{genomes_final}.q_filt.srt.nodup.mit_filt.{ss_PE}{ss_SR}/"
                "temp_clip_1bp"
            )
        ),
    wildcard_constraints:
        flag16="inc_16|exc_16",
    params:
        flag = smkf.clip_1bp_param,
    envmodules:
        config["samtools"],
        config["bamutil"],
    benchmark:
        "benchmarks/clip_1bp/{genomes_not_fused}/"
        "{sample}.{genomes_final}.q_filt.srt.nodup.mit_filt.{ss_PE}{ss_SR}.tsv",
    shell:
        """
        set -x
        temp_dir="{output.tmp_dir}"

        mkdir -p $temp_dir

        # Get header and sam file
        samtools view -H {input} > ${{temp_dir}}/header

        # Clip reads
        samtools view {input} -{params.flag} 16 | \\
        awk -v dir="${{temp_dir}}" '{{read_length=length($10)}} \\
          {{print $0 >> dir "/" read_length ".sam" }} '

        sam_files=$(find ./${{temp_dir}} -regextype awk -iregex ".*sam$")

        for sam in ${{sam_files}}; do
          read_length=$(basename $sam .sam)
          clip_length=$(( $read_length - 1 ))
          echo -e "clippping reads of readlength = $read_length \\n"
          cat ${{temp_dir}}/header $sam |
          bam trimBam - - -R $clip_length --clip | \
          samtools sort -n -o - - | \
          samtools fixmate - - | \
          samtools sort - -o ${{temp_dir}}/${{read_length}}_clipped_temp.bam
        done

        # Merge them
        bam_files=$(find ./${{temp_dir}} -regextype awk -iregex ".*bam$" -type f \\
            | tr "\\n" " ")
        samtools merge -o {output.bam} $bam_files
        """
