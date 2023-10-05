# TODO: add module docstring

import pandas as pd
import re
import sys


def check_sample_table_format(samples_table):
    """Check if samples_table.csv has been properly filled.

    The following conditions should be met:
        * 'samples_names' do not contain '/', '.', ' ', or finish in '_MERGED'.
        * 'samples_names':'reference_genome' combinations should not be
            repeated.
        *  'merge_with' cannot be the path for a deduplicated/filterted bam
            file if the sample has drosophila spike-in (). This is because
            deduplication happens after separation of mm10 and drosohpila
            reads.
        * 'merge_with' must either be a sample_name contained within the
            samples_table, or an absolute path.
        * There must be 2 FASTQs when sample is PE and 1 when it is not.
        * FASTQ1 != FASTQ2
        * If 'dros_equalization_group' is not empty, then
            'dros_spike_in' == True
        * Booleans columns should have boolean values only
        * 'peak_ctrl_file_alias' should not contain '/', '.' or ' '

    Parameters:
    -----------
        samples_table : dataframe produced by reading samples_table.csv
            containing samples information.

    Return:
    -------
        None
    """

    exit_message = ("\n\nThe following error/s have been detected on your "
                    "'samples_table.csv' file:\n"
                    )
    exit_script = False

    # %% Boolean columns
    # Check that all columns with false-true are actually boolean
    cols_to_check = [
        "PE", "dros_spike_in", "get_single_strand",
        "Clip_reads_to_1bp_on_5_prime", "top5000_HS_heatmap"
    ]
    booleans_check = [
        col for col in cols_to_check if samples_table[col].dtype != "bool"
    ]

    if (len(booleans_check) > 0):
        exit_bool_message = (
            "* The following column/s has/ve at least one"
            " row  filled with something different than 'T', 'F', 'True',"
            " 'False', 'TRUE' or 'FALSE':\n"
            f"\t({booleans_check})"
        )
        sys.exit(exit_bool_message)
    # If any of this booleans is not boolean I need to exit now because
    # following checks use the boolean value of some of this columns

    # %% Samples names check
    # Check names of samples don't contain / or . or finish in "_MERGED"
    if (samples_table['sample_name'].str.contains(r"\.") |
            samples_table['sample_name'].str.contains("/") |
            samples_table['sample_name'].str.contains(" ") |
            samples_table['sample_name'].str.contains("_MERGED$")).any():
        exit_message += ("* One or more sample names in samples table contain"
                         " a dot ('.'), a slash ('/'), a space (' ') or ends"
                         " with the word '_MERGED'. These are not allowed in"
                         " sample names. Modify names and try again\n"
                         )
        exit_script = True
    # %% sample_name:genome unique
    # Check that there is no sample_name:genome combination repeated
    if samples_table.loc[:, ['sample_name',
                             'reference_genome'
                             ]
                         ].duplicated().any():
        exit_message += (
            "* One or more sample names in samples table is repeated.\n"
            "Choose different names for all your samples.\n"
        )
        exit_script = True
    # %% Merged samples
    if pd.notnull(samples_table['merge_with']).any():
        # Check that "merge_with" reference is not a deduplicated/filtered
        # file for spiked-in samples (it has to be raw bam)
        if (samples_table['merge_with'].str.contains(".nodup.") &
                samples_table['dros_spike_in']).any():
            exit_message += (
                "* One or more samples' 'merge_with' parameter is a "
                "deduplicated/filtered bam file AND has drosophila spike in. "
                "When sample is spiked, the 'merge_with' file has to be the "
                "raw bam, otherwise you'll loose the drosophila reads.\n"
            )
            exit_script = True

        # Check that to merge samples are either an absolute path or another
        # sample from samples_table
        samples_names_not_in_library = []
        for sample in samples_table['sample_name']:
            if not pd.isnull(samples_table.loc[sample, 'merge_with']):
                merge_samples = samples_table.loc[sample, 'merge_with'].split()
                for merge_sample in merge_samples:
                    absolute_path = merge_sample[0] == "/"
                    in_library = merge_sample in samples_table['sample_name']
                    if not absolute_path and not in_library:
                        samples_names_not_in_library.append(
                            f"\t{merge_sample}")
                        exit_script = True
        format_samples_not_in_lib = "\n".join(samples_names_not_in_library)
        exit_message += (
            "* The following samples are not found at courrent samples table:"
            f"\n{format_samples_not_in_lib}.\n"
        ) if len(samples_names_not_in_library) > 0 else ""

    # %% # FASTQs and PE
    # Check that there are 2 FASTQs when sample is PE and 1 when is not
    if (
        samples_table['PE'] &
            ((samples_table['fastq1'] == "") |
             (samples_table['fastq2'] == ""))).any():
        exit_message += (
            "* At least one sample is set as PE but only contains"
            " one FASTQ path. Interleaved FASTQs are not supported yet.\n")
        exit_script = True
        # Above code doesn't work without the extra parentheses on each
        # condition
    if (~samples_table['PE'] & ~samples_table['fastq2'].isna()).any():
        exit_message += (
            "* At least one sample is set as SR (PE == False) but"
            " contains a FASTQ path at column 'fastq2'. Please put it at "
            "column 'fastq1'\n")
        exit_script = True
    # %% FASTQ1 != FASTQ2
    if (samples_table['PE'] &
            (samples_table['fastq1'] == samples_table['fastq2'])).any():
        fastq_comp = samples_table['fastq1'] == samples_table['fastq2']
        bad_samples = samples_table.loc[fastq_comp]['sample_name'].tolist()
        exit_message += (
            "* The following samples have identical paths for fastq1 and 2:\n"
        )
        for i in bad_samples:
            exit_message += "\t" + i + "\n"
        exit_script = True

    # %% Dros_eq and dros_spike_in
    # Check that the samples with a "dros_equalization_group" have
    # "dros_spike_in" == T (I don't do the complementary becuase you might
    # have the sample to equalize on another library)
    if (~samples_table['dros_equalization_group'].isnull() &
            ~samples_table['dros_spike_in']).any():
        exit_message += (
            "* At least one of your samples has a 'dros_equalization_group' "
            "assigned to it but has the field 'dros_spike_in' set as False.\n"
        )
        exit_script = True

    # %% peak_ctrl_file_alias
    # Convert peak_ctrl_file_alias series to str (in case it's all NaN values)
    peak_ctrls = samples_table['peak_ctrl_file_alias'].astype(str)
    # Get list of wrongly filled samples
    bad_samples = samples_table.loc[
        peak_ctrls.str.contains("/|\\.| ", regex=True, na=False)
    ]
    # If necessary, assemble exit message
    if len(bad_samples) > 0:
        exit_message += (
            "* The following samples have incorrect values on column "
            "'peak_ctrl_file_alias:'\n"
            )
        for i in bad_samples.index:
            exit_message += "\t" + i + "\n"
        exit_script = True

    # %% Exit
    # Exit if any previous condition is met
    if exit_script:
        sys.exit(exit_message)

# %% Functions to get inputs/params


def align_fastq_input(w):
    """Get FASTQ paths for align_fastq rule

    wildcards
    ----------
    genomes_all : mm10|d6|hg19|hg38|mm10_f_d6|hg19_f_d6|hg38_f_d6|
                      mm10_x_CAST_EiJ_f_d6|mm10_x_CAST_EiJ
    sample : [^./ ]+

    Returns
    -------
    input_: dictionary
        Fastq/s path/s, taken from samples_table.csv
    """
    input_ = {"genome_path": config['genomes'][w.genomes_all]}
    seq_tech = ""

    if (samples_table.loc[w.sample, "library_technology"] == "adaptase"):
        seq_tech = ".adaptase_trimmed"

    if samples_table.loc[w.sample, "PE"]:
        fq1 = (f"Results/{w.sample}.adap_trimmed{seq_tech}.R1."
               "PE.fq.gz")
        fq2 = (f"Results/{w.sample}.adap_trimmed{seq_tech}.R2."
               "PE.fq.gz")
    else:
        fq1 = (f"Results/{w.sample}.adap_trimmed{seq_tech}.R1."
               "SE.fq.gz")
        fq2 = []

        # Cannot use "" here because snakemake will look for "" file
        # and give a "missing input" error. Instead, it returns
        # an empty list
    input_ |= {"fq1": fq1, "fq2": fq2}
    return input_


def call_peaks_macs2_input(w):
    """Get input for call_peaks_macs2 rule.

    Wildcards
    ----------
    peak_params : (bco_[0-9]+_)?qv_[0-9]+__[^/]*
        The first part (before the double underscore (__) of the wildcard are
            bco (if peak is broad) and qv parameters. The second part is the
            value of `peak_ctrl_file_alias` column for that sample.

    Returns
    -------
    input_ : dictionary
        Dictionary keys are "treat_bam", "treat_bai" and, unless
        peak_ctrl_file_alias == "no_input", "ctrl_bam" and "ctrl_bai"
    """
    bam = samples_table_2.loc[w.sample, 'dedup_flt_both_strds_bam']
    bai = bam + ".bai"
    input_ = {"treat_bam": bam, "treat_bai": bai}
    # Get ctrl bam and bai
    ctrl_alias = w.peak_params.split("__")  # "bco_1_qv_05__alias"
    ctrl_alias = ctrl_alias[1]  # "alias"
    if ctrl_alias != "no_input":
        ctrl_bam = config['MACS2']['control'][ctrl_alias]
        ctrl_bai = ctrl_bam + ".bai"
        input_ |= {
            "ctrl_bam": ctrl_bam,
            "ctrl_bai": ctrl_bai
        }
    return input_


def clip_1bp_input(w):
    """Get input for clip_1bp rule

    Wildcards
    ----------
    genomes_not_fused : "mm10|d6|hg19|hg38|mm10_x_CAST_EiJ"
    sample_with_extensions : "[^/]+"
    ss_PE: ((83-163|99-147).)?
    ss_SR: (inc|exc)_16

    Returns
    -------
    bam : string
        path of input bam file

    """
    if (w.ss_PE == ""):
        bam = (
            f"Results/{w.genomes_not_fused}/Bams/Single_strand/"
            f"Full_length_reads/{w.sample}.{w.genomes_final}.q_filt.srt."
            f"nodup.mit_filt.{w.ss_SR}.bam"
        )
    else:
        bam = (
            f"Results/{w.genomes_not_fused}/Bams/Single_strand/"
            f"Full_length_reads/{w.sample}.{w.genomes_final}.q_filt.srt."
            f"nodup.mit_filt.{w.ss_PE}bam"
        )
    return bam


def compute_matrix_outfiles_hs_input(w):
    """Get input for rule compute_matrix_outfiles_hs.

    Wildcards
    ----------
    sample : [^./ ]+
    hs_region : B6xCAST_(top_5000_pm_1000bp|
                         PRDM9_assymetric_hs_((invading|receiving)_strand|
                                              mm10_aligned)
                         )|
                top_5000_plus_minus_2000|asymetric_(watson|crick)_strong|
                x_non_par|autosomal_x_non_par_ctrl
    genomes_not_fused : mm10|d6|hg19|hg38|mm10_x_CAST_EiJ

    Returns
    -------
    dict
        'region', path to bed file with hotspots coordinates
        'bigwig', path to bigwig file

    """
    hotspot_protein = "spo11"
    if samples_table.loc[w.sample, 'B6xCAST']:
        hotspot_protein = "dmc1"

    hotspots = config['references']['mm10'][hotspot_protein][w.hs_region]
    bigwig = (f"Results/{w.genomes_not_fused}/Bigwigs/Coverage/{w.strands}/"
              f"{w.cov_params}/{w.sample}.{w.genomes_final}{w.extension}"
              f"{w.strand}.bw"
              )

    return {"region": hotspots,
            "bigwig": bigwig}


def dros_normalization_input(w):
    r"""Get input for dros_normalization rule

    Wildcards
    ----------
    strand : (\.(83-163|99-147|inc_16|exc_16))?
    sample : [^./ ]+
    dros_eq_group : [^./ ]+

    Returns
    -------
    Dictionary
        'report', path to a csv file containing the drosophila normalization
            factors.
        'bam', bam file to normalize
        'bai', index file of bam to normalize
    """
    if w.strand == "":
        bam = samples_table_2.loc[w.sample, 'dedup_flt_both_strds_bam']
    elif ("83-163" in w.strand) | ("inc_16" in w.strand):
        bam = samples_table_2.loc[w.sample, 'ss_83_or_i16_bam']
    elif ("99-147" in w.strand) | ("exc_16" in w.strand):
        bam = samples_table_2.loc[w.sample, 'ss_99_or_e16_bam']

    bai = bam + ".bai"
    return {
        "report": ("Results/d6/Analysis/drosophila_normalization/"
                   f"{w.dros_eq_group}/drosophila_equalization_report.tsv"),
        "bam": bam,
        "bai": bai
    }


def dros_normalization_report_input(w):
    """Get input for rule dros_normalization_report.

    Wildcards
    ----------
    dros_eq_group : [^./ ]+

    Returns
    -------
    List
        Paths to outputs of samtools_flagstat of samples to be normalized on a
            given 'dros_equalization_group'.
    """
    df = samples_table_2.loc[
        samples_table_2['dros_equalization_group'].str.match(
            w.dros_eq_group, na=False
        ), :]
    return df['processed_flagstat_dros'].tolist()


def filter_bam_params(w):
    r"""Assemble part of the bash commands as params in rule filter_bam

    Wildcards
    ----------
    sample: [^./ ]+
    genomes_not_fused: "mm10|d6|hg19|hg38|mm10_x_CAST_EiJ"
    genomes_final:(
        "(?<=\.)(?P<interest_genome>mm10|hg19|hg38|mm10_x_CAST_EiJ)+"
        "(_f_d6\.((?P=interest_genome)|d6))?(?=\.)")

    Returns
    -------
    string
        Filtering options for samtools view (-F and -f), according to wether
        the sample is SE or PE.
    """
    if samples_table.loc[w.sample, "PE"]:
        filter_ = "-F 3852 -f 3"
    else:
        filter_ = "-F 3844"
    return filter_


def filter_peaks_blk_grey_list_input(w):
    """Get input for rule filter_peaks_blk_grey_list

    Output of rule:
    ----------
    ("Results/{genomes_not_fused}/Peaks/MACS2/{peak_type}/"
         "{peak_params}/Black-grey_filtered/{sample}."
         "{genomes_final}{extension}.{peak_type}Peak"
    )

    Widcards
    ----------
    extension = r"(\.q_filt\.srt\.nodup\.mit_filt)?"

    Returns
    -------
    dict:
        "peaks": list with paths of all peaks from that reference genome
        "blck_gry_lst": path of black-greylist bed file if available, if not
            path of black-list

    """
    peaks = (f"Results/{w.genomes_not_fused}/Peaks/MACS2/{w.peak_type}/"
             f"{w.peak_params}/{w.sample}.{w.genomes_final}{w.extension}_"
             f"peaks.{w.peak_type}Peak")

    blck_gry_lst = config['references'][w.genomes_not_fused]['black_grey']
    if (config['references'][w.genomes_not_fused]['black_grey'] == ""):
        blck_gry_lst = config['references'][w.genomes_not_fused]['blacklist']

    return {"peaks": peaks,
            "blck_gry_lst": blck_gry_lst}


def FRIP_input(w):
    r"""Get input for FRIP rule.

    Wildcards
    ----------
    sample : [^./ ]+
    genomes_not_fused : mm10|d6|hg19|hg38|mm10_x_CAST_EiJ
    genomes_final : (?<=\.)(?P<interest_genome>mm10|hg19|hg38|mm10_x_CAST_EiJ)+
                        (_f_d6\.((?P=interest_genome)|d6))?(?=\.)"
    extension : (\.q_filt\.srt\.nodup\.mit_filt)?
    peak_type : narrow|broad
    peak_params : (bco_[0-9]+_)?qv_[0-9]+__[^/]*
    extension : (\.q_filt\.srt\.nodup\.mit_filt)?

    Returns
    -------
    dict
        'peak', path to peaks bed file
        'bam', path to bam file
        'bai', path to bai file
    """
    sample = samples_table_2.loc[w.sample, "dedup_flt_both_strds_bam"]
    bam = (f"Results/{w.genomes_not_fused}/Bams/Both_strands/" +
           f"{w.sample}.{w.genomes_final}{w.extension}.bam")

    return {"peak": (f"Results/{w.genomes_not_fused}/Peaks/MACS2/{w.peak_type}"
                     f"/{w.peak_params}/Black-grey_filtered/{w.sample}."
                     f"{w.genomes_final}{w.extension}.{w.peak_type}Peak"),
            "bam": bam,
            "bai": bam + ".bai"
            }


def get_rv_fw_strand_input(w):
    """Get input for rule get_rv_fw_strand

    Wildcards
    ----------
    sample_with_extensions: [^/]+

    Returns
    -------
    str
        Input's path

    # This function and the rule it affects need to be checked
    """
    if re.match(".*(83-163|99-147)", w.sample_with_extensions):
        strand = "Single_strand/Full_length_reads"
    else:
        strand = "Both_strands"
    return (f"Results/{w.genomes_not_fused}/Bams/{strand}/"
            f"{w.sample_with_extensions}.bam")


def get_call_peaks_macs2_params(w):
    """Get parameters for rule call_peaks_macs2.

    Wildcards
    ----------
    genomes_not_fused : mm10|d6|hg19|hg38|mm10_x_CAST_EiJ
        Genomes not fused to spike-in genome
    peak_params : (bco_[0-9]+_)?qv_[0-9]+__[^/]*
        The first part (before the double underscore (__) of the wildcard are
            bco (if peak is broad) and qv parameters. The second part is the
            value of `peak_ctrl_file_alias` column for that sample.

    Returns
    -------
    dict
        Contains strings that will be used to assemble the shell script.
        More specifically: genome, PE, extension, peak_type_options, qv_bco
        and ctrl
    """
    if re.match(".*mm.*", w.genomes_not_fused):
        genome = "mm"
    elif re.match(".*hg.*", w.genomes_not_fused):
        genome = "hs"
    PE = ""
    ext = ""
    parameters = w.peak_params.split("__")
    bco_qv = parameters[0].split("_")
    ctrl_alias = parameters[-1]
    ctrl_bam = config['MACS2']['control'][ctrl_alias]
    if parameters[-1] == "no_input":
        ctrl = ""
    else:
        ctrl = f"-c {ctrl_bam}"
    if samples_table.loc[w.sample, "PE"]:
        PE = "--format BAMPE"
    else:
        ext = f"--extsize {config['MACS2']['extension']}"
    if w.peak_type == "narrow":
        peak_type = "--call-summits"
        q = f"-q 0.{bco_qv[1]}"
    else:
        peak_type = "--broad"
        q = f"-q 0.{bco_qv[3]} --broad-cutoff 0.{bco_qv[1]}"
    return {
        "genome": genome,
        "PE": PE,
        "extension": ext,
        "peak_type_options": peak_type,
        "qv_bco": q,
        "ctrl": ctrl}


def get_strand_sep_bams_params(w):
    """Get parameters for rule get_strand_sep_bams

    Wildcads
    ----------
    flag : "(83|163|99|147|inc_16|exc_16)"

    Returns
    -------
    dict : string
        (-f 83|99|163|99|16) | (-F 16)
    """
    filter_ = "-F " if "exc" in f"{w.flag}" else "-f "
    flag = w.flag.removeprefix("inc_").removeprefix("exc_")
    argument = filter_ + flag

    return argument


def intersect_peaks_HSs_list_input(w):
    r"""Get input for intersect_peaks_HSs_list rule.

    Parameters
    ----------
    sample : [^./ ]+
    genomes_not_fused : mm10|d6|hg19|hg38|mm10_x_CAST_EiJ
    peak_type: narrow|broad
    peak_params: (bco_[0-9]+_)?qv_[0-9]+__[^/]*
    genomes_final : (
        "(?<=\.)(?P<interest_genome>mm10|hg19|hg38|mm10_x_CAST_EiJ)+"
        "(_f_d6\.((?P=interest_genome)|d6))?(?=\.)")
    extension : (\.q_filt\.srt\.nodup\.mit_filt)?

    Returns
    -------
    dict
        'hotspots', path to hotspots (from either mm10 or B6xCAST) bed file
        'sample_peaks', path to peaks bed file
    """
    if samples_table.loc[w.sample, 'B6xCAST']:
        hotspots = config[
            'references']['mm10']['dmc1']['B6xCAST_top_5000_pm_1000bp']
    else:
        hotspots = config[
            'references']['mm10']['spo11']['top_5000_plus_minus_2000']
    return {
        "hotspots": hotspots,
        "sample_peaks": (
            f"Results/{w.genomes_not_fused}/Peaks/MACS2/"
            f"{w.peak_type}/{w.peak_params}/Black-grey_filtered/"
            f"{w.sample}.{w.genomes_final}{w.extension}.{w.peak_type}Peak"
        )
    }


def markdown_report_aggregate_profiles_input(w):
    """Get input for rule markdown_report_aggregate_profiles.

    Wildcards
    ----------
    genomes_not_fused: mm10|d6|hg19|hg38|mm10_x_CAST_EiJ


    Returns
    -------
    Dictionary
        'peaks_summary': string or list. Path to table with summary of peaks
            for all samples with the same reference genome. If no peaks were
            asked for, it delivers an empty list.
        'bamfiles_reads': list. Paths to samtools_flagstat outputs (for
            samples with the same reference genome as the markdown report).
            These files will be used in the report to get the number of reads
            of the sample in the filtered bam file.
        'mm10_top5000_asmtric_auto_XnonPAR_ag_profs'|'mm10_top_5000_ag_profs':
            string. Only for mm10. Path to '.RData' object containing the
            aggregate profiles in either "mm10 top 5000 hotspots", or in those
            same hotspots plus XnonPAR, autosomal and assymetric (left vs
            right of DSB) hotspots. When all hotsopots lists are asked for
            (later case), the rule that provides the R object is the
            "wrangle_X_nonPAR_asymmetric_HS" rule, as opposed to
            "process_aggregate_profiles" when it is only the top 5000.
        'B6xCAST_top_5000_pm_1000bp': string. Only for mm10. Path to ".RData"
            object containing the aggreagte profiles of all samples aligned to
            B6xCAST fused genome, in top 5000 B6xCAST hotspots.
        'B6xCAST_PRDM9_assymetric_hs_invading_strand': string. Path to ".RData"
            object containing the aggreagte profiles of all samples aligned to
            B6xCAST fused genome, in B6xCAST hotspots that bind PRDM9
            asymmetrically, on the invading strand.
        'B6xCAST_PRDM9_assymetric_hs_receiving_strand': string. same as
            previous but receiving/template strand.
        'B6xCAST_PRDM9_assymetric_hs_mm10_aligned': string. path to ".RData"
            object containing the aggreagte profiles of all samples from
            B6xCAST mice, aligned to mm10 genome (as opossed to B6xCAST fused
            genome).
    """
    # Get peak_summary path
    ref_genome_peaks = samples_table.loc[
        samples_table["reference_genome"] == w.genomes_not_fused,
        'peak_ctrl_file_alias'
    ]
    if pd.notna(ref_genome_peaks).any():
        peaks_summary = {
            "peaks_summary": (f"Results/{w.genomes_not_fused}/Analysis/"
                              "Peaks_summary.tsv")
        }
    else:
        peaks_summary = {
            "peaks_summary": []
        }

    # Get processed_flagstat paths
    selection_criteria = (
        samples_table_2['reference_genome'] == w.genomes_not_fused
    )
    reads = samples_table_2.loc[selection_criteria, 'processed_flagstat']
    reads = reads.values.tolist()
    reads = [x for x in reads if not re.search("_MERGED\\.", x)]
    bamfiles_reads = {"bamfiles_reads": reads}

    # Get fastp reports (fastq # of reads)
    selection_criteria_pe = (
        (samples_table_2['reference_genome'] == w.genomes_not_fused) &
        samples_table_2['PE']
    )
    samples_fastp_pe = samples_table_2.loc[selection_criteria_pe,
                                           'sample_name']
    samples_fastp_pe = samples_fastp_pe.values.tolist()
    samples_fastp_pe = [x for x in samples_fastp_pe if not
                        re.search("_MERGED$", x)]
    samples_fastp_pe = [f"Results/FASTQ_reports/{sample}.PE.fastp.json" for
                        sample in samples_fastp_pe]

    selection_criteria_se = (
        (samples_table_2['reference_genome'] == w.genomes_not_fused) &
        ~samples_table_2['PE']
    )
    samples_fastp_se = samples_table_2.loc[selection_criteria_se, 
                                           'sample_name']
    samples_fastp_se = samples_fastp_se.values.tolist()
    samples_fastp_se = [x for x in samples_fastp_se if not
                        re.search("_MERGED$", x)]
    samples_fastp_se = [f"Results/FASTQ_reports/{sample}.SE.fastp.json" for
                        sample in samples_fastp_se]

    samples_fastp = {"fastp": samples_fastp_pe + samples_fastp_se}

    mm10_agg_profiles = {}
    B6xCAST_agg_profiles = {}
    if (samples_table['top5000_HS_heatmap'] &
        samples_table['get_single_strand'] &
            ~samples_table['B6xCAST']).any():
        # 'Lib_name_aggregate_profiles_data.RData' is the result of processing
        # all top5000, asymmetrix and XnonPAR data (output of
        # 'wrangle_X_nonPAR_asymmetric_HS' rule).
        mm10_agg_profiles = {
            "mm10_top5000_asmtric_auto_XnonPAR_ag_profs": (
                f"Results/{w.genomes_not_fused}/Analysis/"
                "Heatmaps_and_aggregate_profiles/Hotspots/"
                f"{config['library']['name']}_aggregate_profiles_data.RData"
            )
        }
    elif ((samples_table['top5000_HS_heatmap'] &
           ~samples_table['B6xCAST']).any()):
        mm10_agg_profiles = {
            "mm10_top_5000_ag_profs": (
                f"Results/{w.genomes_not_fused}/Analysis"
                "/Heatmaps_and_aggregate_profiles/Hotspots/"
                f"{config['library']['name']}_top_5000_plus_minus_2000.RData"
            )
        }

    if ((samples_table['top5000_HS_heatmap'] &
         samples_table['B6xCAST']).any()):
        B6xCAST_agg_profiles = {
            "B6xCAST_top_5000_pm_1000bp": (
                f"Results/{w.genomes_not_fused}/Analysis"
                "/Heatmaps_and_aggregate_profiles/Hotspots/"
                f"{config['library']['name']}_B6xCAST_top_5000_pm_1000bp.RData"
            ),
            "B6xCAST_PRDM9_assymetric_hs_invading_strand": (
                "Results/mm10_x_CAST_EiJ/Analysis/"
                "Heatmaps_and_aggregate_profiles/Hotspots/"
                f"{config['library']['name']}_B6xCAST_PRDM9_assymetric_hs_"
                "invading_strand.RData"
            ),
            "B6xCAST_PRDM9_assymetric_hs_receiving_strand": (
                "Results/mm10_x_CAST_EiJ/Analysis/"
                "Heatmaps_and_aggregate_profiles/Hotspots/"
                f"{config['library']['name']}_B6xCAST_PRDM9_assymetric_hs_"
                "receiving_strand.RData"
            ),
            "B6xCAST_PRDM9_assymetric_hs_mm10_aligned": (
                "Results/mm10/Analysis/Heatmaps_and_aggregate_profiles/"
                f"Hotspots/{config['library']['name']}_"
                "B6xCAST_PRDM9_assymetric_hs_mm10_aligned.RData"),
        }
    # Return
    if w.genomes_not_fused == "mm10":
        return (peaks_summary | bamfiles_reads | samples_fastp |
                mm10_agg_profiles | B6xCAST_agg_profiles)
    else:
        return peaks_summary | bamfiles_reads | samples_fastp


def merge_bams_input(w):
    r"""Get input for merge_bams rule.
    Information about what to merge is on the 'merge_with' column from
    samples table, and it should be a space-separated list of sample
    names found the same samples table, or absolute paths to files (or both).

    Wildcards
    ----------
    sample : [^./ ]+
    genomes_all : mm10|d6|hg19|hg38|mm10_f_d6|hg19_f_d6|hg38_f_d6|
                  mm10_x_CAST_EiJ_f_d6|mm10_x_CAST_EiJ


    Returns
    -------
    samples : List
        List with paths to the raw bam files to be merged.
    """
    samples_merge = samples_table.loc[w.sample, 'merge_with'].split()
    samples = [f"Results/{w.sample}.{w.genomes_all}.bam"]
    # If sample to merge with is a path, use it as is, otherwise look for the
    # raw bam file on samples_table_2
    for sample in samples_merge:
        if re.search("/", sample):
            samples += [sample]
        else:
            second_bam = f"Results/{sample}.{w.genomes_all}.bam"
            samples += [second_bam]
    return samples


def process_aggregate_profiles_inputs(w):
    """Get inputs for rule process_aggregate_profiles.

    Wildcards
    ----------
    genomes_not_fused : mm10|d6|hg19|hg38|mm10_x_CAST_EiJ
    hs_region : B6xCAST_(top_5000_pm_1000bp|
                         PRDM9_assymetric_hs_((invading|receiving)_strand|
                                              mm10_aligned))|
                top_5000_plus_minus_2000|asymetric_(watson|crick)_strong|
                x_non_par|autosomal_x_non_par_ctrl

    Returns
    -------
    matrixes_list : List
        Contains all matrixes' paths for all the
        sample/strand combinations possible for one particular hostpots list
        (top_5000_plus_minus_2000, x_non_par, etc.).

    Matrixes paths are taken from samples_table_2. So, the actual decition of
    which sample will have a matrix in which hotspot list and in which strands
    is actually taking place during samples_table_2 generation
    (suffixes.generate_samples_table_2())
    """

    matrixes_df = samples_table_2.filter(
        regex=f"^{w.hs_region}.*matrix$"
    )
    matrixes_array = matrixes_df.to_numpy().ravel()
    matrixes_list = matrixes_array[~pd.isnull(matrixes_array)].tolist()
    return matrixes_list


def sumarize_peak_count_input(w):
    """Get input for rule sumarize_peak_count.

    Wildcards
    ----------
    genomes_not_fused: mm10|d6|hg19|hg38|mm10_x_CAST_EiJ

    Returns
    -------
    peaks : dictionary
        'narrow_all': contains paths for all bed files with black/grey-list-
                    filtered narrow peaks (for that particular reference
                    genome),
        'broad_all': same but broad peaks,
        'narrow_hs': bed files with black/grey-list-filtered peaks intersected
                    with hs list (this files only exist for mm10),
        'broad_hs': same as above but broad peaks
    """
    genome_filtered = samples_table_2[
        samples_table_2['reference_genome'] == w.genomes_not_fused
    ]

    peak_types = ["narrow", "broad"]
    peaks = {}
    for peak_type in peak_types:
        selection_criteria_all = (
            genome_filtered[f'{peak_type}_peak_bl_gr_flt'].notnull()
        )
        selection_criteria_hs = (
            genome_filtered[f'{peak_type}_peak_bl_gr_flt_hs_int'].notnull()
        )

        peaks |= {
            f"{peak_type}_all": genome_filtered.loc[
                selection_criteria_all,
                f'{peak_type}_peak_bl_gr_flt'
                ].values.tolist(),
            f"{peak_type}_hs": genome_filtered.loc[
                selection_criteria_hs,
                f'{peak_type}_peak_bl_gr_flt_hs_int'
                ].values.tolist(),
        }
    if not (len(peaks['narrow_all']) > 1):
        exit
    return peaks


def samstats_samtools_flagstat_input(w):
    """Get input for rules samstats and samtools_flagstat.

    Wildcards
    ----------
    sample : [^./ ]+
    genomes_not_fused : mm10|d6|hg19|hg38|mm10_x_CAST_EiJ
    bam_type : Raw_bam|Processed_bam

    Returns
    -------
    dictionary
        Bam and bai paths.
    """
    genome = samples_table.loc[w.sample, 'reference_genome']
    if f"{w.bam_type}" == 'Raw_bam':
        if w.genomes_not_fused != "d6":
            bam = samples_table_2.loc[w.sample, 'raw_bam']
        else:
            bam = f"Results/{w.w.sample}.{genome}_f_d6.d6.bam"
    else:
        if w.genomes_not_fused != "d6":
            bam = samples_table_2.loc[w.sample, 'dedup_flt_both_strds_bam']
        else:
            bam = (f"Results/d6/Bams/Both_strands/{w.sample}.{genome}_f_d6.d6."
                   "q_filt.srt.nodup.mit_filt.bam")
    bai = f"{bam}.bai"
    return {
        "bam": bam,
        "bai": bai
    }


def trim_adapters_PE_input(w):
    """Get FASTQ paths for trim_adapters_PE rule

    wildcards
    ----------
    sample : [^./ ]+

    Returns
    -------
    input_: dictionary
        Fastq/s path/s, taken from samples_table.csv
    """
    return {"fastq1": samples_table.loc[w.sample, "fastq1"],
            "fastq2": samples_table.loc[w.sample, "fastq2"]}
