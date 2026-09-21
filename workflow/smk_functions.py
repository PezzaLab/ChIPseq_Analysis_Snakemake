# TODO: add module docstring

import pandas as pd
import re
import sys
import os

samples_table = None
samples_table_2 = None
config = None
samples_table_no_merged_samples = None
include_hotspots = True

def expand_sample_table(raw_df):
    r"""Expands sample table rows where 'reference_genome' contains multiple genomes.

    Supports:
    - Multiple genomes specified via comma or semicolon (e.g. 'mm10;mm39' or 'mm10, mm39').
    - Broadcasting: any genome-dependent column with a single value is broadcast
      to all genomes specified for that sample.
    - Sets of values: any column with multiple values (delimited by ';' or ',')
      must match the number of reference genomes exactly (1-to-1 mapping).
    - 'merge_with' uses ';' as delimiter between genomes to allow space-separated
      BAM/sample lists for each genome.
    - Safe boolean conversion for 'True'/'False'/'T'/'F'/'1'/'0'.
    - Assigns unique index 'f"{sample_name}.{reference_genome}"'.
    """
    import numpy as np
    expanded_rows = []

    sample_intrinsic_cols = {"sample_name", "fastq1", "fastq2", "library_technology", "PE"}
    boolean_cols = {
        "PE", "dros_spike_in", "get_single_strand",
        "Clip_reads_to_1bp_on_5_prime", "top5000_HS_heatmap",
        "Size_DNA_top_5000_HS", "B6xCAST"
    }
    bool_map = {
        "true": True, "t": True, "1": True,
        "false": False, "f": False, "0": False
    }

    for idx, row in raw_df.iterrows():
        sample_name = str(row["sample_name"]).strip()
        ref_val = row.get("reference_genome", "")
        if pd.isna(ref_val) or not str(ref_val).strip():
            sys.exit(f"\n\n* Sample '{sample_name}' has no reference_genome defined.\n")

        # Split reference_genome by comma or semicolon
        genomes = [g.strip() for g in re.split(r"[,;]", str(ref_val)) if g.strip()]
        num_genomes = len(genomes)

        col_values = {}
        for col in raw_df.columns:
            val = row[col]

            if col == "sample_name":
                col_values[col] = [sample_name] * num_genomes
            elif col in {"fastq1", "fastq2"}:
                val_clean = np.nan if (pd.isna(val) or str(val).strip() in ["-", ""]) else str(val).strip()
                col_values[col] = [val_clean] * num_genomes
            elif col == "reference_genome":
                col_values[col] = genomes
            elif col == "merge_with":
                if pd.isna(val) or str(val).strip() in ["-", ""]:
                    col_values[col] = [np.nan] * num_genomes
                else:
                    parts = [p.strip() for p in str(val).split(";")]
                    if len(parts) == 1:
                        p_val = np.nan if parts[0] in ["-", ""] else parts[0]
                        col_values[col] = [p_val] * num_genomes
                    elif len(parts) == num_genomes:
                        col_values[col] = [np.nan if p in ["-", ""] else p for p in parts]
                    else:
                        sys.exit(
                            f"\n\n* Sample '{sample_name}': column 'merge_with' has {len(parts)} values (semicolon-separated), "
                            f"but {num_genomes} reference genomes were specified ({genomes}).\n"
                        )
            elif col in boolean_cols:
                if isinstance(val, (bool, np.bool_)):
                    col_values[col] = [bool(val)] * num_genomes
                elif pd.isna(val) or str(val).strip() in ["-", ""]:
                    col_values[col] = [False] * num_genomes
                else:
                    parts = [p.strip() for p in re.split(r"[,;]", str(val)) if p.strip()]
                    converted = []
                    for p in parts:
                        p_low = p.lower()
                        if p_low in bool_map:
                            converted.append(bool_map[p_low])
                        else:
                            sys.exit(
                                f"\n\n* Sample '{sample_name}': invalid boolean value '{p}' in column '{col}'.\n"
                            )
                    if len(converted) == 1:
                        col_values[col] = converted * num_genomes
                    elif len(converted) == num_genomes:
                        col_values[col] = converted
                    else:
                        sys.exit(
                            f"\n\n* Sample '{sample_name}': column '{col}' has {len(converted)} values, "
                            f"but {num_genomes} reference genomes were specified ({genomes}).\n"
                        )
            else:
                # Other string / metadata columns
                if pd.isna(val) or str(val).strip() in ["-", ""]:
                    col_values[col] = [np.nan] * num_genomes
                else:
                    parts = [p.strip() for p in re.split(r"[,;]", str(val)) if p.strip()]
                    if len(parts) == 1:
                        p_val = np.nan if parts[0] == "-" else parts[0]
                        col_values[col] = [p_val] * num_genomes
                    elif len(parts) == num_genomes:
                        col_values[col] = [np.nan if p == "-" else p for p in parts]
                    else:
                        sys.exit(
                            f"\n\n* Sample '{sample_name}': column '{col}' has {len(parts)} values, "
                            f"but {num_genomes} reference genomes were specified ({genomes}).\n"
                        )

        # Check sample-intrinsic consistency
        for col in sample_intrinsic_cols:
            if col in col_values and len(set(col_values[col])) > 1:
                sys.exit(
                    f"\n\n* Sample '{sample_name}': sample-intrinsic column '{col}' has conflicting values "
                    f"{col_values[col]} across genomes for the same sample.\n"
                )

        for i, g in enumerate(genomes):
            rec = {col: col_values[col][i] for col in raw_df.columns}
            rec["sample_entry_id"] = f"{sample_name}.{g}"
            expanded_rows.append(rec)

    df_expanded = pd.DataFrame(expanded_rows)
    df_expanded.set_index("sample_entry_id", drop=False, inplace=True, verify_integrity=True)
    return df_expanded


def check_sample_table_format(samples_table):
    r"""Check if samples_table.csv has been properly filled.

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
        * Booleans columns should have boolean values only
        * 'peak_ctrl_file_alias' should not contain '/', '.' or ' '
        * If 'Clip_reads_to_1bp_on_5_prime' is TRUE, 'PE' must be FALSE (must be SE).

    Parameters:
    -----------
        samples_table : dataframe produced by reading and expanding samples_table.csv
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
        col for col in cols_to_check
        if col in samples_table.columns and samples_table[col].dtype != "bool"
    ]

    if (len(booleans_check) > 0):
        exit_bool_message = (
            "* The following column/s has/ve at least one"
            " row filled with something different than 'T', 'F', 'True',"
            " 'False', 'TRUE' or 'FALSE':\n"
            f"\t({booleans_check})"
        )
        sys.exit(exit_bool_message)

    # %% Samples names check
    # Check names of samples don't contain / or . or finish in "_MERGED"
    sample_names = samples_table['sample_name']
    offending_names = sample_names[sample_names.str.contains(r'[./ ]')].to_list()
    offending_names += sample_names[sample_names.str.contains('_MERGED$')].tolist()

    if len(offending_names) > 0:
        offending_names_str = "\n".join([f"- {s}" for s in set(offending_names)])
        exit_message += ("\n\n* The following sample names contain"
                         " a dot ('.'), a slash ('/'), a space (' ') or ends"
                         " with the word '_MERGED'. These are not allowed in"
                         " sample names. Modify names and try again\n"
                         f"{offending_names_str}"
                         )
        exit_script = True

    # %% Check that (sample_name, reference_genome) pairs are not repeated
    sample_genome_pairs = samples_table['sample_name'] + " : " + samples_table['reference_genome']
    duplicated_pairs = sample_genome_pairs[sample_genome_pairs.duplicated()].unique().tolist()
    if len(duplicated_pairs) > 0:
        dup_str = "\n".join([f"\t- {p}" for p in duplicated_pairs])
        exit_message += (
            "\n\n* The following sample_name : reference_genome combinations are "
            f"repeated in your samples table:\n{dup_str}\n"
        )
        exit_script = True

    # %% Merged samples
    if pd.notnull(samples_table['merge_with']).any():
        # Check that "merge_with" reference is not a deduplicated/filtered
        # file for spiked-in samples (it has to be raw bam)
        mask_merge = (
            samples_table['merge_with'].astype(str).str.contains(r"\.nodup\.") &
            samples_table['dros_spike_in']
        )
        offending_merge_df = samples_table[mask_merge]

        if mask_merge.any():
            offending_merge = "\n".join(
                [f"\t- {s}" for s in offending_merge_df["sample_name"].unique()]
            )
            exit_message += (
                "\n\n* In the following samples the 'merge_with' parameter is "
                "a deduplicated/filtered bam file AND has drosophila spike in."
                " When sample is spiked, the 'merge_with' file has to be the "
                "raw bam, otherwise you'll loose the drosophila reads.\n"
                f"{offending_merge}"
            )
            exit_script = True

        # Check that merge samples are either an absolute path or another
        # sample from samples_table
        samples_names_not_in_library = []
        samples_merging_themselves = []
        for idx, row in samples_table.iterrows():
            sample = row['sample_name']
            if pd.notnull(row['merge_with']):
                merge_samples = str(row['merge_with']).split()
                for merge_sample in merge_samples:
                    absolute_path = merge_sample.startswith("/")
                    in_library = merge_sample in samples_table['sample_name'].values
                    merge_itself = merge_sample == sample
                    if not absolute_path and not in_library:
                        samples_names_not_in_library.append(f"\t{merge_sample}")
                        exit_script = True
                    if merge_itself:
                        samples_merging_themselves.append(f"\t{merge_sample}")
                        exit_script = True

        if len(samples_names_not_in_library) > 0:
            format_samples_not_in_lib = "\n".join(set(samples_names_not_in_library))
            exit_message += (
                "\n\n* The following samples are not found at current samples "
                f"table:\n{format_samples_not_in_lib}\n"
            )
        if len(samples_merging_themselves) > 0:
            format_samples_merging_themselves = "\n".join(set(samples_merging_themselves))
            exit_message += (
                "\n\n* The following samples are merging with themselves "
                f"\n{format_samples_merging_themselves}\n"
            )

    # %% FASTQs and PE
    mask_PE = (
        samples_table['PE'] & (
            samples_table['fastq1'].isna() | (samples_table['fastq1'] == "") |
            samples_table['fastq2'].isna() | (samples_table['fastq2'] == "")
        )
    )
    offending_PE = samples_table[mask_PE]["sample_name"].unique().tolist()
    if len(offending_PE) > 0:
        offending_PE_strg = "\n".join([f"\t- {s}" for s in offending_PE])
        exit_message += (
            "\n\n* The following samples are set as PE but only contain"
            " one FASTQ path. Interleaved FASTQs are not supported yet.\n"
            f"{offending_PE_strg}"
        )
        exit_script = True

    mask_SR = ~samples_table['PE'] & ~samples_table['fastq2'].isna()
    offending_SR = samples_table[mask_SR]["sample_name"].unique().tolist()
    if len(offending_SR) > 0:
        offending_SR_strg = "\n".join([f"\t- {s}" for s in offending_SR])
        exit_message += (
            "\n\n* The following samples are set as SR (PE == False) but"
            " contain a FASTQ path at column 'fastq2'. Please check and fix\n"
            f"{offending_SR_strg}"
        )
        exit_script = True

    # %% FASTQ1 != FASTQ2
    mask_FQs = (
        samples_table['PE'] &
        (samples_table['fastq1'] == samples_table['fastq2'])
    )
    offending_FQs = samples_table[mask_FQs]["sample_name"].unique().tolist()
    if len(offending_FQs) > 0:
        offending_FQs_strg = "\n".join([f"\t- {s}" for s in offending_FQs])
        exit_message += (
            "\n\n* The following samples have identical paths for "
            "fastq1 and 2:\n"
            f"{offending_FQs_strg}"
        )
        exit_script = True

    # %% peak_ctrl_file_alias
    peak_ctrls = samples_table['peak_ctrl_file_alias'].astype(str)
    bad_samples = samples_table.loc[
        peak_ctrls.str.contains(r"/|\.| ", regex=True, na=False)
    ]
    if len(bad_samples) > 0:
        exit_message += (
            "\n\n* The following samples have incorrect values on column "
            "'peak_ctrl_file_alias:'\n"
        )
        for i in bad_samples['sample_name'].unique():
            exit_message += "\t" + i + "\n"
        exit_script = True

    # %% Clip_reads and get_single_strand
    clip_Y_get_ss_N = [
        row['sample_name'] for idx, row in samples_table.iterrows()
        if row["Clip_reads_to_1bp_on_5_prime"] and not row["get_single_strand"]
    ]
    if len(clip_Y_get_ss_N) > 0:
        exit_message += (
            "\n\n* The following samples are set TRUE for "
            "`Clip_reads_to_1bp_on_5_prime` but FALSE for `get_single_strand`\n"
        )
        for i in set(clip_Y_get_ss_N):
            exit_message += "\t" + i + "\n"
        exit_message += ("Either set `Clip_reads_to_1bp_on_5_prime` to FALSE "
                         "or `get_single_strand` to TRUE\n"
                         )
        exit_script = True

    # %% Clip_reads and PE
    clip_Y_PE_Y = [
        row['sample_name'] for idx, row in samples_table.iterrows()
        if row["Clip_reads_to_1bp_on_5_prime"] and row["PE"]
    ]
    if len(clip_Y_PE_Y) > 0:
        exit_message += (
            "\n\n* The following samples are set TRUE for "
            "`Clip_reads_to_1bp_on_5_prime` but are also set TRUE for `PE`:\n"
        )
        for i in set(clip_Y_PE_Y):
            exit_message += "\t" + i + "\n"
        exit_message += (
            "Samples with `Clip_reads_to_1bp_on_5_prime` must be Single-End (PE=False). "
            "Either set `PE` to FALSE (and provide only fastq1) or set "
            "`Clip_reads_to_1bp_on_5_prime` to FALSE.\n"
        )
        exit_script = True

    # %% library_technology
    allowed_values = {'adaptase', 'regular'}
    mask = ~samples_table['library_technology'].isin(allowed_values)
    non_matching_rows = samples_table.loc[mask, ['sample_name', 'library_technology']]

    if not non_matching_rows.empty:
        exit_message += (
            f"\n\n* 'library_technology' can only be one of: {allowed_values}.\n"
            "The following samples have other values:\n" +
            non_matching_rows.drop_duplicates().to_string()
        )
        exit_script = True

    # %% Exit
    if exit_script:
        sys.exit(exit_message)

# %% Helper functions for safe sample and genome lookups

def get_sample_fastq1(sample):
    r"""Get FASTQ1 path for a sample name."""
    rows = samples_table[samples_table['sample_name'] == sample]
    if rows.empty:
        sys.exit(f"\n\n* Sample '{sample}' not found in samples_table\n")
    return rows['fastq1'].iloc[0]


def get_sample_fastq2(sample):
    r"""Get FASTQ2 path for a sample name."""
    rows = samples_table[samples_table['sample_name'] == sample]
    if rows.empty:
        sys.exit(f"\n\n* Sample '{sample}' not found in samples_table\n")
    return rows['fastq2'].iloc[0]


def get_processed_bam(sample, genome):
    r"""Get deduplicated, filtered BAM for a given sample and reference genome."""
    rows = samples_table_2[
        (samples_table_2['sample_name'] == sample) &
        (samples_table_2['reference_genome'] == genome)
    ]
    if rows.empty:
        sys.exit(f"\n\n* Sample '{sample}' for genome '{genome}' not found in samples_table_2\n")
    return rows['dedup_flt_both_strds_bam'].iloc[0]


def get_raw_bam(sample, genome):
    r"""Get raw BAM for a given sample and reference genome."""
    rows = samples_table_2[
        (samples_table_2['sample_name'] == sample) &
        (samples_table_2['reference_genome'] == genome)
    ]
    if rows.empty:
        sys.exit(f"\n\n* Sample '{sample}' for genome '{genome}' not found in samples_table_2\n")
    return rows['raw_bam'].iloc[0]


# %% Functions to get inputs/params


def align_fastq_input(w):
    r"""Get FASTQ paths for align_fastq rule

    wildcards
    ----------
    genomes_all : mm10|d6|hg19|hg38|mm10_f_d6|hg19_f_d6|hg38_f_d6|
                      mm10_x_CAST_EiJ_f_d6|mm10_x_CAST_EiJ
    sample : [^./ ]+

    Returns
    -------
    input_: dictionary
        fq1: string with fastq 1 path, taken from samples_table.csv
        fq2: string with fastq 2 path (if SE empty list), taken from samples_table.csv
    """
    # Get reference genome
    genomes_dir = config['genomes'][w.genomes_all]
    
    try:
        all_files = os.listdir(genomes_dir)
    except FileNotFoundError:
        sys.exit(f"Directory not found: {genomes_dir}")

    required_exts = {"sa", "pac", "bwt", "ann", "amb"}
    pattern = re.compile(r"^(?P<base>.+)\.(?P<ext>sa|pac|bwt|ann|amb)$")

    matched = []
    bases = []
    exts_found = set()

    for f in all_files:
        m = pattern.match(f)
        if m:
            matched.append(os.path.join(genomes_dir, f))
            bases.append(m.group("base"))
            exts_found.add(m.group("ext"))

    missing = required_exts - exts_found
    if missing:
        sys.exit(f"Missing BWA index files with extensions: {', '.join(missing)} in {genomes_dir}")

    if not matched:
        sys.exit(f"No BWA index files (.sa, .pac, .bwt, .ann, .amb) found in {genomes_dir}")

    # sanity check: all bases must be identical
    unique_bases = set(bases)
    if len(unique_bases) > 1:
        sys.exit(f"Multiple reference bases found: {unique_bases}")

    input_ = {
        "reference_genome_indexed_files": matched
    }

    # Get fastq files
    sample_rows = samples_table[samples_table['sample_name'] == w.sample]
    if sample_rows.empty:
        sys.exit(f"Sample {w.sample} not found in samples_table")
    lib_tech = sample_rows['library_technology'].iloc[0]
    is_pe = sample_rows['PE'].iloc[0]

    seq_tech = ""
    if lib_tech == "adaptase":
        seq_tech = ".adaptase_trimmed"

    if is_pe:
        fq1 = (f"Results/{w.sample}.adap_trimmed{seq_tech}.R1."
               "PE.fq.gz")
        fq2 = (f"Results/{w.sample}.adap_trimmed{seq_tech}.R2."
               "PE.fq.gz")
    else:
        fq1 = (f"Results/{w.sample}.adap_trimmed{seq_tech}.R1."
               "SE.fq.gz")
        fq2 = []

    input_ |= {"fq1": fq1, "fq2": fq2}
    return input_


def call_peaks_macs2_input(w):
    r"""Get input for call_peaks_macs2 rule.

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
    bam = get_processed_bam(w.sample, w.genomes_not_fused)
    bai = bam + ".bai"
    input_ = {"treat_bam": bam, "treat_bai": bai}
    # Get ctrl/input bam and bai
    ctrl_alias = w.peak_params.split("__")[1]
    if ctrl_alias != "no_input":
        if ctrl_alias in config['MACS2']['control']:
            ctrl_bam = config['MACS2']['control'][ctrl_alias]
            ctrl_bai = ctrl_bam + ".bai"
            input_ |= {
                "ctrl_bam": ctrl_bam,
                "ctrl_bai": ctrl_bai
            }
        elif ctrl_alias in samples_table_2['sample_name'].values:
            ctrl_bam = get_processed_bam(ctrl_alias, w.genomes_not_fused)
            ctrl_bai = ctrl_bam + ".bai"
            input_ |= {
                "ctrl_bam": ctrl_bam,
                "ctrl_bai": ctrl_bai
            }
        else:
            sys.exit(f"Control alias '{ctrl_alias}' not found in config or samples table.")
    return input_


def call_peaks_macs2_params(w):
    r"""Get parameters for rule call_peaks_macs2.

    Wildcards
    ----------
    genomes_not_fused : mm10|mm39|d6|hg19|hg38|mm10_x_CAST_EiJ
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
    if ctrl_alias == "no_input":
        ctrl = ""
    elif ctrl_alias in config['MACS2']['control'].keys():
        ctrl_bam = config['MACS2']['control'][ctrl_alias]
        ctrl = f"-c {ctrl_bam}"
    elif ctrl_alias in samples_table_2['sample_name'].values:
        ctrl_bam = get_processed_bam(ctrl_alias, w.genomes_not_fused)
        ctrl = f"-c {ctrl_bam}"
    else:
        sys.exit(f"Control alias '{ctrl_alias}' not found in config or samples table.")

    rows = samples_table[
        (samples_table['sample_name'] == w.sample) &
        (samples_table['reference_genome'] == w.genomes_not_fused)
    ]
    is_pe = rows['PE'].iloc[0] if not rows.empty else False
    if is_pe:
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


def clip_1bp_input(w):
    r"""Get input for clip_1bp rule

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


def clip_1bp_param(w):
    if re.search("inc", w.ss_SR):
        return "f"
    else:
        return "F"


def compute_matrix_outfiles_hs_input(w):
    r"""Get input for rule compute_matrix_outfiles_hs.

    Rule's output
    -------------
    "Results/{genomes_not_fused}/Analysis/"
                "Heatmaps_and_aggregate_profiles/Hotspots/{hs_region}/{strands}/"
                "{cov_params}/filenames/"
                "{sample}.{genomes_final}.q_filt.srt.nodup.mit_filt."
                "{ss_condit}{clip_strand}filename")

    Wildcards
    ----------
    hs_region : B6xCAST_(top_5000_pm_1000bp|
                         PRDM9_assymetric_hs_((invading|receiving)_strand|
                                              mm10_aligned)
                         )|
                top_5000_plus_minus_2000|asymetric_(watson|crick)_strong|
                x_non_par|autosomal_x_non_par_ctrl

    genomes_not_fused = mm10|mm39|d6|hg19|hg38|mm10_x_CAST_EiJ
    strands = "Both_strands|Single_strand/(1bp_clipped_reads|Full_length_reads)"
    cov_params = "((CPM|RPKM|none|drosNormalized)_)?bs[0-9]+_sm[0-9]+_ex[0-9]+"
    sample = [^./ ]+
    genomes_final =
    ss_condit = r"((83-163|99-147|inc_16|exc_16)\.)?"
    clip_condit = r"(((inc|exc)_16\.)?clipped_1_bp\.)?",

    Returns
    -------
    dict
        'region', path to bed file with hotspots coordinates
        'bigwig', path to bigwig file

    """
    hotspots = config['references'][w.genomes_not_fused][w.hs_region]
    bigwig = (f"Results/{w.genomes_not_fused}/Bigwigs/Coverage/{w.strands}/"
              f"{w.cov_params}/{w.sample}.{w.genomes_final}"
              ".q_filt.srt.nodup.mit_filt."
              f"{w.ss_condit}{w.clip_condit}bw"
              )

    return {"region": hotspots,
            "bigwig": bigwig}


def dros_normalization_input(w):
    r"""Get input for dros_normalization rule

    Wildcards
    ----------
    strand : (\.(83-163|99-147|inc_16|exc_16))?
    sample : [^./ ]+

    Returns
    -------
    Dictionary
        'report', path to a csv file containing the drosophila normalization
            factors.
        'bam', bam file to normalize
        'bai', index file of bam to normalize
    """
    rows = samples_table_2[
        (samples_table_2['sample_name'] == w.sample) &
        (samples_table_2['reference_genome'] == w.genomes_not_fused)
    ]
    if rows.empty:
        sys.exit(f"Sample '{w.sample}' for genome '{w.genomes_not_fused}' not found in samples_table_2")

    if w.strand == "":
        bam = rows['dedup_flt_both_strds_bam'].iloc[0]
    else:
        final_genome = rows['final_genome'].iloc[0]
        strand_clean = w.strand.lstrip(".")
        bam = (f"Results/{w.genomes_not_fused}/Bams/Single_strand/Full_length_reads/"
               f"{w.sample}.{final_genome}.q_filt.srt.nodup.mit_filt.{strand_clean}.bam")

    bai = bam + ".bai"
    return {
        "report": ("Results/d6/Analysis/drosophila_normalization/"
                   "drosophila_100K_reads/drosophila_equalization_report.tsv"),
        "bam": bam,
        "bai": bai
    }


def dros_normalization_report_input(w):
    r"""Get input for rule dros_normalization_report.

    Wildcards
    ----------
    No wildcards in this rule

    Returns
    -------
    List
        Paths to outputs of samtools_flagstat of samples to be normalized using
        drosophila'.
    """
    return samples_table_2.loc[
    samples_table_2["dros_spike_in"], "processed_flagstat_dros"
    ].tolist()


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
    rows = samples_table[
        (samples_table['sample_name'] == w.sample) &
        (samples_table['reference_genome'] == w.genomes_not_fused)
    ]
    is_pe = rows['PE'].iloc[0] if not rows.empty else False
    if is_pe:
        filter_ = "-F 3852 -f 3"
    else:
        filter_ = "-F 3844"
    return filter_


def filter_peaks_blk_grey_list_input(w):
    r"""Get input for rule filter_peaks_blk_grey_list

    Output of rule:
    ----------
    ("Results/{genomes_not_fused}/Peaks/MACS2/{peak_type}/"
         "{peak_params}/black_list/{sample}."
         "{genomes_final}{extension}.{peak_type}Peak"
    )

    Widcards
    ----------
    extension = r"(\.q_filt\.srt\.nodup\.mit_filt)?"

    Returns
    -------
    dict:
        "peaks": list with paths of all peaks from that reference genome
        "black_list": path of path of black-list bed file

    """
    peaks = (f"Results/{w.genomes_not_fused}/Peaks/MACS2/{w.peak_type}/"
             f"{w.peak_params}/{w.sample}.{w.genomes_final}{w.extension}_"
             f"peaks.{w.peak_type}Peak")

    black_list = config['references'][w.genomes_not_fused]['blacklist']

    return {"peaks": peaks,
            "black_list": black_list}


def FRIP_input(w):
    r"""Get input for FRIP rule.

    Wildcards
    ----------
    sample : [^./ ]+
    genomes_not_fused : mm10|mm39|d6|hg19|hg38|mm10_x_CAST_EiJ
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
    bam = (f"Results/{w.genomes_not_fused}/Bams/Both_strands/" +
           f"{w.sample}.{w.genomes_final}{w.extension}.bam")

    return {"peak": (f"Results/{w.genomes_not_fused}/Peaks/MACS2/{w.peak_type}"
                     f"/{w.peak_params}/blacklist_filtered/{w.sample}."
                     f"{w.genomes_final}{w.extension}.{w.peak_type}Peak"),
            "bam": bam,
            "bai": bam + ".bai"
            }


def get_rv_fw_strand_input(w):
    r"""Get input for rule get_rv_fw_strand

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


def get_strand_sep_bams_params(w):
    r"""Get parameters for rule get_strand_sep_bams

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
    genomes_not_fused : mm10|mm39|d6|hg19|hg38|mm10_x_CAST_EiJ
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
    rows = samples_table[
        (samples_table['sample_name'] == w.sample) &
        (samples_table['reference_genome'] == w.genomes_not_fused)
    ]
    is_b6xcast = rows['B6xCAST'].iloc[0] if not rows.empty else False
    if is_b6xcast:
        hotspots = config['references'][w.genomes_not_fused]['B6xCAST_pm_2000bp']
    else:
        hotspots = config['references'][w.genomes_not_fused]['all_plus_minus_2000']

    return {
        "hotspots": hotspots,
        "sample_peaks": (
            f"Results/{w.genomes_not_fused}/Peaks/MACS2/"
            f"{w.peak_type}/{w.peak_params}/blacklist_filtered/"
            f"{w.sample}.{w.genomes_final}{w.extension}.{w.peak_type}Peak"
        )
    }


def markdown_report_aggregate_profiles_input(w):
    r"""Get input for rule markdown_report_aggregate_profiles.

    Wildcards
    ----------
    genomes_not_fused = "mm10|d6|hg19|hg38|mm10_x_CAST_EiJ|mm39"
    smooth = "(smoothed|not_smoothed)+"


    Returns
    -------
    Dictionary
        'peaks_summary':
            string or list. Path to table with summary of peaks
            for all samples with the same reference genome. If no peaks were
            asked for, it delivers an empty list.
        'fastp':
            list. Path to fatp reports (output of trim_adapters_SE or
            trim_adapters_PE rule) for samples with the same reference genome
            as the the markdown report. These files will be used to get the
            number of reads of the raw FASTQ file.
        'bamfiles_reads':
            list. Paths to samtools_flagstat outputs (for samples with the same
            reference genome as the markdown report).
            These files will be used in the report to get the number of reads
            of the sample in the filtered bam file.
        'mm_top5000_asmtric_auto_XnonPAR_ag_profs'|'mm_top_5000_ag_profs':
            string. Only for mm10 or mm39. Path to '.RData' object containing the
            aggregate profiles in either "mm{XX} top 5000 hotspots", or in those
            same hotspots plus XnonPAR, autosomal and assymetric (left vs
            right of DSB) hotspots. When all hotsopots lists are asked for
            (later case), the rule that provides the R object is the
            "wrangle_X_nonPAR_asymmetric_HS" rule, as opposed to
            "process_aggregate_profiles" when it is only the top 5000.
        'top_5000_plus_minus_2000_clipped':
            string. Only for mm10. Path to '.RData' object containing the
            processed agggregate profiles of 1-bp clipped samples. If there is
            no clipped samples it delivers and empty list.
        'B6xCAST_top_5000_pm_2000bp':
            string. Only for mm10. Path to ".RData" object containing the
            aggreagte profiles of all samples aligned to B6xCAST fused genome,
            in top 5000 B6xCAST hotspots.
        'B6xCAST_PRDM9_assymetric_hs_invading_strand':
            string. Path to ".RData" object containing the aggreagte profiles
            of all samples aligned to B6xCAST fused genome, in B6xCAST hotspots
            that bind PRDM9 asymmetrically, on the invading strand.
        'B6xCAST_PRDM9_assymetric_hs_receiving_strand':
            string. same as previous but receiving/template strand.
        'B6xCAST_PRDM9_assymetric_hs_mm10_aligned':
            string. Path to ".RData" object containing the aggreagte profiles
            of all samples from B6xCAST mice, aligned to mm10 genome (as
            opossed to B6xCAST fused genome).
    """
    # Get peak_summary path
    ref_genome_peaks = samples_table.loc[
        samples_table["reference_genome"] == w.genomes_not_fused,
        'peak_ctrl_file_alias'
    ]
    if (ref_genome_peaks.notna() & (ref_genome_peaks != "-")).any():
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

    # Get processed aggregate profiles
    mm_agg_profiles = {}
    B6xCAST_agg_profiles = {}
    B6xCAST_agg_profiles_clipped = {}
    top_5000_plus_minus_2000_clipped = {}
    
    # Subset samples_table first
# Subset samples_table first
    st = samples_table[samples_table['reference_genome'] == w.genomes_not_fused]
    
    if (st['top5000_HS_heatmap'] &
        st['get_single_strand'] &
        ~st['B6xCAST']).any():
        mm_agg_profiles = {
            "mm_top5000_asmtric_auto_XnonPAR_ag_profs": (
                f"Results/{w.genomes_not_fused}/Analysis/"
                "Heatmaps_and_aggregate_profiles/Hotspots/"
                f"{config['library']['name']}_aggregate_profiles_data."
                f"{w.smooth}.RData"
            )
        }
    
    elif (st['top5000_HS_heatmap'] & ~st['B6xCAST']).any():
        mm_agg_profiles = {
            "mm_top_5000_ag_profs": (
                f"Results/{w.genomes_not_fused}/Analysis"
                "/Heatmaps_and_aggregate_profiles/Hotspots/"
                f"{config['library']['name']}_top_5000_plus_minus_2000."
                f"{w.smooth}.RData"
            )
        }
    
    if w.genomes_not_fused == "mm10" and (
       (st['top5000_HS_heatmap'] & st['B6xCAST']).any()):
        B6xCAST_agg_profiles = {
            "B6xCAST_top_5000_pm_2000bp": (
                f"Results/{w.genomes_not_fused}/Analysis"
                "/Heatmaps_and_aggregate_profiles/Hotspots/"
                f"{config['library']['name']}_B6xCAST_top_5000_pm_2000bp."
                f"{w.smooth}.RData"
            ),
            "B6xCAST_PRDM9_assymetric_hs_invading_strand": (
                "Results/mm10_x_CAST_EiJ/Analysis/"
                "Heatmaps_and_aggregate_profiles/Hotspots/"
                f"{config['library']['name']}_B6xCAST_PRDM9_assymetric_hs_"
                f"invading_strand.{w.smooth}.RData"
            ),
            "B6xCAST_PRDM9_assymetric_hs_receiving_strand": (
                "Results/mm10_x_CAST_EiJ/Analysis/"
                "Heatmaps_and_aggregate_profiles/Hotspots/"
                f"{config['library']['name']}_B6xCAST_PRDM9_assymetric_hs_"
                f"receiving_strand.{w.smooth}.RData"
            ),
            f"B6xCAST_PRDM9_assymetric_hs_{w.genomes_not_fused}_aligned": (
                f"Results/{w.genomes_not_fused}/Analysis/Heatmaps_and_aggregate_profiles/"
                f"Hotspots/{config['library']['name']}_"
                f"B6xCAST_PRDM9_assymetric_hs_{w.genomes_not_fused}_aligned."
                f"{w.smooth}.RData"),
        }
    
    # Clipped profiles
    if ((st['top5000_HS_heatmap']
         & st['B6xCAST']
         & st['get_single_strand']
         & st['Clip_reads_to_1bp_on_5_prime']).any()):
        B6xCAST_agg_profiles_clipped = {
            "B6xCAST_top_5000_pm_2000bp_clipped": (
                f"Results/{w.genomes_not_fused}/Analysis"
                "/Heatmaps_and_aggregate_profiles/Hotspots/"
                f"{config['library']['name']}"
                "_B6xCAST_top_5000_pm_2000bp_clipped."
                f"{w.smooth}.RData"
            )}
    
    if ((st['top5000_HS_heatmap'] 
         & ~st['B6xCAST']
         & st['get_single_strand']
         & st['Clip_reads_to_1bp_on_5_prime']).any()):
        top_5000_plus_minus_2000_clipped = {
            "top_5000_plus_minus_2000_clipped": (
                f"Results/{w.genomes_not_fused}/Analysis"
                "/Heatmaps_and_aggregate_profiles/Hotspots/"
                f"{config['library']['name']}"
                "_top_5000_plus_minus_2000_clipped."
                f"{w.smooth}.RData"
            )}



    # Return
    if (w.genomes_not_fused == "mm10") | (w.genomes_not_fused == "mm39"):
        return (peaks_summary | bamfiles_reads | samples_fastp |
                mm_agg_profiles | B6xCAST_agg_profiles |
                B6xCAST_agg_profiles_clipped |
                top_5000_plus_minus_2000_clipped)

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
    base_sample = w.sample.removesuffix("_MERGED")
    rows = samples_table[samples_table['sample_name'] == base_sample]
    if rows.empty:
        rows = samples_table[samples_table['sample_name'] == w.sample]
    if rows.empty:
        sys.exit(f"Sample '{w.sample}' not found in samples_table")

    merge_with_val = rows['merge_with'].dropna()
    if merge_with_val.empty:
        sys.exit(f"Sample '{w.sample}' has no merge_with value specified.")

    samples_merge = str(merge_with_val.iloc[0]).split()
    samples = [f"Results/{base_sample}.{w.genomes_all}.bam"]
    # If sample to merge with is a path, use it as is, otherwise look for the
    # raw bam file on samples_table_2
    for sample in samples_merge:
        if "/" in sample:
            samples.append(sample)
        else:
            second_bam = f"Results/{sample}.{w.genomes_all}.bam"
            samples.append(second_bam)
    return samples


def multiqc_input(w):
    r"""Get inputs for rule multiqc

    Rule output
    -----------
    directory(
            "Results/{genomes_not_fused}/Qctrl/multiqc_report_{library_name}_data"
            ),
    html = "Results/{genomes_not_fused}/Qctrl/multiqc_report_{library_name}.html"

    Wildcards used
    --------------
    genomes_not_fused : mm10|mm39|d6|hg19|hg38|mm10_x_CAST_EiJ
    library_name = config['library']['name']

    Returns
    -------
    List
        List containing file paths for FASTQ_reports, crosscorrlelation,
        insert_size_picard, picard_library_complexity, samstats and flagstat
        output files.
    """
    df = samples_table_no_merged_samples.loc[
        samples_table_no_merged_samples["reference_genome"] == w.genomes_not_fused
    ]

    PE = ["PE" if x else "SE" for x in df['PE']]

    fastp = [
        f"Results/FASTQ_reports/{nme}.{pe}.fastp.json"
        for nme, pe in zip(df['sample_name'], PE)
    ]

    insert_size = [
        f"Results/{w.genomes_not_fused}/Qctrl/{nme}/Processed_bam/{nme}"
        ".insert_size_picard.tab"
        for nme, pe in zip(df['sample_name'], df['PE']) if pe
    ]  # insert size only applies for PE samples

    library_complexity = [
        f"Results/{w.genomes_not_fused}/Qctrl/{nme}/Raw_bam/{nme}."
        "picard_library_complexity.tab"
        for nme in df['sample_name']
    ]

    samstat = [
        f"Results/{w.genomes_not_fused}/Qctrl/{nme}/Processed_bam/{nme}"
        ".samstats.txt"
        for nme in df['sample_name']
    ]

    flagstat = [
        f"Results/{w.genomes_not_fused}/Qctrl/{nme}/Processed_bam/{nme}"
        ".flagstat.txt"
        for nme in df['sample_name']
    ]

    res = (fastp + insert_size + library_complexity + samstat + flagstat)
    return res


def process_aggregate_profiles_clipped_input(w):
    r"""Get inputs for rule process_aggregate_profiles_clipped.
    Rule output
    -----------
    "Results/{genomes_not_fused}/Analysis/Heatmaps_and_aggregate_profiles/"
        "Hotspots/{libary}_{hs_region}_clipped.{smooth}.RData"

    Wildcards used
    --------------
    genomes_not_fused : mm10|mm39|d6|hg19|hg38|mm10_x_CAST_EiJ
    hs_region : B6xCAST_(top_5000_pm_1000bp|
                         PRDM9_assymetric_hs_((invading|receiving)_strand|
                                              mm10_aligned))|
                top_5000_plus_minus_2000|asymetric_(watson|crick)_strong|
                x_non_par|autosomal_x_non_par_ctrl
    smooth : (smoothed|not_smoothed)+

    Returns
    -------
    Only need ss profiles.
    """
    cov_params = (f"{config['coverage']['normalization']}_"
                  f"bs{config['coverage']['bin_size']}_"
                  f"sm{config['coverage']['smooth']}_"
                  f"ex{config['coverage']['extend_reads']}")

    strand_se = ["inc_16", "exc_16"]

    strand_pe = [a + "." + b
                 for a in ["83-163", "99-147"]
                 for b in strand_se
                 ]

    df = samples_table_2[samples_table_2['reference_genome'] == w.genomes_not_fused]
    files = set()
    for _, row in df.iterrows():
        if row.get('Clip_reads_to_1bp_on_5_prime', False):
            strands = strand_pe if row['PE'] else strand_se
            for s in strands:
                files.add(
                    f"Results/{w.genomes_not_fused}/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/"
                    f"{w.hs_region}/Single_strand/1bp_clipped_reads/"
                    f"{cov_params}/matrixes/{row['sample_name']}."
                    f"{row['final_genome']}.q_filt.srt.nodup.mit_filt."
                    f"{s}.clipped_1_bp.matrix"
                )
    return sorted(list(files))


def process_aggregate_profiles_inputs(w):
    r"""Get inputs for rule process_aggregate_profiles.

    Wildcards
    ----------
    genomes_not_fused : mm10|mm39|d6|hg19|hg38|mm10_x_CAST_EiJ
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

    Matrixes paths are taken from samples_table_2. So, the actual decision of
    which sample will have a matrix in which hotspot list and in which strands
    is actually taking place during samples_table_2 generation
    (suffixes.generate_samples_table_2()). For now, this process is only
    available for 'mm10' genome.
    """
    w.genomes_not_fused = "mm10" if w.genomes_not_fused == "mm10_x_CAST_EiJ" else w.genomes_not_fused
    mask_genome = samples_table_2['reference_genome'] == f"{w.genomes_not_fused}"
    df = samples_table_2.loc[mask_genome]
    hs_modified = re.sub(r"mm\d+(?=_aligned)", "mm", w.hs_region)
    matrixes_df = df.filter(
        regex=f"^{hs_modified}.*matrix$"
    )
    matrixes_array = matrixes_df.to_numpy().ravel()
    matrixes_list = matrixes_array[~pd.isnull(matrixes_array)].tolist()
    return matrixes_list


def summarize_peak_count_input(w):
    r"""Get input for rule summarize_peak_count.

    Wildcards
    ----------
    genomes_not_fused: mm10|d6|hg19|hg38|mm10_x_CAST_EiJ

    Returns
    -------
    peaks : dictionary
        'narrow_all': contains paths for all bed files with blacklist-
                    filtered narrow peaks (for that particular reference
                    genome),
        'broad_all': same but broad peaks,
        'narrow_hs': bed files with blacklist-filtered peaks intersected
                    with hs list (this files only exist for mm10),
        'broad_hs': same as above but broad peaks
    """
    genome_filtered = samples_table_2[
        samples_table_2['reference_genome'] == w.genomes_not_fused
    ]

    peak_types = ["narrow", "broad"]
    peaks = {}
    for peak_type in peak_types:
        col_all = f'{peak_type}_peak_bl_gr_flt'
        if col_all in genome_filtered:
            selection_criteria_all = (
                genome_filtered[col_all].notnull()
            )
            peaks[f"{peak_type}_all"] = genome_filtered.loc[
                selection_criteria_all,
                col_all
            ].values.tolist()
        else:
            peaks[f"{peak_type}_all"] = []

        col_hs = f'{peak_type}_peak_bl_gr_flt_hs_int'
        if include_hotspots and col_hs in genome_filtered:
            selection_criteria_hs = (
                genome_filtered[col_hs].notnull()
            )
            peaks[f"{peak_type}_hs"] = genome_filtered.loc[
                selection_criteria_hs,
                col_hs
            ].values.tolist()
        else:
            peaks[f"{peak_type}_hs"] = []

    return peaks


# Backwards compatibility alias
sumarize_peak_count_input = summarize_peak_count_input


def samstats_samtools_flagstat_input(w):
    r"""Get input for rules samstats and samtools_flagstat.

    Wildcards
    ----------
    sample : [^./ ]+
    genomes_not_fused : mm10|mm39|d6|hg19|hg38|mm10_x_CAST_EiJ
    bam_type : Raw_bam|Processed_bam

    Returns
    -------
    dictionary
        Bam and bai paths.
    """
    if w.genomes_not_fused != "d6":
        rows = samples_table_2[
            (samples_table_2['sample_name'] == w.sample) &
            (samples_table_2['reference_genome'] == w.genomes_not_fused)
        ]
        if rows.empty:
            sys.exit(f"Sample '{w.sample}' for genome '{w.genomes_not_fused}' not found in samples_table_2")
        if f"{w.bam_type}" == 'Raw_bam':
            bam = rows['raw_bam'].iloc[0]
        else:
            bam = rows['dedup_flt_both_strds_bam'].iloc[0]
    else:
        rows = samples_table[samples_table['sample_name'] == w.sample]
        if rows.empty:
            sys.exit(f"Sample '{w.sample}' not found in samples_table")
        genome = rows['reference_genome'].iloc[0]
        if f"{w.bam_type}" == 'Raw_bam':
            bam = f"Results/{w.sample}.{genome}_f_d6.d6.bam"
        else:
            bam = (f"Results/d6/Bams/Both_strands/{w.sample}.{genome}_f_d6.d6."
                   "q_filt.srt.nodup.mit_filt.bam")
    bai = f"{bam}.bai"
    return {
        "bam": bam,
        "bai": bai
    }


def trim_adapters_PE_input(w):
    r"""Get FASTQ paths for trim_adapters_PE rule

    wildcards
    ----------
    sample : [^./ ]+

    Returns
    -------
    input_: dictionary
        Fastq/s path/s, taken from samples_table.csv
    """
    return {
        "fastq1": get_sample_fastq1(w.sample),
        "fastq2": get_sample_fastq2(w.sample)
    }
