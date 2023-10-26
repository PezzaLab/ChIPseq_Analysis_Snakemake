#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 15:04:24 2023

@author: quio

This module contains only one function called generate_samples_table_2() that
returns a modified samples_table dataframe containing new columns with
files' full paths (so they can be called upon by snakemake pipeline
using the sample and column name).
"""

import pandas as pd
import numpy as np
import re


def generate_samples_table_2(samples_table, config):
    """Generates different file paths to be used by snakemake and
    smk_functions.py

    Parameters
    -----------
    samples_table: dataframe
        dataframe with samples information taken from file samples_table.csv

    config: snakemake object
        snakemake object with all info from config file. It is accessed in
        snakemake as 'config'

    Returns
    ------
    samples_table_2: dataframe
        samples_table_2 is the same table as samples_table, but with extra
        columns that contain paths to different files for each sample.
    """
    samples_table_2 = samples_table.copy()

    # %% Suffixes to make file names
    nodup_filt = ".q_filt.srt.nodup.mit_filt"

    aligned_genomes = []
    for i in samples_table.index:
        genome = samples_table_2.loc[i, 'reference_genome']
        if samples_table_2.loc[i, 'B6xCAST']:
            genome = "mm10_x_CAST_EiJ"
        if samples_table_2.loc[i, 'dros_spike_in']:
            genome += "_f_d6"
        aligned_genomes.append(genome)

    final_genomes = []
    for ag, dsi in zip(
            aligned_genomes,
            samples_table_2['dros_spike_in']):
        finalg = ag
        if dsi:
            finalg = ag + "." + ag.removesuffix("_f_d6")
        final_genomes.append(finalg)

    samples_nodup_filt = (
        samples_table_2['sample_name']
            + "."
            + final_genomes
            + nodup_filt)

    cov_params = (f"{config['coverage']['normalization']}_"
                  f"bs{config['coverage']['bin_size']}_"
                  f"sm{config['coverage']['smooth']}_"
                  f"ex{config['coverage']['extend_reads']}"
                  "/"
    )

    dros_norm_cov_config_params_string = (
        f"drosNormalized_"
        f"bs{config['coverage']['bin_size']}_"
        f"sm{config['coverage']['smooth']}_"
        f"ex{config['coverage']['extend_reads']}"
        "/"
    )

    # %% Create columns with file names
    # %%% Dros bam
    samples_table_2['dros_dedup_bam'] = (
        np.where(samples_table_2['dros_spike_in'],
                 "Results/d6/Bams/Both_strands/"
                     + samples_table_2['sample_name']
                     + aligned_genomes
                     + ".d6"
                     + nodup_filt
                     + ".bam",
                 np.NaN)
    )

    # Both strands bam
    samples_table_2['raw_bam'] = (
        "Results/" + samples_table_2['sample_name']
        + "." + aligned_genomes + ".bam"
    )
    samples_table_2['dedup_flt_both_strds_bam'] = (
        "Results/" + samples_table_2['reference_genome']
        + "/Bams/Both_strands/" + samples_nodup_filt + ".bam"
    )

    # %%% Peak files
    if pd.notna(samples_table['peak_ctrl_file_alias']).any():
        peak_types = ["narrow", "broad"]
        peak_params = (
            "qv_" + config['MACS2']['qvalue'].split(".")[1] +
            "__" + samples_table_2['peak_ctrl_file_alias']
        )

        for peak_type in peak_types:
            # Modify default values accordingly
            if peak_type == "broad":
                peak_params = (
                    "bco_" + config['MACS2']['broad_cutoff'].split(".")[1]
                    + "_qv_" + config['MACS2']['qvalue'].split(".")[1] + "__"
                    + samples_table_2['peak_ctrl_file_alias']
                )

            # Generate column name and content
            samples_table_2[f'{peak_type}_peak_bl_gr_flt'] = np.where(
                samples_table['peak_ctrl_file_alias'] != "-",
                "Results/"
                    + samples_table_2['reference_genome']
                    + f"/Peaks/MACS2/{peak_type}/"
                    + peak_params
                    + "/Black-grey_filtered/"
                    + samples_nodup_filt
                    + f".{peak_type}Peak",
                np.NaN
            )
            samples_table_2[f'{peak_type}_peak_bl_gr_flt_hs_int'] = np.where(
                (samples_table_2["reference_genome"] == "mm10") &
                    (samples_table['peak_ctrl_file_alias'] != "-"),
                "Results/" 
                    + samples_table_2['reference_genome']
                    + f"/Peaks/MACS2/{peak_type}/"
                    + peak_params
                    + "/Black-grey_filtered/Intersect_HSs_plus_minus_2000_bp/"
                    + samples_nodup_filt + f".{peak_type}Peak",
                np.NaN
            )

            # FRIP
            samples_table_2[f'{peak_type}_blk_gr_flt_FRIP'] = np.where(
                samples_table['peak_ctrl_file_alias'] != "-",
                "Results/" 
                    + samples_table_2['reference_genome']
                    + "/Qctrl/" + samples_table_2['sample_name']
                    + f"/Processed_bam/FRIP/MACS2_{peak_type}_"
                    + peak_params
                    + "_bl-gr_flt/"
                    + samples_nodup_filt
                    + ".FRIP.txt",
                np.NaN
            )

            # Annotated peaks
            samples_table_2[
                f'{peak_type}Peak_blk_gr_flt_annotated'
            ] = np.where(
                    samples_table['peak_ctrl_file_alias'] != "-",
                    "Results/"
                        + samples_table_2['reference_genome']
                        + f"/Peaks/MACS2/{peak_type}/"
                        + peak_params
                        + "/Black-grey_filtered/Annotated_peaks/"
                        + samples_nodup_filt
                        + f".{peak_type}Peak.annotated.tsv",
                    np.NaN)

    # %%% samtools flagstat
    samples_table_2['raw_flagstat'] = (
        "Results/"
        + samples_table_2['reference_genome']
        + "/Qctrl/"
        + samples_table_2['sample_name']
        + "/Raw_bam/"
        + samples_table_2['sample_name']
        + ".flagstat.txt"
    )
    samples_table_2['processed_flagstat'] = (
        "Results/"
        + samples_table_2['reference_genome']
        + "/Qctrl/"
        + samples_table_2['sample_name']
        + "/Processed_bam/"
        + samples_table_2['sample_name']
        + ".flagstat.txt"
    )
    samples_table_2['processed_flagstat_dros'] = np.where(
        samples_table["dros_spike_in"],
        "Results/d6/Qctrl/"
        + samples_table_2['sample_name']
        + "/Processed_bam/"
        + samples_table_2['sample_name']
        + ".flagstat.txt",
        np.NaN
    )

    # %%% Hotspot heatmaps (only for top 5000 hs, both strands)
    samples_table_2['heatmap_top_5000_hs'] = np.where(
        samples_table_2["top5000_HS_heatmap"] & ~samples_table_2["B6xCAST"],
        "Results/"
            + samples_table_2['reference_genome']
            + "/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/"
            + "top_5000_plus_minus_2000/Both_strands/"
            + cov_params
            + "Heatmaps/"
            + samples_nodup_filt
            + ".png",
        np.NaN
    )
    samples_table_2['heatmap_B6xCAST_top_5000_hs'] = np.where(
        samples_table_2["top5000_HS_heatmap"] & samples_table_2["B6xCAST"],
        "Results/"
            + samples_table_2['reference_genome']
            + "/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/"
            + "B6xCAST_top_5000_pm_1000bp/Both_strands/"
            + cov_params
            + "Heatmaps/"
            + samples_nodup_filt
            + ".png",
        np.NaN
    )

    # %%% Drosophila normalized bigwigs
    samples_table_2['dros_normalized_both_strands_coverage_bw'] = np.where(
        ~samples_table["dros_equalization_group"].isnull(),
        "Results/"
            + samples_table['reference_genome']
            + "/Bigwigs/Coverage/Both_strands/"
            + dros_norm_cov_config_params_string
            + samples_table['dros_equalization_group']
            + "/"
            + samples_nodup_filt
            + ".dros_norm.bw",
        np.NaN)

    # %%% Coverage BIGWIG files
    # These are all the potential coverage files:
    # snwe == "sample name with extensions", flag extensions are excluded
        # Results/{genome}/Bigwigs/Coverage/
            # Both_strands/
                # {cov_params}/{snwe}.bw
                # {dros_cov_params}/{dros_eq_group}/{snwe}.bw
            # Single_strand/
                # Full_length_reads/
                    # {dros_cov_params}/{snwe}.
                        # 83-163.bw  --> for PE
                        # 99-147.bw  --> for PE
                        # inc_16.bw  --> for SE
                        # exc_16.bw  --> for SE
                # 1bp_clipped_reads/
                    # {dros_cov_params}/{snwe}.
                        # 83-163.inc_16.clipped_1_bp.bw  --> for PE
                        # 99-147.inc_16.clipped_1_bp.bw  --> for PE
                        # 83-163.exc_16.clipped_1_bp.bw  --> for PE
                        # 99-147.exc_16.clipped_1_bp.bw  --> for PE
                        # inc_16.clipped_1_bp.bw  --> for SE
                        # exc_16.clipped_1_bp.bw  --> for SE
                
    strands = ["83-163", "99-147", "inc_16", "exc_16", "Both_strands"]
    
    for strand in strands:
        # Strings
        strand_string = "Both_strands/"
        genome = samples_table["reference_genome"]
        strand2 = "." + strand
        if strand == "Both_strands":
            strand2 = ""
            # "Both_strands" is not included in the name of the file

        # Booleans
            # For B6xCAST samples both coverage files will be asked (mm10 and
            # B6_f_CAST aligned). For clipped profiles is the same, both
            # clipped and unclipped will be asked. Clipped will only be asked
            # for ssDNA.
        PE = True
        single_strand = True
        clipped = False

        # Modify default values accordingly
        if re.match("83|99|inc_16|exc_16", strand):
            strand_string = "Single_strand/Full_length_reads/"
            single_strand = samples_table["get_single_strand"]
            clipped = samples_table['Clip_reads_to_1bp_on_5_prime']  # Only do
            # clipped for ssDNA

        if re.match("83|99", strand):
            PE = samples_table["PE"]
            single_strand = samples_table["get_single_strand"]

        elif re.match("inc_16|exc_16", strand):
            PE = ~samples_table["PE"]
            single_strand = samples_table["get_single_strand"]

        # Generate column name and content for non-clipped non-B6xCAST
        samples_table_2[f"{strand}_coverage_bw"] = np.where(
            PE & single_strand,
            "Results/" + genome + "/Bigwigs/Coverage/"
                + f"{strand_string}{cov_params}"
                + samples_nodup_filt
                + strand2
                + ".bw",
            np.NaN)

        # Generate column name and content for non-clipped B6xCAST
        samples_table_2[f"B6xCAST_{strand}_coverage_bw"] = np.where(
            PE & single_strand & samples_table['B6xCAST'],
            "Results/mm10_x_CAST_EiJ"
                + "/Bigwigs/Coverage/"
                + f"{strand_string}{cov_params}"
                + samples_nodup_filt
                + strand2
                + ".bw",
            np.NaN)
            # If I don't have a series within the objects to concat, the
            # following error arises:
            # "can only concatenate str (not "list") to str"

        # Generate column name and content for clipped profiles
        ## non-B6xCAST
        if re.match("83|99", strand):
            for strand3 in ["inc_16", "exc_16"]:
                samples_table_2[f"{strand}_{strand3}_clipped_coverage_bw"] = (
                    np.where(
                        PE & single_strand & clipped,
                        "Results/"
                            + genome
                            + "/Bigwigs/Coverage/"
                            + f"Single_strand/1bp_clipped_reads/{cov_params}"
                            + samples_nodup_filt
                            + strand2
                            + "."
                            + strand3
                            + ".clipped_1_bp.bw",
                        np.NaN
                    )
                )
        elif re.match("(inc|exc)_16", strand):
            samples_table_2[f"{strand}_clipped_coverage_bw"] = (
                np.where(
                    PE
                        & single_strand 
                        & clipped,
                    "Results/"
                        + genome
                        + "/Bigwigs/Coverage/"
                        + f"Single_strand/1bp_clipped_reads/{cov_params}"
                        + samples_nodup_filt
                        + strand2
                        + ".clipped_1_bp.bw",
                    np.NaN
                )
            )


        ## B6xCAST
        if re.match("83|99", strand):
            for strand3 in ["inc_16", "exc_16"]:
                samples_table_2[f"B6xCAST_{strand}_{strand3}_clipped_coverage_bw"] = (
                    np.where(
                        PE
                            & single_strand
                            & clipped
                            & samples_table['B6xCAST'],
                        "Results/mm10_x_CAST_EiJ"
                            + "/Bigwigs/Coverage/"
                            + f"Single_strand/1bp_clipped_reads/{cov_params}"
                            + samples_nodup_filt
                            + strand2
                            +"."
                            + strand3
                            + ".clipped_1_bp.bw",
                        np.NaN
                    )
                )
        elif re.match("(inc|exc)_16", strand):
            samples_table_2[f"B6xCAST_{strand}_clipped_coverage_bw"] = (
                np.where(
                    PE
                        & single_strand
                        & clipped
                        & samples_table['B6xCAST'],
                    "Results/mm10_x_CAST_EiJ" 
                        + "/Bigwigs/Coverage/"
                        + f"Single_strand/1bp_clipped_reads/{cov_params}"
                        + samples_nodup_filt
                        + strand2
                        + ".clipped_1_bp.bw",
                    np.NaN
                )
            )

    # %%% Coverage MATRIX files
    # This files are generated by deeptools using a bw and a bed file
    hotspots = ["top_5000_plus_minus_2000",
                "x_non_par",
                "autosomal_x_non_par_ctrl",
                "asymetric_watson_strong",
                "asymetric_crick_strong",
                "B6xCAST_top_5000_pm_1000bp",
                "B6xCAST_PRDM9_assymetric_hs_invading_strand",
                "B6xCAST_PRDM9_assymetric_hs_receiving_strand",
                "B6xCAST_PRDM9_assymetric_hs_mm10_aligned"
                ]

    strands = ["83-163",
               "99-147",
               "inc_16",
               "exc_16",
               "Both_strands"]

    for hs in hotspots:
        for strand in strands:
            # Strings
            strand_string = "Both_strands/"
            genome = samples_table["reference_genome"]
            strand2 = "." + strand
            samples_nodup_filt2 = samples_nodup_filt
            if strand == "Both_strands":
                strand2 = ""  # "Both_strands" is not included in the name of
                # the file

            # Booleans
            PE = True
            single_strand = True
            b6xcast = ~samples_table_2["B6xCAST"]
            # Modify default values accordingly
            if re.search("B6xCAST", hs):
                b6xcast = samples_table_2["B6xCAST"]

            if re.match("83|99|inc_16|exc_16", strand):
                strand_string = "Single_strand/Full_length_reads/"
                single_strand = samples_table_2["get_single_strand"]

            if re.match("watson|crick|x_non_par|", hs):
                single_strand = samples_table_2["get_single_strand"]
                # Only profiles in which ssDNA is asked for we evaluate the
                # XnonPAR and watson/crick assymetric hotspots (mm10 samples)

            if re.match("83|99", strand):
                PE = samples_table_2["PE"]
            elif re.match("inc_16|exc_16", strand):
                PE = ~samples_table_2["PE"]

            if re.match("B6xCAST_PRDM9_assymetric_hs_(invading|receiving)",
                        hs):
                genome = "mm10_x_CAST_EiJ"
                samples_nodup_filt2 = samples_nodup_filt

            # Generate column name and content
            samples_table_2[f"{hs}_{strand}_matrix"] = np.where(
                samples_table_2["top5000_HS_heatmap"] 
                    & b6xcast 
                    & PE 
                    & single_strand,
                "Results/"
                    + genome
                    + "/Analysis/"
                    + "Heatmaps_and_aggregate_profiles/Hotspots/"
                    + hs
                    + "/"
                    + f"{strand_string}{cov_params}matrixes/"
                    + samples_nodup_filt2 + strand2 + ".matrix",
                np.NaN)
    # %% Return
    return samples_table_2
