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
import os
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

    aligned_genomes = (
        np.where(samples_table_2['dros_spike_in'],
                 samples_table_2['reference_genome'] + "_f_d6",
                 samples_table_2['reference_genome'])
    )

    final_genomes = (
        np.where(
            samples_table_2['dros_spike_in'],
           aligned_genomes + "." + samples_table_2['reference_genome'],
           aligned_genomes)
    )
    samples_table_2['final_genome'] = final_genomes

    final_genomes_b6xcast = []
    for row in samples_table_2.itertuples():
        ref = row.reference_genome
        if row.dros_spike_in:
            final_genomes_b6xcast.append(
                f"{ref}_x_CAST_EiJ_f_d6.{ref}_x_CAST_EiJ"
            )
        else:
            final_genomes_b6xcast.append(f"{ref}_x_CAST_EiJ")
    
    samples_table_2["final_genome_b6xcast"] = final_genomes_b6xcast

    samples_nodup_filt = (
        samples_table_2['sample_name'] + "." + final_genomes + nodup_filt
    )

    samples_nodup_filt_b6xcast = (
        samples_table_2['sample_name']
            + "."
            + final_genomes_b6xcast
            + nodup_filt
    )

    cov_params = (
        f"{config['coverage']['normalization']}_"
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
                     + "." + aligned_genomes
                     + ".d6"
                     + nodup_filt
                     + ".bam",
                 np.nan)
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
    peak_types = ["narrow", "broad"]
    for peak_type in peak_types:
        samples_table_2[f'{peak_type}_peak_bl_gr_flt'] = np.nan
        samples_table_2[f'{peak_type}_peak_bl_gr_flt_hs_int'] = np.nan
        samples_table_2[f'{peak_type}_blk_gr_flt_FRIP'] = np.nan

    has_ctrl = (
        samples_table_2['peak_ctrl_file_alias'].notna() &
        (samples_table_2['peak_ctrl_file_alias'] != "-")
    )
    if has_ctrl.any():
        qv  = config['MACS2']['qvalue'].split(".")[1]
        bco = config['MACS2']['broad_cutoff'].split(".")[1]
        is_hs_genome = samples_table_2['reference_genome'].isin(["mm10", "mm39"])

        for peak_type in peak_types:
            # per-type peak_params
            if peak_type == "narrow":
                peak_params_pt = ("qv_" + qv + "__" 
                + samples_table_2['peak_ctrl_file_alias']
                )
            else:  # broad
                peak_params_pt = ("bco_" + bco + "_qv_" + qv + "__" 
                + samples_table_2['peak_ctrl_file_alias']
                )

            # Blacklist filtered peak path
            samples_table_2[f'{peak_type}_peak_bl_gr_flt'] = np.where(
                has_ctrl,
                "Results/" + samples_table_2['reference_genome']
                + f"/Peaks/MACS2/{peak_type}/"
                + peak_params_pt
                + "/blacklist_filtered/"
                + samples_nodup_filt
                + f".{peak_type}Peak",
                np.nan
            )

            # HS intersect path (now for mm10 OR mm39)
            samples_table_2[f'{peak_type}_peak_bl_gr_flt_hs_int'] = np.where(
                has_ctrl & is_hs_genome,
                "Results/" + samples_table_2['reference_genome']
                + f"/Peaks/MACS2/{peak_type}/"
                + peak_params_pt
                + "/blacklist_filtered/Intersect_HSs_plus_minus_2000_bp/"
                + samples_nodup_filt
                + f".{peak_type}Peak",
                np.nan
            )

            # FRIP
            samples_table_2[f'{peak_type}_blk_gr_flt_FRIP'] = np.where(
                has_ctrl,
                "Results/" + samples_table_2['reference_genome']
                + "/Qctrl/" + samples_table_2['sample_name']
                + f"/Processed_bam/FRIP/MACS2_{peak_type}_"
                + peak_params_pt
                + "_bl-gr_flt/"
                + samples_nodup_filt
                + ".FRIP.txt",
                np.nan
            )


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
        np.nan
    )

    # %%% Hotspot heatmaps (only for top 5000 hs, both strands)
    samples_table_2['heatmap_top_5000_hs'] = np.where(
        samples_table_2["top5000_HS_heatmap"] & ~samples_table_2["B6xCAST"],
        "Results/" + samples_table_2['reference_genome'] +
        "/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/"
        "top_5000_plus_minus_2000/Both_strands/"
        + cov_params
        + "Heatmaps/"
        + samples_nodup_filt
        + ".png",
        np.nan
    )
    samples_table_2['heatmap_B6xCAST_top_5000_hs'] = np.where(
        samples_table_2["top5000_HS_heatmap"] & samples_table_2["B6xCAST"],
        "Results/"
            + samples_table_2['reference_genome']
            + "/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/"
            + "B6xCAST_top_5000_pm_2000bp/Both_strands/"
            + cov_params
            + "Heatmaps/"
            + samples_nodup_filt
            + ".png",
        np.nan
    )

    # %%% Drosophila normalized bigwigs
    samples_table_2['dros_normalized_both_strands_coverage_bw'] = np.where(
        samples_table["dros_spike_in"],
        "Results/"
            + samples_table['reference_genome']
            + "/Bigwigs/Coverage/Both_strands/"
            + dros_norm_cov_config_params_string
            + "drosophila_100K_reads/"
            + samples_nodup_filt
            + ".dros_norm.bw",
        np.nan)

    # %%% Coverage BIGWIG files
    # These are all the potential coverage files:
    # snwe == "sample name with extensions", flag extensions are excluded
        # Results/{genome}/Bigwigs/Coverage/
            # Both_strands/
                # {cov_params}/{snwe}.bw
                # {dros_cov_params}/drosophila_100K_reads/{snwe}.bw
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
        B6xCAST = samples_table['B6xCAST']

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
            np.nan)

        # Generate column name and content for non-clipped B6xCAST
        samples_table_2[f"B6xCAST_{strand}_coverage_bw"] = np.where(
            PE 
                & single_strand 
                & samples_table['B6xCAST'] 
                & (samples_table['reference_genome'] == "mm10"),
            "Results/""mm10_x_CAST_EiJ" + "/Bigwigs/Coverage/"
                + f"{strand_string}{cov_params}"
                + samples_nodup_filt_b6xcast
                + strand2
                + ".bw",
            np.nan)

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
                            +"."
                            + strand3
                            + ".clipped_1_bp.bw",
                        np.nan
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
                    np.nan
                )
            )


        ## B6xCAST
        if re.match("83|99", strand):
            for strand3 in ["inc_16", "exc_16"]:
                samples_table_2[
                    f"B6xCAST_{strand}_{strand3}_clipped_coverage_bw"
                ] = (
                    np.where(
                        PE
                            & single_strand
                            & clipped
                            & samples_table['B6xCAST']
                            & (samples_table['reference_genome'] == "mm10"),
                        "Results/mm10_x_CAST_EiJ" + "/Bigwigs/Coverage/"
                            + f"Single_strand/1bp_clipped_reads/{cov_params}"
                            + samples_nodup_filt
                            + strand2
                            +"."
                            + strand3
                            + ".clipped_1_bp.bw",
                        np.nan
                    )
                )
        elif re.match("(inc|exc)_16", strand):
            samples_table_2[f"B6xCAST_{strand}_clipped_coverage_bw"] = (
                np.where(
                    PE
                        & single_strand
                        & clipped
                        & samples_table['B6xCAST'],
                    "Results/mm10_x_CAST_EiJ" + "/Bigwigs/Coverage/"
                        + f"Single_strand/1bp_clipped_reads/{cov_params}"
                        + samples_nodup_filt
                        + strand2
                        + ".clipped_1_bp.bw",
                    np.nan
                )
            )

   # %%% Coverage MATRIX files
    # Matrix files are generated by deeptools using a bw and a bed file
    hotspots = [
        "top_5000_plus_minus_2000",
        "x_non_par",
        "autosomal_x_non_par_ctrl",
        "asymetric_watson_strong",
        "asymetric_crick_strong",
        "B6xCAST_top_5000_pm_2000bp",
        "B6xCAST_PRDM9_assymetric_hs_invading_strand",   # only mm10
        "B6xCAST_PRDM9_assymetric_hs_receiving_strand",  # only mm10
        "B6xCAST_PRDM9_assymetric_hs_mm_aligned"         # aligned to either mm10 or mm39
    ]
    
    strands = ["83-163", "99-147", "inc_16", "exc_16", "Both_strands"]
    
    for hs in hotspots:
        for strand in strands:
            # Strings
            strand_string = "Both_strands/"
            strand2 = "." + strand
            samples_nodup_filt2 = samples_nodup_filt
            if strand == "Both_strands":
                strand2 = ""  # "Both_strands" is not included in the file name
    
            # Booleans
            PE = True
            single_strand = True
            b6xcast = ~samples_table_2["B6xCAST"]
    
            # Modify defaults
            if re.search("B6xCAST", hs):
                b6xcast = samples_table_2["B6xCAST"]
    
            if re.search("83|99|inc_16|exc_16", strand):
                strand_string = "Single_strand/Full_length_reads/"
                single_strand = samples_table_2["get_single_strand"]
    
            if re.search("watson|crick|x_non_par", hs):
                # only evaluate these when ssDNA is requested
                single_strand = samples_table_2["get_single_strand"]
    
            if re.search("83|99", strand):
                PE = samples_table_2["PE"]
            elif re.search("inc_16|exc_16", strand):
                PE = ~samples_table_2["PE"]
    
            # Default base directory for Results (row-wise genome)
            results_root = samples_table_2["reference_genome"]
    
            # Special case: PRDM9 asymmetric (invading/receiving) → use fused genome paths + b6xcast files
            if re.search(r"B6xCAST_PRDM9_assymetric_hs_(invading|receiving)", hs):
                b6xcast = (
                    (samples_table_2["reference_genome"] == "mm10")
                    & (samples_table_2["B6xCAST"])
                )
                results_root = "mm10_x_CAST_EiJ"
                samples_nodup_filt2 = samples_nodup_filt_b6xcast
    
            # Adjust hs only for path building (mm_aligned → mm10_aligned or mm39_aligned)
            if hs == "B6xCAST_PRDM9_assymetric_hs_mm_aligned":
                hs_path = "B6xCAST_PRDM9_assymetric_hs_" + samples_table_2["reference_genome"] + "_aligned"
            else:
                hs_path = hs
    
            # Generate column name and content
            samples_table_2[f"{hs}_{strand}_matrix"] = np.where(
                samples_table_2["top5000_HS_heatmap"] & b6xcast & PE & single_strand,
                "Results/"
                + results_root
                + "/Analysis/"
                + "Heatmaps_and_aggregate_profiles/Hotspots/"
                + hs_path
                + "/"
                + f"{strand_string}{cov_params}matrixes/"
                + samples_nodup_filt2
                + strand2
                + ".matrix",
                np.nan
            )

            

    # %% Return
    return samples_table_2
