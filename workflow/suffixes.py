#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 21 15:04:24 2023

@author: quio
"""

import pandas as pd
import numpy as np
import os
import re

def generate_samples_table_2(samples_table, config):
    mm10_samples_table = samples_table.loc[samples_table['reference_genome'] == "mm10"]

    # For procedures that are only applied to mm10 (such as HS intersection)
    number_samples = len(samples_table)
    
    samples_table_2 = samples_table.copy()
    
    ### Suffixes to make file names ###
    nodup_filt = ".q_filt.srt.nodup.mit_filt"
    
    aligned_genomes = (
        np.where(samples_table_2['dros_spike_in'],
                 "." + samples_table_2['reference_genome'] + "_f_d6",
                 "." + samples_table_2['reference_genome'])
    )
    
    final_genomes = (
        np.where(samples_table_2['dros_spike_in'],
                 aligned_genomes + "." + samples_table_2['reference_genome'],
                 aligned_genomes)
    )
    
    samples_nodup_filt = (
        samples_table_2['sample_name'] + final_genomes +
        nodup_filt
    )
    
    cov_params = f"{config['coverage']['normalization']}_bs{config['coverage']['bin_size']}_sm{config['coverage']['smooth']}/"
    
    peak_nrw_params = (
        "qv_" + (config['MACS2']['qvalue'].split("."))[1] +
        "__" + samples_table_2['peak_ctrl_file_alias']
    )
    
    peak_brd_params = (
        "bco_" + (config['MACS2']['broad_cutoff'].split("."))[1] + "_qv_" +
        (config['MACS2']['qvalue'].split("."))[1] + "__" +
        samples_table_2['peak_ctrl_file_alias']
    )
    
    #### Create columns with file names #####
    # Dros bam
    samples_table_2['dros_dedup_bam'] = (
        np.where(samples_table_2['dros_spike_in'],
                 ("Results/d6/Bams/Both_strands/" + samples_table_2['sample_name'] +
                  aligned_genomes + ".d6" + nodup_filt + ".bam"),
                 np.NaN)
    )
    
    # Both strands bam
    samples_table_2['raw_bam'] = (
        "Results/" + samples_table_2['sample_name'] +
        aligned_genomes + ".bam"
    )
    samples_table_2['dedup_flt_both_strds_bam'] = (
        "Results/" + samples_table_2['reference_genome'] +
        "/Bams/Both_strands/" + samples_nodup_filt + ".bam"
    )
    
    # Single strand bam files, unclipped
    conditions_ss_83_or_i16_bam = [
        samples_table_2['PE'] & samples_table_2['get_single_strand'],
        ~samples_table_2['PE'] & samples_table_2['get_single_strand']
    ]
    choices_ss_83_or_i16_bam = [
        ("Results/" + samples_table_2['reference_genome'] +
         "/Bams/Single_strand/Full_length_reads/" +
         samples_nodup_filt + ".83-163.bam"),
        "Results/" + samples_table_2['reference_genome'] +
        "/Bams/Single_strand/Full_length_reads/" +
        samples_nodup_filt + ".inc_16.bam"
    ]
    
    samples_table_2['ss_83_or_i16_bam'] = (
        np.select(conditions_ss_83_or_i16_bam,
                  choices_ss_83_or_i16_bam,
                  default=np.NaN)
    )
    
    samples_table_2['ss_99_or_e16_bam'] = (
        np.where(samples_table_2['PE'],
                 "Results/" + samples_table_2['reference_genome'] +
                 "/Bams/Single_strand/Full_length_reads/" +
                 samples_nodup_filt + ".99-147.bam",
                 "Results/" + samples_table_2['reference_genome'] +
                 "/Bams/Single_strand/Full_length_reads/" +
                 samples_nodup_filt + ".exc_16.bam"
                 )
    )
    
    # Single strand bam, clipped profiles
    # PE data
    conditions_ss_not_clpd_bam = [
        samples_table_2['PE'] & samples_table_2['get_single_strand'],
        ~samples_table_2['PE'] & samples_table_2['get_single_strand']
    ]
    choices_ss_83_or_i16_bam = [
        ("Results/" + samples_table_2['reference_genome'] +
         "/Bams/Single_strand/Full_length_reads/" +
         samples_nodup_filt + ".83-163.bam"),
        "Results/" + samples_table_2['reference_genome'] +
        "/Bams/Single_strand/Full_length_reads/" +
        samples_nodup_filt + ".inc_16.bam"
    ]
    choices_ss_99_or_e16_bam = [
        ("Results/" + samples_table_2['reference_genome'] +
         "/Bams/Single_strand/Full_length_reads/" +
         samples_nodup_filt + ".99-147.bam"),
        "Results/" + samples_table_2['reference_genome'] +
        "/Bams/Single_strand/Full_length_reads/" +
        samples_nodup_filt + ".exc_16.bam"
    ]
    
    samples_table_2['ss_83_or_i16_bam'] = (
        np.select(conditions_ss_not_clpd_bam,
                  choices_ss_83_or_i16_bam,
                  default=np.NaN)
    )
    samples_table_2['ss_99_or_e16_bam'] = (
        np.select(conditions_ss_not_clpd_bam,
                  choices_ss_99_or_e16_bam,
                  default=np.NaN)
    )
    
    # Clipped
    samples_table_2['ssPE_83_i16_clpd_bam'] = np.where(
        samples_table_2['PE'] & samples_table_2['Clip_reads_to_1bp_on_5_prime'],
        ("Results/" + samples_table_2['reference_genome'] +
         "/Bams/Single_strand/1bp_clipped_reads/" + samples_nodup_filt +
         ".83-163.inc_16.clipped_1_bp.bam"),
        np.NaN)
    
    samples_table_2['ssPE_83_e16_clpd_bam'] = np.where(
        samples_table_2['PE'] & samples_table_2['Clip_reads_to_1bp_on_5_prime'],
        ("Results/" + samples_table_2['reference_genome'] +
         "/Bams/Single_strand/1bp_clipped_reads/" + samples_nodup_filt +
         ".83-163.exc_16.clipped_1_bp.bam"),
        np.NaN)
    
    samples_table_2['ssPE_99_i16_clpd_bam'] = np.where(
        samples_table_2['PE'] & samples_table_2['Clip_reads_to_1bp_on_5_prime'],
        ("Results/" + samples_table_2['reference_genome'] +
         "/Bams/Single_strand/1bp_clipped_reads/" + samples_nodup_filt +
         ".99-147.inc_16.clipped_1_bp.bam"),
        np.NaN)
    
    
    samples_table_2['ssPE_99_e16_clpd_bam'] = np.where(
        samples_table_2['PE'] & samples_table_2['Clip_reads_to_1bp_on_5_prime'],
        ("Results/" + samples_table_2['reference_genome'] +
         "/Bams/Single_strand/1bp_clipped_reads/" + samples_nodup_filt +
         ".99-147.exc_16.clipped_1_bp.bam"),
        np.NaN)
    
    # SR data
    samples_table_2['ssSR_i16_clp_bam'] = np.where(
        ~samples_table_2['PE'] & samples_table_2['Clip_reads_to_1bp_on_5_prime'],
        ("Results/" + samples_table_2['reference_genome'] +
         "/Bams/Single_strand/1bp_clipped_reads/" + samples_nodup_filt +
         ".inc_16.clipped_1_bp.bam"),
        np.NaN)
    samples_table_2['ssSR_e16_clp_bam'] = np.where(
        ~samples_table_2['PE'] & samples_table_2['Clip_reads_to_1bp_on_5_prime'],
        ("Results/" + samples_table_2['reference_genome'] +
         "/Bams/Single_strand/1bp_clipped_reads/" + samples_nodup_filt +
         ".exc_16.clipped_1_bp.bam"),
        np.NaN)
    
    # Coverage files
    # Both strands bw
    samples_table_2['dedup_flt_both_strds_bw'] = (
        "Results/" + samples_table_2['reference_genome'] +
        "/Bigwigs/Coverage/Both_strands/" + cov_params + samples_nodup_filt + ".bw"
    )
    
    # Single strand bw files, unclipped
    conditions_ss_83_or_i16_bw = [
        samples_table_2['PE'] & samples_table_2['get_single_strand'],
        ~samples_table_2['PE'] & samples_table_2['get_single_strand']
    ]
    choices_ss_83_or_i16_bw = [
        ("Results/" + samples_table_2['reference_genome'] +
         "/Bigwigs/Coverage/Single_strand/Full_length_reads/" + cov_params +
         samples_nodup_filt + ".83-163.bw"),
        "Results/" + samples_table_2['reference_genome'] +
        "/Bigwigs/Coverage/Single_strand/Full_length_reads/" + cov_params +
        samples_nodup_filt + ".inc_16.bw"
    ]
    
    samples_table_2['ss_83_or_i16_bw'] = (
        np.select(conditions_ss_83_or_i16_bw,
                  choices_ss_83_or_i16_bw,
                  default=np.NaN)
    )
    
    samples_table_2['ss_99_or_e16_bw'] = (
        np.where(samples_table_2['PE'],
                 "Results/" + samples_table_2['reference_genome'] +
                 "/Bigwigs/Coverage/Single_strand/Full_length_reads/" + cov_params +
                 samples_nodup_filt + ".99-147.bw",
                 "Results/" + samples_table_2['reference_genome'] +
                 "/Bigwigs/Coverage/Single_strand/Full_length_reads/" + cov_params +
                 samples_nodup_filt + ".exc_16.bw"
                 )
    )
    
    # Single strand bw, clipped profiles
    # PE data
    conditions_ss_not_clpd_bw = [
        samples_table_2['PE'] & samples_table_2['get_single_strand'],
        ~samples_table_2['PE'] & samples_table_2['get_single_strand']
    ]
    choices_ss_83_or_i16_bw = [
        ("Results/" + samples_table_2['reference_genome'] +
         "/Bigwigs/Coverage/Single_strand/Full_length_reads/" + cov_params +
         samples_nodup_filt + ".83-163.bw"),
        "Results/" + samples_table_2['reference_genome'] +
        "/Bigwigs/Coverage/Single_strand/Full_length_reads/" + cov_params +
        samples_nodup_filt + ".inc_16.bw"
    ]
    choices_ss_99_or_e16_bw = [
        ("Results/" + samples_table_2['reference_genome'] +
         "/Bigwigs/Coverage/Single_strand/Full_length_reads/" + cov_params +
         samples_nodup_filt + ".99-147.bw"),
        "Results/" + samples_table_2['reference_genome'] +
        "/Bigwigs/Coverage/Single_strand/Full_length_reads/" + cov_params +
        samples_nodup_filt + ".exc_16.bw"
    ]
    
    samples_table_2['ss_83_or_i16_bw'] = (
        np.select(conditions_ss_not_clpd_bw,
                  choices_ss_83_or_i16_bw,
                  default=np.NaN)
    )
    samples_table_2['ss_99_or_e16_bw'] = (
        np.select(conditions_ss_not_clpd_bw,
                  choices_ss_99_or_e16_bw,
                  default=np.NaN)
    )
    
    # Clipped
    samples_table_2['ssPE_83_i16_clpd_bw'] = np.where(
        samples_table_2['PE'] & samples_table_2['Clip_reads_to_1bp_on_5_prime'],
        ("Results/" + samples_table_2['reference_genome'] +
         "/Bigwigs/Coverage/Single_strand/1bp_clipped_reads/" + cov_params + samples_nodup_filt +
         ".83-163.inc_16.clipped_1_bp.bw"),
        np.NaN)
    
    samples_table_2['ssPE_83_e16_clpd_bw'] = np.where(
        samples_table_2['PE'] & samples_table_2['Clip_reads_to_1bp_on_5_prime'],
        ("Results/" + samples_table_2['reference_genome'] +
         "/Bigwigs/Coverage/Single_strand/1bp_clipped_reads/" + cov_params + samples_nodup_filt +
         ".83-163.exc_16.clipped_1_bp.bw"),
        np.NaN)
    
    samples_table_2['ssPE_99_i16_clpd_bw'] = np.where(
        samples_table_2['PE'] & samples_table_2['Clip_reads_to_1bp_on_5_prime'],
        ("Results/" + samples_table_2['reference_genome'] +
         "/Bigwigs/Coverage/Single_strand/1bp_clipped_reads/" + cov_params + samples_nodup_filt +
         ".99-147.inc_16.clipped_1_bp.bw"),
        np.NaN)
    
    
    samples_table_2['ssPE_99_e16_clpd_bw'] = np.where(
        samples_table_2['PE'] & samples_table_2['Clip_reads_to_1bp_on_5_prime'],
        ("Results/" + samples_table_2['reference_genome'] +
         "/Bigwigs/Coverage/Single_strand/1bp_clipped_reads/" + cov_params + samples_nodup_filt +
         ".99-147.exc_16.clipped_1_bp.bw"),
        np.NaN)
    
    # SR data
    samples_table_2['ssSR_i16_clp_bw'] = np.where(
        ~samples_table_2['PE'] & samples_table_2['Clip_reads_to_1bp_on_5_prime'],
        ("Results/" + samples_table_2['reference_genome'] +
         "/Bigwigs/Coverage/Single_strand/1bp_clipped_reads/" + cov_params + samples_nodup_filt +
         ".inc_16.clipped_1_bp.bw"),
        np.NaN)
    samples_table_2['ssSR_e16_clp_bw'] = np.where(
        ~samples_table_2['PE'] & samples_table_2['Clip_reads_to_1bp_on_5_prime'],
        ("Results/" + samples_table_2['reference_genome'] +
         "/Bigwigs/Coverage/Single_strand/1bp_clipped_reads/" + cov_params + samples_nodup_filt +
         ".exc_16.clipped_1_bp.bw"),
        np.NaN)
    
    # Peak files
    samples_table_2['peak_bl_gr_flt_nrw'] = np.where(
        samples_table['peak_ctrl_file_alias'] != "-",
        "Results/" + samples_table_2['reference_genome'] + "/Peaks/MACS2/Narrow/" +
        peak_nrw_params + "/Black-grey_filtered/" + samples_nodup_filt +
        ".narrowPeak",
        np.NaN
    )
    
    samples_table_2['peak_bl_gr_flt_brd'] = np.where(
        samples_table['peak_ctrl_file_alias'] != "-",
        "Results/" + samples_table_2['reference_genome'] + "/Peaks/MACS2/Broad/" +
        peak_brd_params + "/Black-grey_filtered/" + samples_nodup_filt +
        ".broadPeak",
        np.NaN
    )
    
    samples_table_2['peak_bl_gr_flt_hs_int_brd'] = np.where(
        (samples_table_2["reference_genome"] == "mm10") & (
            samples_table['peak_ctrl_file_alias'] != "-"),
        ("Results/" + samples_table_2['reference_genome'] + "/Peaks/MACS2/Broad/" +
         peak_brd_params + "/Black-grey_filtered/Intersect_HSs_plus_minus_2000_bp/" +
         samples_nodup_filt + ".broadPeak"),
        np.NaN)
    
    samples_table_2['peak_bl_gr_flt_hs_int_nrw'] = np.where(
        (samples_table_2["reference_genome"] == "mm10") & (
            samples_table['peak_ctrl_file_alias'] != "-"),
        ("Results/" + samples_table_2['reference_genome'] + "/Peaks/MACS2/Narrow/" +
         peak_nrw_params + "/Black-grey_filtered/Intersect_HSs_plus_minus_2000_bp/" +
         samples_nodup_filt +
         ".narrowPeak"),
        np.NaN)
    
    # Qctrl files
    # FRIP
    samples_table_2['FRIP_nrw_blk_gr_flt'] = np.where(
        samples_table['peak_ctrl_file_alias'] != "-",
        ("Results/" + samples_table_2['reference_genome'] +
         "/Qctrl/" + samples_table_2['sample_name'] + "/Processed_bam/FRIP/MACS2_Narrow_" +
         peak_nrw_params + "_bl-gr_flt/" + samples_nodup_filt + ".FRIP.txt"),
        np.NaN
    )
    samples_table_2['FRIP_brd_blk_gr_flt'] = np.where(
        samples_table['peak_ctrl_file_alias'] != "-",
        ("Results/" + samples_table_2['reference_genome'] +
         "/Qctrl/" + samples_table_2['sample_name'] + "/Processed_bam/FRIP/MACS2_Broad_" +
         peak_brd_params + "_bl-gr_flt/" + samples_nodup_filt + ".FRIP.txt"),
        np.NaN
    )
    # samtools flagstat
    samples_table_2['raw_flagstat'] = (
        "Results/" + samples_table_2['reference_genome'] +
        "/Qctrl/" + samples_table_2['sample_name'] + "/Raw_bam/" +
        samples_table_2['sample_name'] + ".flagstat.txt"
    )
    samples_table_2['processed_flagstat'] = (
        "Results/" + samples_table_2['reference_genome'] +
        "/Qctrl/" + samples_table_2['sample_name'] + "/Processed_bam/" +
        samples_table_2['sample_name'] + ".flagstat.txt"
    )
    samples_table_2['processed_flagstat_dros'] = np.where(
        samples_table["dros_spike_in"],
        ("Results/d6/Qctrl/" + samples_table_2['sample_name'] + "/Processed_bam/" +
         samples_table_2['sample_name'] + ".flagstat.txt"),
        np.NaN)
    
    
    # Merged files (same experiment sequenced more than once)
    samples_table_2['merged_dedup_bam'] = np.where(
        ~((samples_table_2['merge_with'] == "") | (
            samples_table_2['merge_with'] == "-")),
        "Results/" + samples_table_2['reference_genome'] +
        "/Bams/Both_strands/" + samples_table_2['sample_name'] + "_MERGED" + final_genomes +
        nodup_filt + ".bam",
        np.NaN)
    # Hotspot heatmaps (only for top 5000 hs, both strands)
    samples_table_2['heatmap_top_5000_hs'] = np.where(
        samples_table_2["top5000_HS_heatmap"],
        ("Results/" + samples_table_2['reference_genome'] +
         "/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/top_5000_plus_minus_2000/Both_strands/" +
         cov_params + "Heatmaps/" + samples_nodup_filt + ".png"),
        np.NaN)
    
    # Dros. normalized/equalized samples (bw)
    samples_table_2['dros_eq_both_strand_cov'] = (
        np.where(~samples_table['dros_equalization_group'].isna(),
                 "Results/" + samples_table['reference_genome'] +
                 "/Bigwigs/Coverage/Both_strands/" + cov_params +
                 "Drosophila_normalized/" + samples_table['dros_equalization_group'] +
                 "/" + samples_nodup_filt + ".dros_norm.bw",
                 np.NaN))
    
    conditions_dros_eq = [
        (~samples_table['dros_equalization_group'].isna()) & (
            samples_table_2['PE']) & (samples_table_2['get_single_strand']),
        (~samples_table['dros_equalization_group'].isna()) & (
            ~samples_table_2['PE']) & (samples_table_2['get_single_strand'])
    ]
    
    choices_ss_83_or_i16_cov = [
        "Results/" + samples_table_2['reference_genome'] +
        "/Bigwigs/Coverage/Single_strand/Full_length_reads/" + cov_params +
        "Drosophila_normalized/" + samples_table['dros_equalization_group'] + "/" +
        samples_nodup_filt + ".83-163.dros_norm.bw",
        "Results/" + samples_table_2['reference_genome'] +
        "/Bigwigs/Coverage/Single_strand/Full_length_reads/" + cov_params +
        "Drosophila_normalized/" + samples_table['dros_equalization_group'] + "/" +
        samples_nodup_filt + ".inc_16.dros_norm.bw"
    ]
    
    choices_ss_99_or_e16_cov = [
        "Results/" + samples_table_2['reference_genome'] +
        "/Bigwigs/Coverage/Single_strand/Full_length_reads/" + cov_params +
        "Drosophila_normalized/" + samples_table['dros_equalization_group'] + "/" +
        samples_nodup_filt + ".99-147.dros_norm.bw",
        "Results/" + samples_table_2['reference_genome'] +
        "/Bigwigs/Coverage/Single_strand/Full_length_reads/" + cov_params +
        "Drosophila_normalized/" + samples_table['dros_equalization_group'] + "/" +
        samples_nodup_filt + ".exc_16.dros_norm.bw"]
    
    samples_table_2['dros_eq_83_or_i16_cov'] = (
        np.select(conditions_dros_eq,
                  choices_ss_83_or_i16_bam,
                  default=np.NaN)
    )
    samples_table_2['dros_eq_99_or_e16_cov'] = (
        np.select(conditions_dros_eq,
                  choices_ss_99_or_e16_cov,
                  default=np.NaN)
    )
    return samples_table_2
