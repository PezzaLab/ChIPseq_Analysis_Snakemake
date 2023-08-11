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

    samples_table_2 = samples_table.copy()
    
    ### Suffixes to make file names ###
    nodup_filt = ".q_filt.srt.nodup.mit_filt"
    
    
    aligned_genomes = (
        np.where(samples_table_2['dros_spike_in'],
                 "." + samples_table_2['reference_genome'] + "_f_d6",
                 "." + samples_table_2['reference_genome'])
    )
    
    final_genomes = (
        np.where(
            samples_table_2['dros_spike_in'],
            aligned_genomes + "." + samples_table_2['reference_genome'],
            aligned_genomes)
    )
    
    final_genomes_b6xcast = []
    for i in samples_table_2.index:
        if samples_table_2.loc[i, "dros_spike_in"]:
            final_genomes_b6xcast.append(".mm10_x_CAST_EiJ_f_d6.mm10_x_CAST_EiJ")
        else:
            final_genomes_b6xcast.append(".mm10_x_CAST_EiJ")
            # For some reason the np.where is not working if I don't concatenate 
            # an object with a string
    
    samples_nodup_filt = (
        samples_table_2['sample_name'] + final_genomes +
        nodup_filt
    )
    
    samples_nodup_filt_b6xcast = (
        samples_table_2['sample_name'] + final_genomes_b6xcast +
        nodup_filt
    )
    
    cov_params = f"{config['coverage']['normalization']}_bs{config['coverage']['bin_size']}_sm{config['coverage']['smooth']}/"
    dros_norm_cov_config_params_string=f"drosNormalized_bs{config['coverage']['bin_size']}_sm{config['coverage']['smooth']}/"
    
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
    
    # Peak files
    if pd.notna(samples_table['peak_ctrl_file_alias']).any(): 
       peak_types=["narrow","broad"]
       peak_params = (
           "qv_" + (config['MACS2']['qvalue'].split("."))[1] +
           "__" + samples_table_2['peak_ctrl_file_alias']
       )
       
       for peak_type in peak_types:
               # Strings

               # Modify default values accordingly
               if peak_type == "broad":
                   peak_params = (
                       "bco_" + (config['MACS2']['broad_cutoff'].split("."))[1] + "_qv_" +
                       (config['MACS2']['qvalue'].split("."))[1] + "__" +
                       samples_table_2['peak_ctrl_file_alias']
                   )
               
               # Generate column name and content   
               samples_table_2[f'{peak_type}_peak_bl_gr_flt'] = np.where(
                   samples_table['peak_ctrl_file_alias'] != "-",
                   "Results/" + samples_table_2['reference_genome'] + 
                   f"/Peaks/MACS2/{peak_type}/" + peak_params + 
                   "/Black-grey_filtered/" + samples_nodup_filt + f".{peak_type}Peak",
                   np.NaN
               )
               samples_table_2[f'{peak_type}_peak_bl_gr_flt_hs_int'] = np.where(
                   (samples_table_2["reference_genome"] == "mm10") & (
                       samples_table['peak_ctrl_file_alias'] != "-"),
                   "Results/" + samples_table_2['reference_genome'] + 
                   f"/Peaks/MACS2/{peak_type}/" + peak_params + 
                   "/Black-grey_filtered/Intersect_HSs_plus_minus_2000_bp/" +
                    samples_nodup_filt + f".{peak_type}Peak", 
               np.NaN)

               # FRIP
               samples_table_2[f'{peak_type}_blk_gr_flt_FRIP'] = np.where(
                   samples_table['peak_ctrl_file_alias'] != "-",
                   "Results/" + samples_table_2['reference_genome'] +
                   "/Qctrl/" + samples_table_2['sample_name'] + 
                   f"/Processed_bam/FRIP/MACS2_{peak_type}_" +
                   peak_params + "_bl-gr_flt/" + samples_nodup_filt + ".FRIP.txt",
                   np.NaN
               )
               
               # Annotated peaks
               samples_table_2[f'{peak_type}Peak_blk_gr_flt_annotated'] = np.where(
                   samples_table['peak_ctrl_file_alias'] != "-",
                   "Results/" + samples_table_2['reference_genome'] +
                   f"/Peaks/MACS2/{peak_type}/" + peak_params + 
                   "/Black-grey_filtered/Annotated_peaks/" +
                   samples_nodup_filt + f".{peak_type}Peak.annotated.tsv", 
                   np.NaN)
        
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
    
    # Hotspot heatmaps (only for top 5000 hs, both strands)
    samples_table_2['heatmap_top_5000_hs'] = np.where(
        samples_table_2["top5000_HS_heatmap"] & ~samples_table_2["B6xCAST"],
        ("Results/" + samples_table_2['reference_genome'] +
         "/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/top_5000_plus_minus_2000/Both_strands/" +
         cov_params + "Heatmaps/" + samples_nodup_filt + ".png"),
        np.NaN)
    samples_table_2['heatmap_B6xCAST_top_5000_hs'] = np.where(
        samples_table_2["top5000_HS_heatmap"] & samples_table_2["B6xCAST"],
        ("Results/" + samples_table_2['reference_genome'] +
         "/Analysis/Heatmaps_and_aggregate_profiles/Hotspots/B6xCAST_top_5000_pm_1000bp/Both_strands/" +
         cov_params + "Heatmaps/" + samples_nodup_filt + ".png"),
        np.NaN)
   
    ##########
    # Coverage MATRIX files (generated by deeptools 
    # ousing a bw and a bed file a list of sites of interest such as hotspots)
    ##########
    hotspots=["top_5000_plus_minus_2000", "x_non_par", 
                    "autosomal_x_non_par_ctrl", "asymetric_watson_strong", 
                    "asymetric_crick_strong", 
                    "B6xCAST_top_5000_pm_1000bp", 
                    "B6xCAST_PRDM9_assymetric_hs_invading_strand",
                    "B6xCAST_PRDM9_assymetric_hs_receiving_strand"]
    strands=["83-163", "99-147", "inc_16", "exc_16", "Both_strands"] 

    for hs in hotspots:
        for strand in strands: 
            # Strings
            strand_string="Both_strands/"
            genome=samples_table["reference_genome"]
            strand2="." + strand
            samples_nodup_filt2=samples_nodup_filt
            if strand == "Both_strands":
                strand2="" # "Both_strands" is not included in the name of the file
            
            # Booleans
            PE=True
            single_strand=True
            b6xcast=~samples_table_2["B6xCAST"]
            # Modify default values accordingly
            if re.search("B6xCAST", hs):
                b6xcast=samples_table_2["B6xCAST"]
            
            if re.match("83|99|inc_16|exc_16", strand):
                strand_string="Single_strand/Full_length_reads/"
                single_strand=samples_table_2["get_single_strand"]

            if re.match("83|99", strand):
                PE=samples_table_2["PE"]
            elif re.match("inc_16|exc_16", strand):
                PE=~samples_table_2["PE"]

            if re.match("B6xCAST_PRDM9_assymetric_hs", hs):
                genome="mm10_x_CAST_EiJ"
                samples_nodup_filt2=samples_nodup_filt_b6xcast
                
            # Generate column name and content   
            samples_table_2[f"{hs}_{strand}_matrix"]= np.where(
                samples_table_2["top5000_HS_heatmap"] & b6xcast & PE & single_strand,
                ("Results/" + genome + "/Analysis/" +
                  "Heatmaps_and_aggregate_profiles/Hotspots/" + hs + "/" +
                  f"{strand_string}{cov_params}matrixes/" + 
                  samples_nodup_filt2 + strand2 + ".matrix"),
                np.NaN)
    
    return samples_table_2
