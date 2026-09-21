#!/usr/bin/env python3
"""
Unit test suite for running the ChIP-seq Snakemake pipeline without peak calling / FRIP.
"""
import os
import sys
import io
import pandas as pd
import numpy as np
import yaml

# Add workflow directory to sys.path
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), "..")))

import smk_functions as smkf
import suffixes as sfxs

CONFIG = {
    "qctrl": True,
    "library": {"name": "Test_lib"},
    "min_mapping_qual": "30",
    "coverage": {
        "bin_size": "4",
        "smooth": "1",
        "extend_reads": "250",
        "normalization": "CPM",
    },
    "MACS2": {
        "gsize": {"hs": 2.7e9, "mm": 1.87e9},
        "broad_cutoff": "0.1",
        "qvalue": "0.05",
        "control": {
            "no_input": "",
            "mock_ctrl": "/path/to/ctrl.bam",
        },
        "extension": "300",
    },
    "heatmaps": {
        "before_region": "3000",
        "after_region": "3000",
        "bin_size": "2",
    },
    "aggregate_profiles": {
        "process_outfile": {
            "overwrite": "TRUE",
            "smooth": True,
        }
    },
}


def test_no_peaks_all_samples():
    """Test that setting peak_ctrl_file_alias to '-' for all samples defines peak & FRIP columns as NaN without KeyError."""
    csv_data = """sample_name,fastq1,fastq2,PE,library_technology,reference_genome,peak_ctrl_file_alias,dros_spike_in,get_single_strand,Clip_reads_to_1bp_on_5_prime,top5000_HS_heatmap,Size_DNA_top_5000_HS,merge_with,B6xCAST
SampleA,R1.fq.gz,R2.fq.gz,True,regular,mm10,-,False,True,False,True,False,-,False
SampleB,R1_se.fq.gz,-,False,regular,mm39,-,False,True,False,True,False,-,False
"""
    raw_df = pd.read_csv(
        io.StringIO(csv_data),
        dtype=str,
        na_values={"merge_with": "-", "fastq2": "-", "peak_ctrl_file_alias": "-"},
        keep_default_na=True,
    )
    expanded = smkf.expand_sample_table(raw_df)
    s2 = sfxs.generate_samples_table_2(expanded, CONFIG)

    # All peak and FRIP columns must exist in samples_table_2
    for col in [
        "narrow_peak_bl_gr_flt",
        "broad_peak_bl_gr_flt",
        "narrow_peak_bl_gr_flt_hs_int",
        "broad_peak_bl_gr_flt_hs_int",
        "narrow_blk_gr_flt_FRIP",
        "broad_blk_gr_flt_FRIP",
    ]:
        assert col in s2.columns, f"Column {col} missing from samples_table_2"
        assert s2[col].isna().all(), f"Expected column {col} to be all NaN, got non-null values"

    # broad_frip_scores must evaluate cleanly to empty list
    broad_frip_scores = [
        f for f in s2["broad_blk_gr_flt_FRIP"].values.tolist()
        if pd.notna(f)
    ]
    assert broad_frip_scores == []

    # basic_peaks must evaluate cleanly to empty list
    basic_peaks = [
        p for p in s2.filter(regex=r"^(narrow|broad)_peak_bl_gr_flt$").values.flatten().tolist()
        if pd.notna(p)
    ] + broad_frip_scores
    assert basic_peaks == []

    # peaks_summary_files must evaluate to empty list
    peaks_summary_files = [
        f"Results/{genome}/Analysis/Peaks_summary.tsv"
        for genome in expanded["reference_genome"].unique()
        if (
            pd.notna(expanded.loc[expanded["reference_genome"] == genome, "peak_ctrl_file_alias"]) &
            (expanded.loc[expanded["reference_genome"] == genome, "peak_ctrl_file_alias"] != "-")
        ).any()
    ]
    assert peaks_summary_files == []


def test_mixed_peaks():
    """Test that when only some samples call peaks, only those samples receive peak and FRIP paths."""
    csv_data = """sample_name,fastq1,fastq2,PE,library_technology,reference_genome,peak_ctrl_file_alias,dros_spike_in,get_single_strand,Clip_reads_to_1bp_on_5_prime,top5000_HS_heatmap,Size_DNA_top_5000_HS,merge_with,B6xCAST
SampleWithCtrl,R1.fq.gz,R2.fq.gz,True,regular,mm10,no_input,False,True,False,True,False,-,False
SampleNoCtrl,R1.fq.gz,R2.fq.gz,True,regular,mm10,-,False,True,False,True,False,-,False
"""
    raw_df = pd.read_csv(
        io.StringIO(csv_data),
        dtype=str,
        na_values={"merge_with": "-", "fastq2": "-", "peak_ctrl_file_alias": "-"},
        keep_default_na=True,
    )
    expanded = smkf.expand_sample_table(raw_df)
    s2 = sfxs.generate_samples_table_2(expanded, CONFIG)

    ctrl_idx = "SampleWithCtrl.mm10"
    no_ctrl_idx = "SampleNoCtrl.mm10"

    # SampleWithCtrl should have peak and FRIP paths
    assert pd.notna(s2.loc[ctrl_idx, "narrow_peak_bl_gr_flt"])
    assert pd.notna(s2.loc[ctrl_idx, "broad_peak_bl_gr_flt"])
    assert pd.notna(s2.loc[ctrl_idx, "broad_blk_gr_flt_FRIP"])
    assert "no_input" in s2.loc[ctrl_idx, "broad_blk_gr_flt_FRIP"]

    # SampleNoCtrl should NOT have peak or FRIP paths
    assert pd.isna(s2.loc[no_ctrl_idx, "narrow_peak_bl_gr_flt"])
    assert pd.isna(s2.loc[no_ctrl_idx, "broad_peak_bl_gr_flt"])
    assert pd.isna(s2.loc[no_ctrl_idx, "broad_blk_gr_flt_FRIP"])

    # Collected FRIP scores should only contain 1 element
    broad_frip_scores = [
        f for f in s2["broad_blk_gr_flt_FRIP"].values.tolist()
        if pd.notna(f)
    ]
    assert len(broad_frip_scores) == 1
    assert "SampleWithCtrl" in broad_frip_scores[0]


def test_markdown_report_input_without_peaks():
    """Test that smkf.markdown_report_aggregate_profiles_input yields empty peaks_summary if no peaks are called."""
    class MockWildcard:
        genomes_not_fused = "mm10"
        smooth = "smoothed"

    csv_data = """sample_name,fastq1,fastq2,PE,library_technology,reference_genome,peak_ctrl_file_alias,dros_spike_in,get_single_strand,Clip_reads_to_1bp_on_5_prime,top5000_HS_heatmap,Size_DNA_top_5000_HS,merge_with,B6xCAST
SampleA,R1.fq.gz,R2.fq.gz,True,regular,mm10,-,False,True,False,True,False,-,False
"""
    raw_df = pd.read_csv(
        io.StringIO(csv_data),
        dtype=str,
        na_values={"merge_with": "-", "fastq2": "-", "peak_ctrl_file_alias": "-"},
        keep_default_na=True,
    )
    expanded = smkf.expand_sample_table(raw_df)
    s2 = sfxs.generate_samples_table_2(expanded, CONFIG)

    smkf.samples_table = expanded
    smkf.samples_table_2 = s2
    smkf.samples_table_no_merged_samples = expanded
    smkf.config = CONFIG

    w = MockWildcard()
    inp = smkf.markdown_report_aggregate_profiles_input(w)
    assert inp["peaks_summary"] == []


def test_summarize_peak_count_input_without_peaks():
    """Test that smkf.summarize_peak_count_input returns empty peak lists when no peaks are present."""
    class MockWildcard:
        genomes_not_fused = "mm10"

    csv_data = """sample_name,fastq1,fastq2,PE,library_technology,reference_genome,peak_ctrl_file_alias,dros_spike_in,get_single_strand,Clip_reads_to_1bp_on_5_prime,top5000_HS_heatmap,Size_DNA_top_5000_HS,merge_with,B6xCAST
SampleA,R1.fq.gz,R2.fq.gz,True,regular,mm10,-,False,True,False,True,False,-,False
"""
    raw_df = pd.read_csv(
        io.StringIO(csv_data),
        dtype=str,
        na_values={"merge_with": "-", "fastq2": "-", "peak_ctrl_file_alias": "-"},
        keep_default_na=True,
    )
    expanded = smkf.expand_sample_table(raw_df)
    s2 = sfxs.generate_samples_table_2(expanded, CONFIG)

    smkf.samples_table_2 = s2
    smkf.include_hotspots = True

    w = MockWildcard()
    peaks_dict = smkf.summarize_peak_count_input(w)
    assert peaks_dict["narrow_all"] == []
    assert peaks_dict["broad_all"] == []
    assert peaks_dict["narrow_hs"] == []
    assert peaks_dict["broad_hs"] == []
