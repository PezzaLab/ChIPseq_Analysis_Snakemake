#!/usr/bin/env python3
"""
Unit test suite for multi-genome support per sample in ChIP-seq Snakemake pipeline.
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

def test_single_genome_backwards_compatible():
    """Test that existing single-genome rows expand identically and preserve all attributes."""
    csv_data = """sample_name,fastq1,fastq2,PE,library_technology,reference_genome,peak_ctrl_file_alias,dros_spike_in,get_single_strand,Clip_reads_to_1bp_on_5_prime,top5000_HS_heatmap,Size_DNA_top_5000_HS,merge_with,B6xCAST
SampleA,R1.fq.gz,R2.fq.gz,True,regular,mm10,no_input,False,True,False,True,False,-,False
SampleB,R1_se.fq.gz,-,False,regular,mm39,no_input,True,True,False,True,False,-,False
"""
    raw_df = pd.read_csv(io.StringIO(csv_data), dtype=str)
    expanded = smkf.expand_sample_table(raw_df)

    assert len(expanded) == 2
    assert "SampleA.mm10" in expanded.index
    assert "SampleB.mm39" in expanded.index

    assert expanded.loc["SampleA.mm10", "PE"] == True
    assert expanded.loc["SampleA.mm10", "dros_spike_in"] == False
    assert expanded.loc["SampleA.mm10", "reference_genome"] == "mm10"
    assert pd.isna(expanded.loc["SampleA.mm10", "merge_with"])

    assert expanded.loc["SampleB.mm39", "PE"] == False
    assert expanded.loc["SampleB.mm39", "dros_spike_in"] == True
    assert expanded.loc["SampleB.mm39", "reference_genome"] == "mm39"
    assert pd.isna(expanded.loc["SampleB.mm39", "fastq2"])
    print("PASS: test_single_genome_backwards_compatible")


def test_multi_genome_broadcasting():
    """Test that a sample with multiple genomes broadcasts single values to all genomes."""
    csv_data = """sample_name,fastq1,fastq2,PE,library_technology,reference_genome,peak_ctrl_file_alias,dros_spike_in,get_single_strand,Clip_reads_to_1bp_on_5_prime,top5000_HS_heatmap,Size_DNA_top_5000_HS,merge_with,B6xCAST
MultiSamp,R1.fq.gz,R2.fq.gz,True,regular,"mm10;mm39",no_input,False,True,False,True,False,-,False
"""
    raw_df = pd.read_csv(io.StringIO(csv_data), dtype=str)
    expanded = smkf.expand_sample_table(raw_df)

    assert len(expanded) == 2
    assert "MultiSamp.mm10" in expanded.index
    assert "MultiSamp.mm39" in expanded.index

    # Both genomes should share sample-level and broadcasted attributes
    for g in ["mm10", "mm39"]:
        idx = f"MultiSamp.{g}"
        assert expanded.loc[idx, "sample_name"] == "MultiSamp"
        assert expanded.loc[idx, "reference_genome"] == g
        assert expanded.loc[idx, "fastq1"] == "R1.fq.gz"
        assert expanded.loc[idx, "fastq2"] == "R2.fq.gz"
        assert expanded.loc[idx, "PE"] == True
        assert expanded.loc[idx, "library_technology"] == "regular"
        assert expanded.loc[idx, "peak_ctrl_file_alias"] == "no_input"
        assert expanded.loc[idx, "dros_spike_in"] == False
        assert expanded.loc[idx, "get_single_strand"] == True
        assert expanded.loc[idx, "Clip_reads_to_1bp_on_5_prime"] == False
        assert expanded.loc[idx, "top5000_HS_heatmap"] == True
        assert pd.isna(expanded.loc[idx, "merge_with"])
    print("PASS: test_multi_genome_broadcasting")


def test_multi_genome_explicit_values():
    """Test specifying per-genome values (1-to-1 mapping)."""
    csv_data = """sample_name,fastq1,fastq2,PE,library_technology,reference_genome,peak_ctrl_file_alias,dros_spike_in,get_single_strand,Clip_reads_to_1bp_on_5_prime,top5000_HS_heatmap,Size_DNA_top_5000_HS,merge_with,B6xCAST
MultiSamp,R1.fq.gz,R2.fq.gz,True,regular,"mm10;mm39","ctrl_mm10;ctrl_mm39",False,True,False,"True;False",False,"mergeA;mergeB",False
"""
    raw_df = pd.read_csv(io.StringIO(csv_data), dtype=str)
    expanded = smkf.expand_sample_table(raw_df)

    assert len(expanded) == 2
    # mm10
    assert expanded.loc["MultiSamp.mm10", "peak_ctrl_file_alias"] == "ctrl_mm10"
    assert expanded.loc["MultiSamp.mm10", "top5000_HS_heatmap"] == True
    assert expanded.loc["MultiSamp.mm10", "merge_with"] == "mergeA"

    # mm39
    assert expanded.loc["MultiSamp.mm39", "peak_ctrl_file_alias"] == "ctrl_mm39"
    assert expanded.loc["MultiSamp.mm39", "top5000_HS_heatmap"] == False
    assert expanded.loc["MultiSamp.mm39", "merge_with"] == "mergeB"
    print("PASS: test_multi_genome_explicit_values")


def test_comma_and_semicolon_delimiters():
    """Test both commas and semicolons as delimiters."""
    csv_data = """sample_name,fastq1,fastq2,PE,library_technology,reference_genome,peak_ctrl_file_alias,dros_spike_in,get_single_strand,Clip_reads_to_1bp_on_5_prime,top5000_HS_heatmap,Size_DNA_top_5000_HS,merge_with,B6xCAST
SampComma,R1.fq.gz,R2.fq.gz,True,regular,"mm10, mm39","c1, c2",False,True,False,"True, False",False,-,False
"""
    raw_df = pd.read_csv(io.StringIO(csv_data), dtype=str)
    expanded = smkf.expand_sample_table(raw_df)

    assert len(expanded) == 2
    assert expanded.loc["SampComma.mm10", "peak_ctrl_file_alias"] == "c1"
    assert expanded.loc["SampComma.mm10", "top5000_HS_heatmap"] == True
    assert expanded.loc["SampComma.mm39", "peak_ctrl_file_alias"] == "c2"
    assert expanded.loc["SampComma.mm39", "top5000_HS_heatmap"] == False
    print("PASS: test_comma_and_semicolon_delimiters")


def test_mismatched_value_count_raises_error():
    """Test that providing M values where M != 1 and M != N raises an informative error."""
    csv_data = """sample_name,fastq1,fastq2,PE,library_technology,reference_genome,peak_ctrl_file_alias,dros_spike_in,get_single_strand,Clip_reads_to_1bp_on_5_prime,top5000_HS_heatmap,Size_DNA_top_5000_HS,merge_with,B6xCAST
ErrSamp,R1.fq.gz,R2.fq.gz,True,regular,"mm10;mm39","ctrl1;ctrl2;ctrl3",False,True,False,True,False,-,False
"""
    raw_df = pd.read_csv(io.StringIO(csv_data), dtype=str)
    try:
        smkf.expand_sample_table(raw_df)
        assert False, "Should have exited or raised error on mismatched value count"
    except SystemExit as e:
        assert "has 3 values" in str(e)
        print("PASS: test_mismatched_value_count_raises_error")


def test_generate_samples_table_2_multi_genome():
    """Test suffixes.generate_samples_table_2 with multi-genome expanded table."""
    config_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../Config/config.yaml"))
    with open(config_path) as f:
        config = yaml.safe_load(f)

    csv_data = """sample_name,fastq1,fastq2,PE,library_technology,reference_genome,peak_ctrl_file_alias,dros_spike_in,get_single_strand,Clip_reads_to_1bp_on_5_prime,top5000_HS_heatmap,Size_DNA_top_5000_HS,merge_with,B6xCAST
Bloom_10_PE,R1.fq.gz,R2.fq.gz,True,regular,"mm10;mm39",no_input,False,True,False,True,False,-,False
"""
    raw_df = pd.read_csv(io.StringIO(csv_data), dtype=str)
    expanded = smkf.expand_sample_table(raw_df)

    st2 = sfxs.generate_samples_table_2(expanded, config)
    assert len(st2) == 2

    # Check mm10 paths
    row_mm10 = st2.loc["Bloom_10_PE.mm10"]
    assert row_mm10["reference_genome"] == "mm10"
    assert row_mm10["raw_bam"] == "Results/Bloom_10_PE.mm10.bam"
    assert row_mm10["dedup_flt_both_strds_bam"] == "Results/mm10/Bams/Both_strands/Bloom_10_PE.mm10.q_filt.srt.nodup.mit_filt.bam"
    assert "Results/mm10/Bigwigs/Coverage/Both_strands/" in row_mm10["Both_strands_coverage_bw"]
    assert "Results/mm10/Analysis/Heatmaps_and_aggregate_profiles/" in row_mm10["heatmap_top_5000_hs"]

    # Check mm39 paths
    row_mm39 = st2.loc["Bloom_10_PE.mm39"]
    assert row_mm39["reference_genome"] == "mm39"
    assert row_mm39["raw_bam"] == "Results/Bloom_10_PE.mm39.bam"
    assert row_mm39["dedup_flt_both_strds_bam"] == "Results/mm39/Bams/Both_strands/Bloom_10_PE.mm39.q_filt.srt.nodup.mit_filt.bam"
    assert "Results/mm39/Bigwigs/Coverage/Both_strands/" in row_mm39["Both_strands_coverage_bw"]
    assert "Results/mm39/Analysis/Heatmaps_and_aggregate_profiles/" in row_mm39["heatmap_top_5000_hs"]

    # Test helper functions
    smkf.samples_table = expanded
    smkf.samples_table_2 = st2
    smkf.config = config

    bam_mm10 = smkf.get_processed_bam("Bloom_10_PE", "mm10")
    bam_mm39 = smkf.get_processed_bam("Bloom_10_PE", "mm39")
    assert "Results/mm10/" in bam_mm10
    assert "Results/mm39/" in bam_mm39

    raw_mm10 = smkf.get_raw_bam("Bloom_10_PE", "mm10")
    raw_mm39 = smkf.get_raw_bam("Bloom_10_PE", "mm39")
    assert raw_mm10 == "Results/Bloom_10_PE.mm10.bam"
    assert raw_mm39 == "Results/Bloom_10_PE.mm39.bam"

    fq1 = smkf.get_sample_fastq1("Bloom_10_PE")
    assert fq1 == "R1.fq.gz"

    print("PASS: test_generate_samples_table_2_multi_genome")


def test_adaptase_multi_genome():
    """Test adaptase trimming flow with multi-genome samples."""
    csv_data = """sample_name,fastq1,fastq2,PE,library_technology,reference_genome,peak_ctrl_file_alias,dros_spike_in,get_single_strand,Clip_reads_to_1bp_on_5_prime,top5000_HS_heatmap,Size_DNA_top_5000_HS,merge_with,B6xCAST
AdaptaseSamp,R1.fq.gz,R2.fq.gz,True,adaptase,"mm10;mm39",no_input,False,True,False,True,False,-,False
"""
    raw_df = pd.read_csv(io.StringIO(csv_data), dtype=str)
    expanded = smkf.expand_sample_table(raw_df)

    assert len(expanded) == 2
    assert expanded.loc["AdaptaseSamp.mm10", "library_technology"] == "adaptase"
    assert expanded.loc["AdaptaseSamp.mm39", "library_technology"] == "adaptase"

    # Conflicting library_technology for the same sample must raise error
    conflicting_csv = """sample_name,fastq1,fastq2,PE,library_technology,reference_genome,peak_ctrl_file_alias,dros_spike_in,get_single_strand,Clip_reads_to_1bp_on_5_prime,top5000_HS_heatmap,Size_DNA_top_5000_HS,merge_with,B6xCAST
ConfSamp,R1.fq.gz,R2.fq.gz,True,"adaptase;regular","mm10;mm39",no_input,False,True,False,True,False,-,False
"""
    raw_conf = pd.read_csv(io.StringIO(conflicting_csv), dtype=str)
    try:
        smkf.expand_sample_table(raw_conf)
        assert False, "Should have raised error on conflicting library_technology"
    except SystemExit as e:
        assert "sample-intrinsic column 'library_technology' has conflicting values" in str(e)

    # Test align_fastq_input returns .adaptase_trimmed FASTQ paths
    class MockWildcards:
        sample = "AdaptaseSamp"
        genomes_all = "mm10"

    smkf.samples_table = expanded
    # Mock config references for mm10
    config_path = os.path.abspath(os.path.join(os.path.dirname(__file__), "../../Config/config.yaml"))
    with open(config_path) as f:
        config = yaml.safe_load(f)
    smkf.config = config

    # We temporarily inject a dummy directory to satisfy os.listdir in align_fastq_input
    import tempfile
    with tempfile.TemporaryDirectory() as tmpdir:
        for ext in [".sa", ".pac", ".bwt", ".ann", ".amb"]:
            open(os.path.join(tmpdir, f"ref.fa{ext}"), "w").close()
        smkf.config['genomes']['mm10'] = tmpdir

        inp = smkf.align_fastq_input(MockWildcards())
        assert "AdaptaseSamp.adap_trimmed.adaptase_trimmed.R1.PE.fq.gz" in inp["fq1"]
        assert "AdaptaseSamp.adap_trimmed.adaptase_trimmed.R2.PE.fq.gz" in inp["fq2"]

    print("PASS: test_adaptase_multi_genome")


if __name__ == "__main__":
    test_single_genome_backwards_compatible()
    test_multi_genome_broadcasting()
    test_multi_genome_explicit_values()
    test_comma_and_semicolon_delimiters()
    test_mismatched_value_count_raises_error()
    test_generate_samples_table_2_multi_genome()
    test_adaptase_multi_genome()
    print("\nALL UNIT TESTS PASSED SUCCESSFULLY!")

