# Run script to set up te pipeline  
On the terminal run the following code (copy all together):  
`ml slurm python/3.10.2 pandas/1.4.2 && \
python /Volumes/Pezza/hpc-nobackup/Agustin/test_folder/ChIPseq_Analysis_Snakemake/workflow/Scripts/Make_samples_table_lab_library.py`

# Fill/check samples_table.csv  
`samples_table.csv` is found at folder `Config` and contains all the metadata related to your samples.  

## Check the following columns to see if content is correct:  
- **sample_name**: str
	It can be any string as long as it does not contain spaces, symbols (other than `-` or `_`) or ends-up with the word `_MERGED`.  
- **PE**: Boolean (`True` or `False`)
	Answers the question: 'was the sample sequenced as Paired-End (PE)?'.    
	If you have different kind of sequencing read-pairing for different samples you should manually fill this column.  
- **library_technology**: `regular` or `adaptase`
	Any library prepared by doing end-repair (overhanging DNA digestion + gap fill followed by dsDNA adaptor ligation) or using tagmentation should be set as `regular`. 
	`Adaptase` refers to the technique in which a low complexity DNA-tail is added to the start of read 2. For more information see this IDT's xGen™ ssDNA & Low-Input DNA Library Prep Kit (at the time of writing cat # 10009859). Procedurally speaking, the only difference in the pipeline is that, after running fastp module and bedore alignement, 10 bases will be prunned from the start of both read 1 and 2.  
- **reference_genome**: `mm10` or `hg38`
	Reference genome build to which align the samples.
- **peak_ctrl_file_alias**: string (optional)
	Alias for control file to be used to call peaks by `MACS2`. Available aliases are found on `Config/config.yaml` file, on `MACS2:control`.  
	If no peak is to be called, it should be left empty or with a dash (`-`). To call peaks with no control put `no_input`.  
- **dros_spike_in**: Boolean (`True` or `False`)
	Answers the question: 'does the sample have drosophila spike-in?'  
- **get_single_strand**: Boolean (`True` or `False`)
	Answers the question: 'should single strand profiles (bam, biwigs and plots) be created?'.  
	**NOTE:** samples with no ss profiles won't have heatmaps or aggregate profiles on any hotspot list other than top-5000 hotspots.  
- **Clip_reads_to_1bp_on_5_prime**: Boolean (`True` or `False`)
	Answers the question: 'should read be clipped on their 5' end to get 1 bp-long reads?'.  
	This is necessary for techniques such as [S1-seq](http://www.genesdev.org/cgi/doi/10.1101/gad.336032.119) or [END-seq](https://doi.org/10.1038/s41467-020-14654-w).  
- **top5000_HS_heatmap**: Boolean (`True` or `False`)
	Answers the question: 'should heatmaps and aggregate profiles at top 5000 hotspots be created?'.  
	**NOTE:** samples set to `False` will not have aggreagte profiles or heatmaps at any other list of hotspots either.  
- Size_DNA_top_5000_HS: Not functional yet.
- **merge_with**: str
	Samples with wich to merge bam files.  
	Space separated list of either other experiment names within the same samples table,  or paths to raw bam files (or combination of both). The name of the merged sample will be the same as the sample of the row you are filling, plus the suffix '_MERGED'. Hence, this filed needs to be filled only in the sample carrying the name that you want to keep for the merged sample (see example below). All the settings for the sample to be merged will be copied to the merged sample. Everything calculated for the original sample will also be calculated for the merged sample.  
- **dros_equalization_group**: str
	Name for the group of samples to be normalized by drosophila spike-in. For example, if you want to normalize `ChIP_X_WT` and `ChIP_X_KO` you could fill this field (on both samples) with `ChIP_X_WT_vs_KO`.  
- **B6xCAST**: Boolean (`True` or `False`)
	Answers the question: 'Was this experiment done using F1 of a cross between B6 and CAST?'.  