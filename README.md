# ChIPseq_Analysis_Snakemake
Snakemake pipeline to analyze ChIP-seq (and cut and run) experiments, with a focus on meiotic processes.  

## Pipeline output  
Pipeline has the capability of generating the following files:
- Deduplicated and q-filtered (dqf) bam files
- Single strand dqf bam files
- 1 bp front-clipped-reads dqf bam files (for experiments such as END-seq or S1-seq
- Bigwig files of all menioned bams
- Quality control metrics:
  * FASTQ metrics calculated by `fastp`
  * Crosscorelation metrics (NSC and RSC) calculated by `phantompeakqualtools` 
  * FRIP (fraction of reads in peaks)
  * General metrics calculated by `samtools stats` and `samtools flagstat`
  * Insert size calculated by `picard CollectInsertSizeMetrics`
  * Library complexity calculated by `picard EstimateLibraryComplexity`
  * Reads duplication status calculated by `picard MarkDuplicates`
- Quality control report produced by `multiQC`
- Narrow and broad peaks (MACS2)
- Above mentioned peaks filtered by black(+grey) list
- Above mentioned peaks annotated by `homer`
- Above mentioned peaks intersected with SPO11 hotspots (mice only)
- Heatmaps and aggregate profiles in the following regions:
  * Top 5000 strongest hotspots (as per Spo11 signal, mice only)
  * X non-PAR hotspots and Spo11-intensity matched autosomal hotspots (mice only)
  * Asymetric hotspots (mice only)
- Library html report.

The pipeline can also merge FASTQ files and generate (drosophila) spike-in normalized bigwig files.

## Running the pipeline  

### With our libraries (add link here)  

### With mined-data (add link here)
