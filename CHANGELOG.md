#Changelog
#Notable changes made to the pipeline will be recorded in this file 

## [Unreleased]
- Placeholder for  upcoming features, bug fixes or improvement. 

## [V1.0.0] - 2025-07-29
- First release of the modular ATAC-seq Pipeline
- Core Processes:
  - Preprocessing: fastp, FastQC
  - Alignment: Bowtie2, Samtools sorting, indexing, and deduplication
  - Post-Alignment QC: mitochondrial read quantification, fragment size analysis, TSS enrichment, PhantomPeakQualTools, Preseq, Qualimap
  - Coverage and normalization: genome coverage, bigwig conversion, CPM normalization
  - Peak Calling and Filtering: MACS2, ENCODE blacklist filtering 
  - Visualization: heatmaps, motif analysis, correlation plots
- Design Features: 
 - Modular Snakefile  with per-rule configuration
 - Configurable via `config.yaml`
 - Optimized for scalability and reproducibility across datasets
 
## [V1.0.1] - YYYY-MM-DD
 - Placeholder for  updates made to the pipeline such as bug fixes, performance improvement, etc.  


