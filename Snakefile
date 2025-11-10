
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>              Modular ATACseq Pipleine                                                                         #>
#>              Author: Himanshu Bhandary          
#>              Mail: 2032ushimanshu@gmail.com                                                              #>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#restartable: True

configfile: "config.yaml"

SAMPLES = config["samples"]

if not SAMPLES:
   raise ValueError("No samples found in config.yaml")


#Include all rule files
include: "rules/fastp_01.smk"
include: "rules/fastqc_02.smk"
include: "rules/bowtie2_03.smk"
include: "rules/samtools_sort_04.smk"
include: "rules/calculate_mito_reads_05.smk"
include: "rules/remove_mito_reads_06.smk"
include: "rules/samtools_index_07.smk"
include: "rules/samtools_fixmate_pre_08_rule.smk"
include: "rules/samtools_markdup_08.smk"
include: "rules/samtools_index_after_markdup_09.smk"
include: "rules/samtools_view_10.smk"
include: "rules/samtools_index_post_filter_11.smk"
include: "rules/tn5_shift_12.smk"
include: "rules/samtools_stats_13.smk"
include: "rules/fragment_size_analysis_14.smk"
include: "rules/picard_CollectAlignmentSummaryMetrics_15.smk"
include: "rules/picard_CollectInsertSizeMetrics_16.smk"
include: "rules/tss_enrichment_17.smk"
include: "rules/bedtools_genomecov_18.smk"
include: "rules/sorted_bedgraph_19.smk"
include: "rules/bigwig_conversion_20.smk"
include: "rules/correlation_analysis_21.smk"
include: "rules/normalize_coverage_22.smk"
include: "rules/macs2_peak_calling_23.smk"
include: "rules/blacklist_region_filter_24.smk"
include: "rules/heatmap_25.smk"
include: "rules/frip_calculation_26.smk"
include: "rules/peak_annotation_27.smk"
include: "rules/motif_analysis_28.smk"
#include: "rules/phantompeakqualtool_29.smk"
include: "rules/preseq_30.smk"
include: "rules/qualimap_bamqc_31.smk"
include: "rules/multiqc_32.smk"

rule all:
    input:
        [
        #=============================================================================================================
        
        # 01_preprocessing

        #=============================================================================================================
        # fastp_01 outputs
        expand("results/fastp/{sample}_R1_trimmed.fastq.gz", sample=SAMPLES), 
        expand("results/fastp/{sample}_R2_trimmed.fastq.gz", sample=SAMPLES), 
        expand("results/fastp/{sample}.html", sample=SAMPLES), 
        expand("results/fastp/{sample}.json", sample=SAMPLES),
        
        # fastqc_02 outputs
        expand("results/fastqc/{sample}_R1_trimmed_fastqc.html", sample=SAMPLES), 
        expand("results/fastqc/{sample}_R1_trimmed_fastqc.zip", sample=SAMPLES), 
        expand("results/fastqc/{sample}_R2_trimmed_fastqc.html", sample=SAMPLES), 
        expand("results/fastqc/{sample}_R2_trimmed_fastqc.zip", sample=SAMPLES),
        
        #==============================================================================================================

        # 02 - Alignment
        
        #==============================================================================================================

        # bowtie2_03 outputs
        expand("results/bowtie2/{sample}.bam", sample=SAMPLES),
        

        #===============================================================================================================        
        
        # 03 - Post-alignment and Mitochondrial Filtering 

        #===============================================================================================================
        
        # samtools_sort_04 outputs
        expand("results/samtools_sort/{sample}.sorted.bam", sample=SAMPLES),
 
        # calculate mito reads outputs
        expand("results/mito-ATAC/{sample}_mito_stats.txt", sample=SAMPLES),
       
        # remove mitochondrial reads counts
        expand("results/remove_mito_reads/{sample}_noMT.sorted.bam", sample=SAMPLES), 

        # samtools index outputs
        expand("results/samtools_index/{sample}_noMT.sorted.bam.bai", sample=SAMPLES),
        
        # samtools fixmate outputs
        expand("results/samtools_fixmate/{sample}_noMT.sorted.fixmate.bam", sample=SAMPLES), 

        #samtools markdup outputs
        expand("results/samtools_markdup/{sample}_noMT.sorted.dedup.bam", sample=SAMPLES), 

        # samtools_index_post_markdup outputs
        expand("results/samtools_index/post_markdup/{sample}_noMT.sorted.dedup.bam.bai", sample=SAMPLES),
     
        # samtools_view (filters dedup bam) outputs
        expand("results/samtools_view/{sample}.filtered.bam", sample=SAMPLES),
        
        #samtools_index_post_filter outputs
        expand("results/samtools_view/{sample}.filtered.bam.bai", sample=SAMPLES), 

        # tn5_shift outputs
        expand("results/tn5_shift/{sample}.filtered.shifted.bam", sample=SAMPLES),
        expand("results/tn5_shift/{sample}.filtered.shifted.bam.bai", sample=SAMPLES),
         
        # samtools_stats outputs
        expand("results/samtools_stats/{sample}_postFiltering.stats.txt", sample=SAMPLES),
        
        #fragment analysis outputs
        expand("results/fragment_size_analysis/{sample}_fragment_sizes.txt", sample=SAMPLES), 
        expand("results/fragment_size_analysis/{sample}_fragment.png", sample=SAMPLES), 
        expand("results/fragment_size_analysis/{sample}_fragment_stats.txt", sample=SAMPLES), 
        

        #=============================================================================================================

        # 04 - Picard Metrics
 
        #=============================================================================================================
        
        # picard CollectAlignmentSummaryMetrics outputs
        expand("results/picard/CollectAlignmentSummaryMetrics/{sample}.alignment_metrics.txt", sample=SAMPLES),
        
        # picard CollectInsertSizeMetrics outputs
        expand("results/picard/CollectInsertSizeMetrics/{sample}.insert_metrics.txt", sample=SAMPLES),
        expand("results/picard/CollectInsertSizeMetrics/{sample}.insert_histogram.pdf", sample=SAMPLES),
        
        #=============================================================================================================

        # 05 - TSS Enrichment  
 
        #=============================================================================================================
 
        #tss enrichment outputs
        expand("results/tss_enrichment/{sample}_tss_enrichment.txt", sample=SAMPLES), 
        expand("results/tss_enrichment/{sample}_tss_enrichment.pdf", sample=SAMPLES),

        
        #=============================================================================================================

        # 06 - Coverage and BigWig 

        #=============================================================================================================


        # bedtools_genomecov outputs
        expand("results/bedtools_genomecov/{sample}.bedGraph", sample=SAMPLES),
        
        # sorted_bedgraph outputs
        expand("results/sorted_bedgraph_file/{sample}.sorted.bedGraph", sample=SAMPLES),
        
        # bigwig_conversion outputs
        expand("results/bigwig/{sample}.bw", sample=SAMPLES),

        #normalized coverage outputs
        expand("results/normalized_coverage/{sample}_CPM.bw", sample=SAMPLES), 
        
        #=============================================================================================================

        # 07 - Correlation Analysis

        #=============================================================================================================

        #Correlation analysis outputs
        "results/correlation_analysis/matrix.npz", 
        "results/correlation_analysis/matrix.tab", 
        "results/correlation_analysis/correlation_heatmap.png", 
        "results/correlation_analysis/correlation_values.tab", 

           
        #=============================================================================================================

        # 08 - Peak Calling and Filtering 

        #=============================================================================================================

        # macs2_peak_calling outputs
        expand("results/macs2_peakcall/{sample}_peaks.narrowPeak", sample=SAMPLES),

        #blacklist region filter outputs
        expand("results/filtered_peaks/{sample}_filtered_peaks.bed", sample=SAMPLES), 
        
            
        #=============================================================================================================

        # 09 - HeatMap and FriP Calculation

        #=============================================================================================================

        #Heatmap visualization and correlation
        expand("results/heatmap/matrix/{sample}_matrix.gz",  sample=SAMPLES), 
        expand("results/heatmap/{sample}_regions.bed", sample=SAMPLES), 
        expand("results/heatmap/plot/{sample}_tss_heatmap.pdf", sample=SAMPLES), 
        
        
        #frip scores
        expand("results/frip_calculation/{sample}_frip.txt", sample=SAMPLES), 

           
        #==============================================================================================================

        # 10 - Peak Annotation and Motif Analysis

        #==============================================================================================================
        
        #Peak annotation
        expand("results/peak_annotation/{sample}_peak_annotation.txt", sample=SAMPLES), 

        #motif analysis
        expand("results/motif_analysis"), 
        
           
        #==============================================================================================================

        # 11  - QC Metrics

        #==============================================================================================================

        # phantompeakqualtools
        #expand("results/phantompeakqualtools/{sample}_qc.txt", sample=SAMPLES),
        #expand("results/phantompeakqualtools/{sample}_qc.pdf", sample=SAMPLES),
        
        # preseq
        expand("results/preseq/{sample}.ccurve.txt", sample=SAMPLES),
        
        # qualimap_bamqc
        expand("results/qualimap/{sample}_qualimap_report", sample=SAMPLES),

           
        #===============================================================================================================

        # 12 - MULTIQC
        #===============================================================================================================

        # multiqc
        "results/multiqc"
        ]
