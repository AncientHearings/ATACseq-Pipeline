#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
SAMPLES: 
    - sample1
    - sample2

configfile: "../config.yaml"

include: "../rule/bowtie2.smk"
include: "../rule/samtools_sort.smk"
include: "../rule/samtools_index.smk"
include: "../rule/samtools_markdup.smk"
include: "../rule/samtools_view.smk"
include: "../rule/samtools_stats.smk"
include: "../rule/picard_CollectAlignmentSummaryMetrics.smk"
include: "../rule/picard_CollectInsertSizerMetrics.smk"
include: "../rule/tn5_shift.smk"
include: "../rule/bedtools_genomecov.smk"



rule all:
    input:
        #02_alignment
        #bowtie2
        expand("results/output/{sample}.bam", sample=SAMPLES),
        
        #post_alignment
        #samtools
        #samtools_sort
        expand("results/samtools_sort/{sample}.sorted.bam", sample=SAMPLES),
        
        #samtools_index
        expand("results/samtools_index/{sample}.sorted.bam.bai", sample=SAMPLES), 
        
        #samtools_markdup
        expand("results/samtools_markdup/{sample}.dedup.bam", sample=SAMPLES), 
        
        #samtools_view
        expand("results/samtools_view/{sample}.filtered.bam", sample=SAMPLES), 
        
        #samtools stats
        expand("results/samtools_stats/{sample}.stats.txt", sample=SAMPLES), 
        
        #picard CollectAlignmentSummaryMetrics
        expand("results/picard/CollectAlignmentSummaryMetrics/{sample}.alignment_metrics.txt", sample=SAMPLES), 
        
        #picard CollectInsertSizeMetrics
        expand("results/picard/CollectInsertSizeMetrics/{sample}.insert_metrics.txt", sample=SAMPLES),
        expand("results/picard/CollectInsertSizeMetrics/{sample}.alignment_metrics.pdf", sample=SAMPLES), 
        
        #tn5 shift
        expand("results/tn5_shift/{sample}.shifted.bam", sample=SAMPLES), 
        expand("results/tn5_shift/{sample}.shifted.bam.bai", sample=SAMPLES),
        
        #bedtools_genomecov
        expand("results/bedtools_genomecov/{sample}.bedgraph", sample=SAMPLES), 
        
         
