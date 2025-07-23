#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
SAMPLES: 
    - sample1
    - sample2

configfile: "../config.yaml"

include: "../rules/bowtie2.smk"
include: "../rules/samtools_sort.smk"
include: "../rules/samtools_index.smk"
include: "../rules/samtools_markdup.smk"
include: "../rules/samtools_view.smk"
include: "../rules/samtools_stats.smk"
include: "../rules/picard_CollectAlignmentSummaryMetrics.smk"
include: "../rules/picard_CollectInsertSizerMetrics.smk"
include: "../rules/tn5_shift.smk"
include: "../rules/samtools_markdup.smk"
include: "../rules/samtools_index_post_markdup.smk"
include: "../rules/bedtools_genomecov.smk"
include: "../rules/sorted_bedgraph.smk"
include: "../rules/

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
        
        #samtools_markdup
        expand("results/samtools_markdup/{sample}.dedup.bam", sample=SAMPLES), 
        
        #samtools_index_post_markdup
        expand("results/samtools_index/post_markdup/{sample}dedup.bam.index", sample=SAMPLES), 
                
        #bedtools_genomecov
        expand("results/bedtools_genomecov/{sample}.bedgraph", sample=SAMPLES), 
        
        #sorted bedgraph
        expand("results/sorted_bedgraph_file/{sample}.sorted.bedGraph", sample=SAMPLES),
        
        #bigwig conversion
        expand("results/bigwig/{sample}.bw", sample=SAMPLES), 
        
        

