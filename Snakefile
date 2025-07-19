SAMPLES: 
    - sample1
    - sample2

configfile: "../config.yaml"

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
        
