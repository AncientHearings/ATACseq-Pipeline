SAMPLES: 
    - sample1
    - sample2

configfile: "../config.yaml"

rule all:
    input:
        #02_alignment
        #bowtie2
        expand("results/output/{sample}.bam", samples=SAMPLES),
        
        #post_alignment
        #samtools
        #samtools_sort
        expand("results/samtools_sort/{sample}.sorted.bam", samples=SAMPLES),
