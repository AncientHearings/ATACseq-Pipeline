rule samtools_sort:
    input:
        bam=lambda wildcards: config['samtools_sort']['input']['bam'] + "/{wildcards.sample}.bam"
        
    output:
        bam_sorted=config['samtools_sort']['output']['sorted_bam'] + "/{sample}.sorted.bam"
        
    benchmark:
        "benchmarks/samtools_sort/{sample}.txt"
        
    log:
        stdout="logs/samtools_sort/{sample}.out", 
        stderr="logs/samtools_sort/{sample}.err"
        
    conda:
        "envs/02_alignment/post_alignment/samtools/samtools.yaml"
        
    threads:
        config["samtools_sort"]["params"]["threads"]["default"]
        
    message:
        "Sorting alignments by genomic coordinates..."
        
    shell:
        """
        samtools sort \
        -@ {threads} \
        -o {output.bam_sorted} \
        {input.bam}
        > {log.stdout} \
        2> {log.stderr} 
        """
        
        
        
