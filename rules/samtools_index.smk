rule samtools_index:
    input:
        sorted_bam=lambda wildcards: config['samtools_index']['input']['sorted'] + "/{wildcards.sample}.sorted.bam"
        
    output:
        indexed_bam=config['samtools_index']['output']['index'] + "/{wildcards.sample}.sorted.bam.bai"
        
    benchmark:
        "benchmarks/samtools_index/{sample}.txt"
        
    log:
        stdout="logs/samtools_index/stdout/{sample}.out"
        stderr="logs/samtools_index/stderr/{sample}.err"
        
    conda:
        "envs/02_alignment/post_alignment/samtools/samtools.yaml"
        
    threads:
        config['samtools_index']['params']['threads']
        
    message:
        "Producing BAM index file..."
        
    shell:
        """
        samtools \
        -@ {threads} \
        {input.sorted_bam}\
        > {log.stdout} \
        2> {log.stderr}
        """
         
