rule samtools_index_postmarkdup:
    input:
        sorted_bam=lambda wildcards: f"{config['samtools_index_post_markdup']['input']['markdup']}/{wildcards.sample}.dedup.bam"
        
    output:
        indexed_bam=config['samtools_index_post_markdup']['output']['index'] + "/{wildcards.sample}.dedup.bam.bai"
        
    benchmark:
        "../benchmarks/samtools_index/post_markdup/{sample}.txt"
        
    log:
        "../logs/samtools_index/post_markdup/{sample}.err"
        
    conda:
        "../envs/02_alignment/post_alignment/samtools/samtools.yaml"
        
    threads:
        config['samtools_index_post_markdup']['params']['threads']
        
    message:
        "Producing BAM index file..."
        
    shell:
        """
        samtools \
        -@ {threads} \
        {input.sorted_bam}\
        2> {log.stderr}
        """
         
