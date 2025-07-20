rule bedtools_genomecov:
    input:
        sorted_bam=lambda wildcards: f"{config['bedtool_genomecov']['input']['sorted_bam']}+"/{wildcards.sample}.sorted.bam"
        
    output:
        bedgraph=config['bedtools_genomecov']['output']['bedgraph']+"/{sample}.bedgraph"
        
    params:
        extra=config['bedtools_genomecov']['params']['extra'],
        genome=config['bedtools_genomecov']['params']['genome']
        
    benchmark:
        "../benchmarks/bedtools_genomecov/{sample}.txt"
        
    log:
        "../logs/bedtools_genomecov/{sample}.err"
        
    conda:
        "../envs/02_alignment/post_alignment/bedtools.yaml"
        
    threads:
        config['bedtools_genomecov']['params']['threads']['default']
        
    message:
        "Running bedtools genomecov for sample {sample}"
        
    shell:
        """
        bedtools genomecov \
         -ibam {input.sorted_bam} \
         -g {params.genome} \
         {params.extra} \
         > {output.bedgraph} \
         2> {log} 
         """
