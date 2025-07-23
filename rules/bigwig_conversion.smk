rule bigwig_conversion:
    input:
        sorted_bedgraph=lambda wildcards:f"{config['bigwig']['input']['sorted_bedgraph']/{wildcards.sample}.sorted.bedGraph"
        
    output:
        bigwig=protected(config['bigwig']['output']['bigwig']+"/{sample}.bw")
        
    params:
        genome=config['bigwig']['params']['genome']
        
    benchmark:
        "../benchmarks/bigwig/{sample}.txt"
        
    log:
        "../logs/bigwig/{sample}.err"
         
    conda:
        "../envs/02_alignment/post_alignment/bedGraph_to_bigwig.yaml"
        
    threads:
        config['bigwig']['params']['threads']['default']
            
    message:
       "Converting to bigwig..."
       
    shell:
        """
        bedGraphToBigWig \
        -@ {threads} \
        {input.sorted_bedgraph} \
        {params.genome} \
        {output.bigwig} \
        2> {log} 
        """
        
