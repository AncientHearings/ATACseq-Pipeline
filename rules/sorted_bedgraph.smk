rule sorted_bedgraph:
    input:
        bedgraph=lambda wildcards: f"{config['sorted_bedGraph']['input']['bedGraph']/{wildcards.sample}.bedGraph"
        
    output:
        sorted_bedgraph=temp(config['sorted_bedGraph']['output']['sorted_bedgraph'] + "/{sample}.sorted.bedGraph")
        
    resources:
        mem_mb=int(config['sorted_bedGraph']['params']['mem_mb'].replace("G", "000"))
    
    benchmark:
        "benchmarks/sorted_bedgraph/{sample}.txt"
            
    log:
        "logs/sorted_bedgraph/{sample}.err"
        
    threads:
        config['sorted_bedGraph']['params']['threads']['default']
        
    message:
        "Sorting bedgraph file ..."
        
    shell:
        """
        sort \
        -k1,1 -k2,2n \
        --parallel {threads} \
        -S {resources.mem_mb}M \
        {input.bedgraph} \
        > {output.sorted_bedgraph} \
        2> {log}
        """
        
