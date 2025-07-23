rule samtools_markdup:
    input:
        shifted_bam=lambda wildcards: config['tn5_shift']['output']['shifted_bam'] + "/{wildcards.sample}.shifted.bam"
        
    output:
        deduplicated_bam=config['samtools_markdup']['output']['markdup_bam'] + "/{wildcards.sample}.dedup.bam"
    
    params:
        remove_duplicates=config['samtools_markdup']['params']['remove_duplicates']    
    
    benchmark:
        "benchmarks/samtools_markdup/{sample}.txt"
        
    log:
        "logs/samtools_markdup/stderr/{sample}.err"
        
    conda:
        "envs/02_alignment/post_alignment/samtools/samtools.yaml"
        
    threads:
        config['samtools_markdup']['params']['threads']
        
    message:
        "Producing dedpulicated bam.."
        
    shell:
        """
        samtools markdup \
        -{params.remove_duplicates} \
        -@ {threads} \
        {input.sorted_bam} \
        {output.deduplicated_bam} \
        2> {log.stderr}
        """
