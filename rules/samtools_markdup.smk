rule samtools_markdup:
    input:
        sorted_bam=lambda wildcards: config['samtools_index']['input']['sorted'] + "/{wildcards.sample}.sorted.bam"
        
    output:
        deduplicated_bam=config['samtools_markdup']['output']['markdup_bam'] + "/{wildcards.sample}.dedup.bam"
    
    params:
        remove_duplicates=config['samtools_markdup']['params']['remove_duplicates']    
    
    benchmark:
        "benchmarks/samtools_markdup/{sample}.txt"
        
    log:
        stdout="logs/samtools_markdup/stdout/{sample}.out"
        stderr="logs/samtools_markdup/stderr/{sample}.err"
        
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
        > {log.stdout} \
        2> {log.stderr}
        """
