rule samtools_view:
    input:
        dedup_bam=lambda wildcards: config['samtools_view']['input']['markdup_bam'] + "/{wildcards.sample}.dedup.bam"
        
    output:
        filtered_bam=config['samtools_view']['output']['filtered_bam'] + "/{wildcards.sample}.filtered.bam"
    
    params:
        minimum_mapq=config['samtools_view']['params']['MAPQ']
        filter_flags=config['samtools_view']['params']['flags']
                
    benchmark:
        "benchmarks/samtools_view/{sample}.txt"
        
    log:
        stdout="logs/samtools_view/stdout/{sample}.out"
        stderr="logs/samtools_view/stderr/{sample}.err"
        
    conda:
        "envs/02_alignment/post_alignment/samtools/samtools.yaml"
        
    threads:
        config['samtools_view']['params']['threads']['default']
        
    message:
        "Filtering  deduplicated bam.."
        
    shell:
        """
        samtools view \
        -@ {threads} \
        -b \
        -q {params.minimum_mapq} \
        -F {params.filter_flags} \
        {input.dedup_bam} \
        -o {output.filtered_bam} \
        > {log.stdout} \
        2> {log.stderr}
        """

   
