rule samtools_stats:
    input:
        filtered_bam=lambda wildcards: config['samtools_stats']['input']['filtered_bam'] + "/{wildcards.sample}.filtered.bam"
        
    output:
        stats=config['samtools_stats']['output']['stats'] + "/{wildcards.sample}.stats.txt"
                        
    benchmark:
        "benchmarks/samtools_stats/{sample}.txt"
        
    log:
        stdout="logs/samtools_stats/stdout/{sample}.out"
        stderr="logs/samtools_stats/stderr/{sample}.err"
        
    conda:
        "envs/02_alignment/post_alignment/samtools/samtools.yaml"
        
    threads:
        config['samtools_stats']['params']['threads']['default']
        
    message:
        "Generating Statistics.."
        
    shell:
        """
        samtools stats \
        -@ {threads} \
        {input.filtered_bam} \
        > {output.stats} \
        > {log.stdout} \
        2> {log.stderr} 
        """

   
