rule picard_CollectInsertSizeMetrics:
    input:
        sorted_bam=lambda wildcards: config['picard_collect_insert_size_metrics']['input']['bam_sorted'] + "/{wildcards.sample}.sorted.bam
        
    output:
        insert_metrics=config['picard_collect_insert_size_metrics']['output']['metrics'] + "/{sample}.insert_metrics.txt",
        insert_histogram=config['picard_collect_insert_size_metrics']['output']['histogram'] +  "/{sample}.insert_histogram.pdf"
        
    params:
        m=config['picard_collect_insert_size_metrics']['params']['M'],
        validation_stringency=config['picard_collect_insert_size_metrics']['params']['validation_stringency']
        
    benchmark:
        "../benchmarks/picard/CollectInsertSizeMetrics/{sample}.txt"        

    log:
        stderr="../logs/picard/CollectInsertSizeMetrics/stderr/{sample}.err"
        
    conda:
        "../envs/02_alignment/post_alignment/picard/picard_metrics.yaml"

    threads:
        config['picard_collect_insert_size_metrics']['params']['threads']['default']
        
    message:
        "Calculating Insert Size Statistics..."
        
    shell:
        """
        picard CollectInsertSizeMetrics \
        I={input.sorted_bam} \
        O={output.insert_metrics} \
        H={output.insert_histogram} \
        M={params.m} \
        VALIDATION_STRINGENCY={params.validation_stringency} \
        2> {log.stderr} 
        """
