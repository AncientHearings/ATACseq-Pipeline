rule tn5_shift:
    input:
        filtered_bam=lambda wildcards: f"{config['tn5_shift']['input']['filtered_bam']}/{wildcards.sample}.filtered.bam"
        
    output:
        shifted_bam=config['tn5_shift']['output']['shifted_bam'] + "{sample}.shifted.bam",
        shifted_bam_index=config['tn5_shift']['output']['shifted_bam_index'] + "{sample}.shifted.bam.bai"
        
    benchmark:
        "../benchmarks/tn5_shift/{sample}.txt"
        
    log:
        "../logs/tn5_shift/{sample}.err"
        
    conda:
        "../envs/06_visualization/deeptools.yaml"
    threads:
        config['tn5_shift']['params']['threads']['default']
        
    message:
        "Adjusting ATAC-seq read positions by +4-5 bp to reflect true Tn5 cut sites"
        
    shell:
        """
        alignmentSieve --ATACshift \
           -b {input.filtered_bam} \
           -o {output.shifted_bam} \
           -p {threads} \
           2> {log}  && \
        samtools index {output.shifted_bam}
        """  
