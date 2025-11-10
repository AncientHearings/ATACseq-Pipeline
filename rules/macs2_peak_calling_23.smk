rule macs2_peak_calling:
    input:
        shifted_bam=lambda wildcards: f"{config['macs_peakcall']['input']['shifted_bam']}/{wildcards.sample}.filtered.shifted.bam"
        
    output:
        peaks=f"{config['macs_peakcall']['output']['peaks']}/{{sample}}_peaks.narrowPeak"
          
    params:
        gsize=config['macs_peakcall']['params']['genome_size'],
        qval=config['macs_peakcall']['params']['qvalue'],
        nomodel=config['macs_peakcall']['params']['nomodel'], 
        format=config['macs_peakcall']['params']['format'], 
        dir=directory(config['macs_peakcall']['output']['peaks'])

    resources:
        mem_mb=config['macs_peakcall']['resources']['mem_mb'], 
        time=config['macs_peakcall']['resources']['time']
            
    benchmark:
        "benchmarks/macs2/{sample}.txt"
        
    log:
        "logs/macs2/{sample}.err"
        
    conda:
        "envs/04_peak_calling/macs2.yaml"

    threads:
        config['macs_peakcall']['threads']
        
    message:
        "[MACS2 PEAKCALLING] SAMPLE:  {wildcards.sample} | Markdup_Bam: {input.shifted_bam} | Peaks: {output.peaks} | Genome Size: {params.gsize} | QVal: {params.qval} | Nomodel: {params.nomodel} | Model: {params.format}]"

    shell: 
        """
        macs2 callpeak \
            -t {input.shifted_bam} \
            -f {params.format} \
            -g {params.gsize} \
            -n {wildcards.sample} \
            --outdir {params.dir} \
            {params.nomodel} \
            -q {params.qval} \
            2> {log}
                
         """
         

