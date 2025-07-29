rule bowtie2_align: 
    input:
        R1_fastp=lambda wildcards: f"{config['bowtie2']['input']}/{wildcards.sample}_R1_trimmed.fastq.gz", 
        R2_fastp=lambda wildcards: f"{config['bowtie2']['input']}/{wildcards.sample}_R2_trimmed.fastq.gz"
        
    output:
        BAM=config['bowtie2']['output'] + "/{sample}.bam"
        
    params:
        index = config['bowtie2']['params']['index']        
    
    benchmark: 
        "../benchmarks/bowtie2/{sample}.txt"
        
    log: 
        "../log/bowtie2/{sample}.err"
        
    conda: 
        "../envs/02_alignment/bowtie2.yaml"
        
    threads: config['bowtie2']['params']['threads']
        
    message:
        "Aligning reads on {wildcards.sample}..."
         
    shell:
        r"""
        bowtie2 -x {params.index} \
        -1 {input.R1_fastp} \
        -2 {input.R2_fastp} \
        --threads {threads} 2> {log.stderr} | \
        samtools view -Sb - > {output.BAM} \
        2>> {log} \
        """
         
         
         
   
