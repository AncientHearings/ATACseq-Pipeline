rule bowtie2_align: 
    input:
        R1_fastp=lambda wildcards: "config['bowtie2']['input']["R1"]" + "{wildcards.sample}_R1.fastq.gz", 
        R2_fastp=lambda wildcards: "config['bowtie2']['input']["R2"]" + "{wildcards.sample}_R2.fastq.gz"
        
    output:
        BAM="config['bowtie2']['output']" + "{sample}.bam"
        
    params:
        index = config['bowtie2']['params']['index']        
    
    benchmark: 
        "benchmarks/bowtie2/{sample}.txt"
        
    log: 
        stdout="log/bowtie2/stdout/{sample}.out", 
        stderr="log/bowtie2/stderr/{sample}.err"
        
    conda: 
        "envs/02_alignment/bowtie2.yaml"
        
    threads:
        config['bowtie2']['params']['threads']
        
     message:
        "Aligning reads..."
         
     shell:
        """
        bowtie2 -x {params.index} \
        -1 {input.R1_fastp} \
        -2 {input.R2_fastp} \
        --threads {threads} \        
        2> {log.stderr} | \
        samtools view -Sb - > {output.BAM} \
        2> {log.stderr} \
        > {log.stdout}
        """
         
         
         
   
