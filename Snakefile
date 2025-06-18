#The config file is loaded containing the sample names.
configfile: "config.yaml" 
#This the final goal of my output. Run all the steps needed to get the desired #outputs. 
rule all:
    input:#These are the files I want at the end of the pipeline.
        expand("results/fastqc/{sample}_fastqc.html", sample=config["sample"])
 #Template for each FastQC output file.       
 #expand fills in {sample} with sample names  in config.yaml["sample"]. 
#A rule is defined that describes how to get the "outputs" from the "inputs" using  #the "shell command/script".
rule fastqc: 
  input:
    "raw_data/{sample}.fastq" 
    #{sample} is a wildcard. It gets automatically filled in by the sample name   
    #when the rule is executed.
 #These are the outputs expected.
  output: 
    html="results/fastqc/{sample}_fastqc.html",
    zip="results/fastqc/{sample}_fastqc.zip"    
#Saves the stdout + stderr of the fastqc command to the log file for debugging. 
  log:
    "logs/{sample}_fastqc.log" 
#Create a conda environment before running the command.    
  conda: 
    "envs/fastqc.yaml" 
#{log} 2>%1 redirects both stdout and stderr to log file    
  shell:
    """
    fastqc {input} --output results/fastqc > {log} 2>&1
    
    """
 
 SAMPLES = config["sample"]
 RAW_DATA = config["dir_raw"]
 TRIMMED_DATA = config["dir_trimmed_data"]
 THREAD = config["threads"]
   
rule fastp_trim:
    input:                                                   
       R1= lambda wildcards:f"{RAW_DATA}/{wildcards.sample}_R1.fastq.gz"
       R2 = lambda wildcards:f"{RAW_DATA}/{wildcards.sample}.fastq.gz"
    output:
       R1_TRIM = f"{TRIMMED_DATA}/{wildcards.sample}_R1.fastq.gz", 
       R2_TRIM = f"{TRIMMED_DATA}/{wildcards.sample}_R2.fastq.gz",
       html = f"{TRIMMED_DATA}/{wildcards.sample}_fastp.html", 
       json= f"{TRIMMED_DATA}/{wildcards.sample}_fastp.json"
    log:
       "logs/fastp/{sample}.log"
    conda:
       "envs/fastp.yaml"
    shell:
         """
         fastp \
         -i {input.R1} -I {input.R2} \
         -o {output.R1_TRIM} -O {output.R2_TRIM} \
         --detect_adapter_for_pe \
         --html {output.html} \
         --json {output.json} \
         > {log} 2>&1 
          """
#{log} >2&1: standard output and standard error to log file       

#............................................................  
#........................bowtie2.............................  
#............................................................  
INPUT_SAMPLES_BOWTIE2: config["inputs_bowtie2"]
INDEXES_BOWTIE2: config["index_bowtie2"]
OUTPUT_BOWTIE2: config["outputs_bowtie2"]
THREADS_BOWTIE2: config["threads"]


rule bowtie2_align:
     input: 
         cleaned_reads_1 =  lambda wc: f"{INPUT_SAMPLES_BOWTIE2}/{wc.sample}_Cleaned_R1.fastq.gz",
         cleaned_reads_2 =  lambda wc: f"{INPUT_SAMPLES_BOWTIE2}/{wc.sample}_Cleaned_R2.fastq.gz",
         bowtie2_indexes =  INDEXES_BOWTIE2
      output:
          bam = f"{OUTPUT_BOWTIE2}/{{sample}}.bam"   

      threads: THREADS_BOWTIE2

      conda: 
          "envs/bowtie2.yaml"
      log:
          "logs/bowtie2/{sample}.log"
      shell:
         """
         bowtie2 -x {input.bowtie2_indexes} \
         -1 {input.cleaned_reads_1} -2 {input.cleaned_reads_2}\
         -p {threads} |
         samtools view -Sb - > {output.bam} \
         2> {log} 
         """
#.............................................................        
#.............................................................        
#.............................................................

#.............................................................
#.......................SAMTOOLS..............................
#.........................VIEW................................
#CONVERT SAM TO BAM.
#SAM(Sequence Alignment/Map): Text-based format


INPUT_SAMTOOLS_VIEW: config["input_samtools_view"]
OUTPUT_SAMTOOLS_VIEW: config["output_samtools_view"]

rule samtool_view:
     input: 
        sam_file = lambda wc: f"{INPUT_SAMTOOLS_VIEW}/{wc.sample}.sam"
     output:
        bam_file = lambda wc: f"{OUTPUT_SAMTOOLS_VIEW}/{wc.sample}.bam"
     log:
        lambda wc: f"logs/{wc.sample}.log"
     conda:
        "envs/Samtools_view.yaml"
     shell:
        """
        samtools view -Sb \
        {input.sam_file} \
        > {output.bam_file} \
        2> {log}
        
        """
#.............................................................
#..........................samtools...........................       
#.............................sort............................

INPUT_SAMTOOLS_SORT = config["input_samtools_sort"]
OUTPUT_SAMTOOLS_SORT = config["output_samtools_sort"]
THREADS = 4

rule samtools_sort:
     input: 
       BAM = lambda wildcards: f"{INPUT_SAMTOOLS_SORT}/{wildcards.sample}.bam"
     output:
       bam_sorted = protected(lambda wildcards: f"{OUTPUT_SAMTOOLS_SORT}/{wildcards.sample}.sorted.bam")
     log:
       lambda wildcards: f"logs/samtools_sort/{wildcards.sample}.log"
     conda:
       "envs/samtools_sort.yaml"
     threads:      
        THREADS
     benchmark: 
        lambda wildcards: f"Benchmarks/samtools_sort/{wildcards.sample}.txt"
     message: 
        "Sorting, in progress, the BAM file for the sample {wildcards.sample} using the samtools sort"            
     shell:
        """
         samtools sort \
         -@{threads} \
         -o {output.bam_sorted} \
         {input.BAM} \
         2> {log} # Redirect the standard error to log file. 
        """ 
 #Space to avod concatenation.The comment is put outside to #avoid passing to the shell that might cause an error. 
  #"Might cause an error" means that the shell can get #confused as the \ means to continue on the next and the #comment interrupts that          
                 
      
       
   















#...............................................................
#...............................................................
#...............................................................







#.............................................................
#.....................samtools................................
#.......................index.................................

INPUT_SAMTOOLS_INDEX = config["input_samtools_sort"]
OUTPUT_SAMTOOLS_INDEX = config["output_samtools_index"]
THREADS =  4

rule samtools_index:
    input:
      bam_sorted = lambda wildcards: f"{INPUT_SAMTOOLS_INDEX}/{wildcards.sample}.sorted.bam"

    output:
      bam_indexed = lambda wildcards: f"{OUTPUT_SAMTOOLS_INDEX}/{wildcards.sample}.sorted.bam.bai"

    log:
      lambda wildcards: f"logs/samtools_index/{wildcards.sample}.log"

    conda:
      "envs/samtools_index.yaml"

    threads:
       THREADS

    shell:
      """
        samtools index \
	{input.bam_sorted}\
	2> {log}

      """
#.............................................................
#.............................................................
#.............................................................

#.............................................................
#......................samtools...............................
#.......................markdup...............................

#INPUT_1_SAMTOOLS_MARKDUP = config["input_1_samtools_markdup"]
#INPUT_2_SAMTOOLS_MARKDUP = config["input_2_samtools_markdup"]
OUTPUT_SAMTOOLS_MARKDUP = config["output_samtools_markdup"]
THREADS = config["threads"]

rule samtools_markdup:

     input:
       bam_sorted_for_samtools_markdup   = rules.samtools_sort.output.bam_sorted,
       bam_indexed_for_samtools_markdup  = rules.samtools_index.output.bam_indexed

     output:
       bam_markdup = lambda wildcards: f"{OUTPUT_SAMTOOLS_MARKDUP}/{wildcards.sample}.bam.markdup"

     benchmark:
       benchmark_markdup = lambda wildcards: f"Benchmarks/samtools_markdup/{wildcards.sample}.txt"

     log:
       samtools_markdup = lambda wildcards: f"logs/post_samtools_markdup_qc_tools_log_files/samtools_markdup/{wildcards.sample}_samtools_markdup.log"

     conda:
       "envs/SamtoolsMarkdup.yaml"

     threads:
       THREADS

     message:
       "Marking duplicates, in progress, by  the samtools markdup"

     shell:
      """
       samtools markdup -r -s \
       -@{threads} \
       {input.bam_sorted_for_samtools_markdup} \
       {output.bam_markdup} \
       2> {log.samtools_markdup}

      """
       
 #..............................................................
 #........................SAMTOOLS..............................
 #........................FLAGSTAT..............................
 
#INPUT_SAMTOOLS_FLAGSTAT = config["input_samtools_flagstat"]
#Snakefile is written in python syntax. Assignment operator is used to assign value to a variable. '=' is used instead of ':' as in #a yaml file.

OUTPUT_SAMTOOLS_FLAGSTAT = config["output_samtools_flagstat"]
THREADS = config["threads'] 
#config[] is to access the value from the variable threads in config.yaml file just like accessing value using key in python
#dictionary.

rule samtools_flagstat:
     input: 
       BAM_MARKDUP = rule.samtools_markdup.output.bam_markdup
     output:
       SAMTOOLS_FLAGSTAT = lambda wildcards: f"{OUTPUT_SAMTOOLS_F.FLAGSTAT}/{wildcards.sample}.flagstat.txt"
     benchmark:
       samtools_flagstat = lambda wildcards: f"Benchmarks/post_samtools_markdup_qc_tools_benchmark_files/samtools_flagstats/{wildcards.sample}.txt"
     log:
       samtools_flagstat = lambda wildcards: f"logs/post_samtools_markdup_qc_tools_log_files/samtools_flagstat/{wildcards.sample}_post_samtools_markdup_qc.log"
     conda:
       "envs/post_samtools_markdup_qc_env_files/samtools_flagstat.yaml"
     threads:
        THREADS
     message:
       "Generating general mapping statistics"
     shell:
       """
       samtools flagstat \
       -@{threads} \
       {input.BAM_MARKDUP} \
       > \
       {output.SAMTOOLS_FLAGSTAT} \
       2> {log}
       """
#................................................................
#................................................................
#...............................................................



#...............................................................
#.........................SAMTOOLS..............................
#...........................STATS...............................

OUTPUT_SAMTOOLS_STATS = config["output_samtools_stats"]
THREADS: threads

rule samtools_stats:
     input:
       INPUT_SAMTOOLS_STATS = rule.samtools_markdup.output.bam_markdup
     output:
       SAMTOOLS_STATS  = lambda wildcards: f"{OUTPUT_SAMTOOLS_STATS/{wildcards.sample}.stats.txt" 
     benchmark:
       samtools_stats = lambda wildcards: f"Benchmarks/post_samtools_markdup_qc_benchmark_files/samtools_stats/{wildcards.sample}.stats.txt"
     log:
       samtools_stats = lambda wildcards: f"logs/post_samtools_markdup_qc_tools_log_files/samtools_stats/{sample}_samtools.sort_post_samtools_markdup.log"
     conda:
       "envs/post_samtools_markdup_qc_tools_env_files/samtools_stats.yaml"
     message:
       "Providing detailed read statistics"
     shell:
       """
       samtools stats \
       -@{threads} \
       {input.INPUT_SAMTOOLS_STATS} \
       > \
       {output.SAMTOOLS_STATS} \
       2> {log}
       """
#...............................................................
#...............................................................
#...............................................................


     
#...............................................................
#........................Picard.................................
#............CollectAlignmentSummaryMterics.....................

OUTPUT_PICARD = config["output_picard"]
THREADS: config["threads"]

rule picard:
     input:
       INPUT_PICARD = rule.samtool_markdup.output.bam_markdup
     output:
       OUTPUT_PICARD = lambda wildcards: f"{OUTPUT_PICARD}/{wildcards.sample}.alignment_metrics.txt"
     benchmark:
       picard = lambda wildcards: f"Benchmarks/post_samtools_markdup_qc_tools_benchmark_files/picard/{wildcards.sample}.alignment_metrics.txt"
     log:
       picard = lambda wildcard: f"logs/post_samtools_markdup_qc_tools_log_files/picard/{wildcards.sample}_alignment_metrics_post_samtools_markdup.log"
     conda:
       "envs/post_samtools_markdup_qc_tools_env_files/Picard_CollectAlignmentSummaryMetrics.yaml"
     message:
       "Proviiding Detailed Alignment QC, paltform-specific stats"
     shell:
       """
       picard CollectAlignmentSummaryMetrics \
       R={params.ref} \ *
       I={input.INPUT_PICARD} \
       O={output.OUTPUT_PICARD} \
       VALIDATION_STRINGENCY=LENIENT \
       2> {log}
       """
#................................................................
#................................................................
#...............................................................



#................................................................
#.....................qualimap...................................
#......................bamqc.....................................

OUTPUT_QUALIMAP = config["output_qualimap"]
THREADS = config["threads"]

rule qualimap_bamqc:
     input:
       INPUT_QUALIMAP_BAMQC=rule.samtools_markdup.output.bam_markdup
     output:
       OUTPUT_QUALIMAP_BAMQC = lambda wildcards: f"{OUTPUT_QUALIMAP}/{wildcards.sample}_qualimapreport.html"
     benchmark:
       qualimap_bamqc = lambda wildcards: f"Benchmarks/post_samtools_markdup_qc_tools_benchmark_files/
       qualimap_bamqc/{wildcards.sample}qualimapreport.txt"
     log:
       qualimap_bamqc = lambda wildcards: f"logs/post_samtools_markdup_qc_tools_log_files/
       qualimap_bamqc/{wildcards.sample}_qualimapreport_post_samtools_markdup.log"
     conda:
       "envs/post_samtools_markdup_qc_tools_env_files/qualimap_bamqc.yaml:"
     message:
       "Generating visual summary, coverage distribution, GC basis"
     shell:
       """
       qualimap bamqc \
       -@{threads} \
       -bam \
       {input.INPUT_QUALIMAP_BAMQC} \
       -outdir \
       {output.OUTPUT_QUALIMAP_BAMQC} \
       -outformat HTML \
       2> {log}
       """
#...............................................................
#...............................................................
#...............................................................



#...............................................................
#.......................picard.................................
#.....................CollectInsertSizeFragements...............................

OUTPUT_INSERTSIZE = config["output"] ["picard_CollectInsertSizeFragements"]

OUTPUT_INSERTSIZE_VISUALIZATION = config ["output"] ["picard_CollectInsertSizeFragements"]

THREADS = config["threads"]

rule picard_CollectInsertSizeFragmentsMetricsAnd histogram:
     input:
       input_picard_CollectInsertSizeFragementMetrics = rule.samtools_markdup.output.bam_markdup
     
     output:  
       metrics = lambda wildcards: f"{OUTPUT_INSERTSIZE}/{wildcards.sample}.txt"
       histogram = lambda wildcards: f"{OUTPUT_INSERTSIZE_VISUALIZATION}/{wildcards.sample}.pdf"
       
     benchmark: 
       lambda wildcards: f"Benchmarks/picard_CollectInsertSizeFragements/{wildcards.sample}.txt"
       
     log: 
       stdout = lambda wildcards: f"logs/picard_CollectInsertSizeFragements/{wildcards.sample}.out"
       stderr = lambda wildcards: f"logs/picard_CollectInsertSizeFragements/{wildcards.sample}.err"
       
     conda:
       "envs/picard_CollectInsertSizeFragementsMetrics.yaml"
      
     threads:
       THREADS
      
     message: 
       " Assessing nucleosome periodicity"
       
     shell:
       """
       picard CollectInsertSizeMetrics \
       I=
       {input.input_picard_CollectInsertSizeFragementMetrics} \
       O={output.metrics} \
       H={output.histogram} \
       M=0.5 \
       > {log.stdout} 2> {log.stderr}
       """
 #By default M=0.5       
#It excludes the tails  in the histogram where the insert size is extremely rare.
#Only the central 50-% of the data will be used for calculating mean,etc.      
#To avoid skewed statistics from contaminant fragements.
#Need to read paper regarding the value. I may make it configurable by adding setting in the configuration file for future-proofing the pipeline.

#lambda patterns is used for dynamic output.
#config[] is for scalability.
#................................................................
#...............................................................
#...............................................................




#..............................................................
#........................preseq................................
#..............................................................
#To estimateif sequencing deeper will yeild more unique reads.
#...............................................................
#..............................................................
#..............................................................

OUTPUT_PRESEQ = config["output"] ["preseq"]
THREADS  = config["threads"]
#config[] is for scalability

rule  preseq_lc_exstrap:
      input: 
        in_preseq = rule.samtools_markdup.output.bam_markdup
      
      output: 
        out_preseq = lambda wildcards: f"{OUTPUT_PRESEQ}/{wildcards.sample}.txt"
        
      benchmark: 
        lambda wildcards: f"Benchmarks/preseq/{wildcards.sample}.txt"
         
      log:
        stdout = lambda wildcards: f"logs/preseq/{wildcards.sample}.out"
        stderr = lambda wildcards: f"logs/preseq/{wildcards.sample}.err"
        
     conda: 
        "envs/preseq.yaml"
        
     message: 
        "Estimating library complexity and unique reads with preseq"
        
     shell:
       """
        preseq lc_exstrap -B \
        -o {output.out_preseq} \
        {input.in_preseq} \
        >{log.stdout} 2> {log.stderr}
       """
#.............................................................
#.............................................................
#.............................................................   
 
 #Message is for makinglogs readable
 #Clean log handling - a hallmark of prodiction-grade pipeline   
 #resource and params parameter are kept on hold for now.   
#.............................................................
#.............................................................
#.............................................................
































#The rule syntax is used to create a rule with this name.
#This is used to manage dependencies internally*. 
   #Used for visualization in DAG. 
   #To run a specific rule.
   #Supports advanced rule management
#The rule parameters define how the rule behaves.


#..............................................................
#.........................deeptools............................
#.......................bamCoverage............................
#Instead of 



OUTPUT_BAMCOVERAGE = config["output_bamCoverage"]
THREADS = config["threads"]


rule bamCoverage: 
     input: 
       INPUT_BAMCOVERAGE_SAMTOOLS_MARKDUP_OUTPUT = rules.samtools_rule.output.bam_markdup, 
       INPUT_BAMCOVERAGE_SAMTOOLS_INDEXED = rules.samtools_indexed.output.bam_indexed
     output: 
       BAMCOVERAGE_OUTPUT = lambda wildcards: f"{OUTPUT_BAMCOVERAGE}/{wildcards.sample}_coverage.bw" 
     benchmark:
       lambda wildcards: f"Benchmark/bamCoverage/{wildcards.sample}.txt"    
     log:
       bamCoverage_stderr = lambda wildcards: f"logs/bamCoverage/{wildcards.sample}.err", 
       bamCoverage_stdout = lambda wildcards: f"logs/bamCoverage/{wildcards.sample}.out"
     conda: 
       "envs/bamCoverage.yaml"
     threads: 
       THREADS
     params: 
       binSize = config["BINSIZE"],
       normalize = config["NORMALIZE"],
     message: 
       "Generating BigWig for genome browsers"
     shell:
       """
       bamCoverage \
       -b {input.INPUT_BAMCOVERAGE_SAMTOOLS_MARKDUP_OUTPUT} \
       -o {output.BAMCOVERAGE_OUTPUT} \
       --binSize {params.binSize} \
       --normalizeUsing {params.normalize} \
       --numberOfProcessors {threads} \
       > {log.bamCoverage_stdout} 2> {log.bamCoverage_stderr}
       """
       
           
#Only one input file is mentioned explicitly but the other #must be present alongside to work efficiently.       
#params is a special snakemake keyword used to pass arbitrary parameters to the shell command within a rule.
#numberOfProcessors or numberOfThreads. confirm using bamCoverage --help

#..............................................................
#..............................................................
#..............................................................


#The rule synatx is used to define a rule. A rule starts with #with the rule syntax, followed by the  rule's name and #colon. The colon  indicates the start of the rule body. Inside the rule body, parameters like input, output, threads, params, log, and shell are provided. These parameters control the behavior and execution of the rule.


#..............................................................
#.........................MACS2................................
#.......................PeakCalling ............................

OUTPUT_MACS2 =  config["output_macs2_peak_calling"]
THREADS = config["threads"]

rule macs_PeakCalling:
     input:
       input_samtools_markdup =  temp(rule.samtools_markdup.output.bam_markdup)
        
     output:
       output_narrowPeak = protected(lambda wildcards: f"{OUTPUT_MACS2}/{wildcards.sample}_peaks.narrowPeak")
     
     params: 
     
       
     benchmark:
       lambda wildcards: f"Benchmarks/MACS2/{wildcards.sample}.txt"
       
     log:
       stdout = lambda wildacards: f"logs/MACS2/{wildcards.sample}.out"
       stedrr = lambda wildcards: f"logs/MACS2/{wildcards.sample}.err"
       
     conda:
       "envs/macs2_peakcall.yaml"
       
     threads: 
        THREADS
     
     resources:
     
         
       
     message:
       "Peak calling, im progress for {wildcards.sample}"
       
     shell:
       """
       macs2 callpeak \
       -t {input.input_samtools_markdup} \
       -f BAMPE \
       -g {params.genome_size} \
       -n {wildcards.sample} \
       -q {params.qval} \
       --outdir {params.outdir} \
       --keep-dup all \
       --call-summits \
       > {log.stdout} 2> {log.stderr}
       """
       
#Need to understand params and resources rule parameters
#...............................................................
#...............................................................
#...............................................................
#...............................................................        

#................................................................
#........................multiqc.................................
#................................................................
#expand() is a snakemake utility function. 

OUTPUT_MULTIQC_1 = config["output_multiqc"]["html"]
OUTPUT_MULTIQC_2 = config["output_multiqc"]["zip"] 
THREADS = config["threads"]

rule multiqc:
    input:
       expand("results/
    output:
       html = protect(lambda wildcards: f"{OUTPUT_MULTIQC_1}/{wildcards.sample}_report.html")
       zip = protect(lambda wildcards: f"{OUTPUT_MULTIQC_2}/{wildcards.sample}+_report.zip")
       
    benchmark:
       multiqc = lambda wildcards: f"Benchmarks/multiqc/multiqc_benchmark.txt"

    log:
      stdout = lambda wildacards: f"{logs/multiqc/{wildcards.sample}.out"}
      stderr = lambda wildcards: f"{logs/multiqc/{wildcards.sample}.err"
      
    params:
    
    resources:  
            
    conda:
      "envs/multiqc.yaml"
      
    threads:
       Threads
       
    message:
      "Generating a single readable HTML reports for all"

    shell:
      """
      multiqc \
      -@{threads} \
      
      """
      
#................................................................
#.....................Some confusion. Sort out during revision........................................................
#.................................................................
      
       
     
     
   
      
      
