#! /usr/bin/env bash
## Snakefile
####################
Samples = ['H9_ESC','NPC_D1','NPC_D2','NPC_D3','NPC_D6']
BOWTIE2_INDEX = "/bowtie2_index/hg38"
MARKDUP="/Pipelines/atac-seq-pipeline-snakemake/dependencies/picard.jar MarkDuplicates"

rule all:
  input: 
    expand("bam/{sample}_1.nodup.bam",sample=Samples),
    expand("bam/{sample}_1.nodup.bam.bai",sample=Samples),
    expand("bigWig/{sample}_1.nodup.bw",sample=Samples),

    expand("bam/{sample}_2.nodup.bam",sample=Samples),
    expand("bam/{sample}_2.nodup.bam.bai",sample=Samples),
    expand("bigWig/{sample}_2.nodup.bw",sample=Samples),

    expand("bam/{sample}_merge.nodup.bam",sample=Samples),
    expand("bam/{sample}_merge.nodup.bam.bai",sample=Samples),
    expand("bigWig/{sample}_merge.nodup.bw",sample=Samples),

    expand("peaks/{sample}_merge/{sample}_peaks.narrowPeak",sample=Samples),
    expand("peaks/{sample}_merge/{sample}_peaks.xls",sample=Samples),
    expand("peaks/{sample}_merge/{sample}_treat_pileup.bdg",sample=Samples),
    expand("peaks/{sample}_merge/{sample}_control_lambda.bdg",sample=Samples)

rule bowtie2_align1:
  output: 
    bam=temp("bam/{sample}_1.sorted.bam"),
    raw_qc = "qc/{sample}_1.raw.flagstat.qc",
    log="logs/{sample}_1.bowtie2.log"
  input:
#    lambda wildcards: FASTQ_DICT[wildcards.sample]
    "{sample}_1/{sample}_1_1.fq.gz",
    "{sample}_1/{sample}_1_2.fq.gz"
  threads: 10 
  run:
    print(input)
    if len(input) == 1:
      middle = "-U " + str(input)
    elif len(input) == 2:
#     middle = "-1 " + str(input[0]).replace("fastq","fastq_trim",1).replace(".fastq","_val_1.fq") + " -2 " + str(input[1]).replace("fastq","fastq_trim",1).replace(".fastq","_val_2.fq") + " -X 2000"
      middle = "-1 " + str(input[0]) + " -2 " + str(input[1])  + " -X 2000"
    shell(
    "bowtie2 -x {BOWTIE2_INDEX} "
    "{middle} "
    "-p {threads} 2> logs/{wildcards.sample}_1.bowtie2.log|"
    "samtools view -bS |"
    "samtools sort -@ {threads} -m 4G > {output.bam};"
    "samtools flagstat {output.bam} > {output.raw_qc};"
    )

rule bowtie2_align2:
  output: 
    bam=temp("bam/{sample}_2.sorted.bam"),
    raw_qc = "qc/{sample}_2.raw.flagstat.qc",
    log="logs/{sample}_2.bowtie2.log"
  input:
#    lambda wildcards: FASTQ_DICT[wildcards.sample]
    "{sample}_2/{sample}_2_1.fq.gz",
    "{sample}_2/{sample}_2_2.fq.gz"
  threads: 10 
  run:
    print(input)
    if len(input) == 1:
      middle = "-U " + str(input)
    elif len(input) == 2:
#     middle = "-1 " + str(input[0]).replace("fastq","fastq_trim",1).replace(".fastq","_val_1.fq") + " -2 " + str(input[1]).replace("fastq","fastq_trim",1).replace(".fastq","_val_2.fq") + " -X 2000"
      middle = "-1 " + str(input[0]) + " -2 " + str(input[1])  + " -X 2000"
    shell(
    "bowtie2 -x {BOWTIE2_INDEX} "
    "{middle} "
    "-p {threads} 2> logs/{wildcards.sample}_2.bowtie2.log|"
    "samtools view -bS |"
    "samtools sort -@ {threads} -m 4G > {output.bam};"
    "samtools flagstat {output.bam} > {output.raw_qc};"
    )

rule bam_rmdup1:
  input:
    bam = "bam/{sample}_1.sorted.bam",
  output:
    bam = "bam/{sample}_1.nodup.bam",
    bai = "bam/{sample}_1.nodup.bam.bai",
    qc = "qc/{sample}_1.dup.qc"
  log:
    "logs/markdup/{sample}_1.markdup.log"
  threads: 3
  shell:
    "java -Xmx12G -XX:ParallelGCThreads=3 -jar {MARKDUP} TMP_DIR=tmp/{wildcards.sample}_1 INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.qc} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true 2> {log};"
    "samtools index {output.bam}"

rule bam_rmdup2:
  input:
    bam = "bam/{sample}_2.sorted.bam",
  output:
    bam = "bam/{sample}_2.nodup.bam",
    bai = "bam/{sample}_2.nodup.bam.bai",
    qc = "qc/{sample}_2.dup.qc"
  log:
    "logs/markdup/{sample}_2.markdup.log"
  threads: 3
  shell:
    "java -Xmx12G -XX:ParallelGCThreads=3 -jar {MARKDUP} TMP_DIR=tmp/{wildcards.sample}_2 INPUT={input.bam} OUTPUT={output.bam} METRICS_FILE={output.qc} VALIDATION_STRINGENCY=LENIENT ASSUME_SORTED=true REMOVE_DUPLICATES=true 2> {log};"
    "samtools index {output.bam}"

rule bam2bigwig1:
  input:
    bam = "bam/{sample}_1.nodup.bam"
  output: 
    bw = "bigWig/{sample}_1.nodup.bw"
  threads: 6
  shell:
    "bamCoverage -b {input.bam} -o {output.bw} --outFileFormat bigwig "
    "-bs 50 --numberOfProcessors {threads} --normalizeUsing RPKM"

rule bam2bigwig2:
  input:
    bam = "bam/{sample}_2.nodup.bam"
  output: 
    bw = "bigWig/{sample}_2.nodup.bw"
  threads: 6
  shell:
    "bamCoverage -b {input.bam} -o {output.bw} --outFileFormat bigwig "
    "-bs 50 --numberOfProcessors {threads} --normalizeUsing RPKM"

rule merge:
    input:
        bam1 = "bam/{sample}_1.nodup.bam",
        bam2 = "bam/{sample}_2.nodup.bam"
    output:
        bam = "bam/{sample}_merge.nodup.bam",
        bai = "bam/{sample}_merge.nodup.bam.bai"
    shell:
        "samtools merge {output.bam} {input.bam1} {input.bam2}; samtools index {output.bam}"

rule bam2bigwig_merge:
  input:
    bam = "bam/{sample}_merge.nodup.bam"
  output: 
    bw = "bigWig/{sample}_merge.nodup.bw"
  threads: 6
  shell:
    "bamCoverage -b {input.bam} -o {output.bw} --outFileFormat bigwig "
    "-bs 50 --numberOfProcessors {threads} --normalizeUsing RPKM"

rule macs2callpeaks:
    input:
        "bam/{sample}_merge.nodup.bam"
    output:
        "peaks/{sample}_merge/{sample}_peaks.narrowPeak",
        "peaks/{sample}_merge/{sample}_peaks.xls",

        "peaks/{sample}_merge/{sample}_treat_pileup.bdg",
        "peaks/{sample}_merge/{sample}_control_lambda.bdg"
    log:
        "logs/{sample}_merge_peaks.log"
    shell:
        "macs2 callpeak -t {input} -n {wildcards.sample} --outdir peaks/{wildcards.sample}_merge -g hs -B --SPMR --nomodel --shift 100 --extsize 200"
