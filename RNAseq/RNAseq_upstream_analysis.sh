#! /bin/bash

module load fastqc/0.11.9
module load star/2.7.8a
module load fastp/0.23.2
module load subread/2.0.2

FASTQ_DIR="/data/RawData" # directory where your FASTQ files are located

FASTQC_DIR="/data/fastqc_results" # directory to store FastQC results

TRIMMED_DIR="/data/trimmed_fastq" # directory to store trimmed FASTQ files

STAR_INDEX="/RefGenome/GRCh38.p13/star_index" # directory containing the STAR genome index

ALIGN_DIR="/data/star_alignments" # directory to store STAR alignment outputs

FASTP_JSON="/data/fastp.json"

FASTP_HTML="/data/fastp.html"

for FASTQ_FILE in $FASTQ_DIR/*_R1_001.fastq.gz; do

    BASENAME=$(basename "$FASTQ_FILE" _R1_001.fastq.gz)   
    
    fastqc -o $FASTQC_DIR "$FASTQ_FILE"
    fastqc -o $FASTQC_DIR "${FASTQ_FILE/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    
    fastp -i "$FASTQ_FILE" \
                         -I "${FASTQ_FILE/_R1_001.fastq.gz/_R2_001.fastq.gz}" \
                        -o "$TRIMMED_DIR/${BASENAME}_R1_trimmed.fastq.gz" \
                        -O "$TRIMMED_DIR/${BASENAME}_R2_trimmed.fastq.gz" \
                        --json "$FASTP_JSON/${BASENAME}_fastp.json" --html "$FASTP_HTML/${BASENAME}_fastp.html"

    STAR --genomeDir $STAR_INDEX \
         --readFilesIn "$TRIMMED_DIR/${BASENAME}_R1_trimmed.fastq.gz" "$TRIMMED_DIR/${BASENAME}_R2_trimmed.fastq.gz" \
         --readFilesCommand zcat \
         --outFileNamePrefix "$ALIGN_DIR/${BASENAME}_" \
         --outSAMtype BAM SortedByCoordinate --runThreadN 8

done

COUNT_DIR="/data/featureCounts_results" # directory to store featureCounts results

ANNOTATION="/RefGenome/GRCh38.p13/GCF_000001405.39_GRCh38.p13_genomic.gtf" # path to the gene annotation file

BAM_LIST=$ALIGN_DIR/*.bam
featureCounts -a $ANNOTATION -t exon -g gene_id --extraAttributes gene_name \
                  -o "$COUNT_DIR/counts.txt" \
                  -T 6 -p -B -C $BAM_LIST
