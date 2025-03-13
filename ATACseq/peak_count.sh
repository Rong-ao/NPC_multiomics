#!/bin/bash

peaks=$(ls peaks/*_merge/*_peaks.narrowPeak)
cat $peaks | cut -f 1-3 |sort -k1,1 -k2,2n -|bedtools merge -i stdin | awk -v OFS='\t' 'BEGIN{num=0}{num++; print $0, "peak"num}' - |Rscript keep_regular_chrom.R > merged_peaks.bed
bedtools intersect -a merged_peaks.bed -b  <(cat $peaks | awk -v OFS="\t" '{if ($10==-1){ print $1,int(($2+$3)/2),int(($2+$3)/2) }  else {print $1,$2+$10,$2+$10} }' ) -loj > merged_peaks.all_summits.txt

awk -v OFS='\t' '{ if (NF >3) {print $4,$1,$2,$3,$6 } 
  else { print $1":"$2"-"$3,$1,$2,$3,"." } }' merged_peaks.bed > merged_peaks.saf

echo "conversion to SAF done"
## count the reads
if [ ! -d counts ]; then mkdir counts; fi 
files=$(ls bam/*_?.nodup.bam)
#featureCounts -a peaks/atac_peaks.saf -o counts/atac.read.counts $files -F SAF -T 8
#featureCounts -a peaks/atac_peaks.saf -o counts/atac.frag.counts $files -F SAF -T 8 -p 
featureCounts -a merged_peaks.saf -o merged_peaks.counts $files -F SAF -T 8 
echo "feature counts done"
