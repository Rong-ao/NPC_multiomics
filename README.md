# NPC_multiomics
This projects aims to profiling cell cycle-related gene expression and regulation in the early time of human neuron development. 

The resporitory contains bioinformatics analysis scripts for cell-cycle related neural progenitor cell differentiation omics data, including __RNA-seq__, __ATAC-seq__ and __single-cell RNA-seq__ data. All requirements of softwares or environment are recorded in scripts.

Citation is still unpublished. 

__Thanks [@silentFUSU](https://github.com/silentFUSU) for creating parts of scripts and contribution to the project.__
## RNA-seq analysis
* RNA-seq rawdata quality control with *FastQC*, trimming with *fastp*, mapping with *STAR* and counting with *featureCounts*: __RNAseq_upstream_analysis.sh__
* Differential expression gene analysis with *DESeq2* and GO enrichment with *clusterProfiler*: __RNAseq1_DEG.R__ and __RNAseq1_DEG.R__


## ATAC-seq analysis
* ATAC-seq rawdata mapping with *Bowtie2*, deduplicating with *picard*, generate bigWig files with *deepTools* and peak calling with *MACS2*: __ATAC_mapping.snakefile__
* Counting peaks from upstream analysis with _featureCounts_: __peak_count.sh__ (__keep_regular_chrom.R__ is used for modifying chromosome names in this step)
* Finding motifs with *HOMER* in different samples of time point: __homer_findmotif.sh__


## Single-cell RNA-seq analysis
Single-cell RNA-seq in this project was generated with scMulti-omics sequencing from 10X Gemonics. Only scRNA-seq data was applied for this project currently. 

Rawdata of Single-cell RNA-seq was processed with _CellRanger-ARC 1.0.1_ and downstream analysis was mainly done with _Seurat (v5.0)_ in _R (4.2.1)_.

Please submit issues or contact with kourongao@westlake.edu.cn if any questions or other code requests.
