# Mouse stroke RNA seq preprocessing and analysis

This repository provides an end to end workflow for RNA seq analysis of wild type mouse brain comparing stroke and control.  
It includes two parts

* Preprocessing in shell  
  downloads SRA runs, builds sample FASTQs, runs QC, aligns to mm10, and produces gene level counts

* Analysis in R Markdown  
  performs differential expression with DESeq2, PCA and volcano plots, GSEA, and generates the concise answer sheet

---

## 1. Data sources

* Public SRA BioProject PRJNA860208 wild type mouse brain single end RNA seq  
* Groups and runs used  
  * Stroke SRR20315306 SRR20315307 SRR20315314 SRR20315326  
  * Control SRR20315316 SRR20315322 SRR20315324
* Reference  
  * Genome mm10 HISAT2 index  
  * Annotation GENCODE vM32 GTF


---

## 2. How to download

Install the tools with conda and fetch FASTQs from SRA (see above)


## 3. Pre processing and subsampling
No subsampling is applied by default. Raw FASTQs are used as is. However only the WT group of mice after stroke and control was used
Quality control is performed with FastQC and can be summarized with MultiQC.


## 4. How the workflow works
The shell script executes the following steps. Save your script and run it from a terminal. It creates and uses ~/0_stroke_immuno_def.

Step 1 Quality control

Purpose check base quality and GC content
Tools FastQC and optional MultiQC
Inputs FASTQs in fastq
Outputs per sample HTML and zip in qc

Step 2 Reference setup

Purpose obtain the mm10 HISAT2 index
Tools curl tar
Outputs index in hisat2_index/mm10

Step 3 Alignment

Purpose align single end reads to mm10 and sort BAM
Tools HISAT2 samtools
Inputs FASTQs in fastq
Outputs BAM and BAI in aligned plus logs in logs

Step 4 Post alignment QC

Purpose review mapping statistics
Tools samtools

Step 5 Quantification

Purpose count reads per gene using gene symbols
Tools featureCounts
Inputs BAM files and GENCODE vM32 GTF
Outputs counts/final_counts_symbols.tsv

Step 6 Analysis

Purpose differential expression and pathway analysis and generation of the answer sheet
Tools R DESeq2 enrichplot clusterProfiler ggplot2
Inputs counts/final_counts_symbols.tsv and sample metadata
Outputs figures tables and a concise results file

In R

run stroke_immunodef.Rmd to load data and perform the main analysis
run Path_to_answers.Rmd to render the short answers once objects are in memory

Quick answers for this dataset
The Path_to_answers.Rmd document prints the following

Notes
You can skip raw FASTQ processing and start from counts/final_counts_symbols.tsv
Works on macOS and Linux with the conda environment above


