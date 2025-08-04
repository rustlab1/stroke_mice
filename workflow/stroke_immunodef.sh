#!/bin/bash

# RNA-seq pipeline for WT Stroke vs WT Control (SINGLE-END)
# Project: PRJNA860208
# Structure: stroke_mice/{raw,fastq,trimmed,aligned,counts,logs,qc}

# -------------------- Setup --------------------

# Conda setup (for Apple Silicon)
# curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-MacOSX-arm64.sh
# bash Miniconda3-latest-MacOSX-arm64.sh
# eval "$(/Users/ruslanrust/miniconda3/bin/conda shell.bash hook)"

# Create environment
conda create -n stroke_rnaseq -c bioconda -c conda-forge \
  sra-tools fastqc multiqc hisat2 samtools trimmomatic subread -y
conda activate stroke_rnaseq

# Folder setup
mkdir -p ~/stroke_mice/{raw,fastq,trimmed,aligned,counts,logs,qc}
cd ~/stroke_mice/raw

# -------------------- Define SRA Runs --------------------

# WT Stroke and Control
STROKE=(SRR20315306 SRR20315307 SRR20315314 SRR20315326)
CTRL=(SRR20315316 SRR20315322 SRR20315324)

# -------------------- Download FASTQ --------------------

# Download .sra files
for r in "${STROKE[@]}" "${CTRL[@]}"; do
  prefetch "$r"
done

# Convert to FASTQ and gzip
for r in "${STROKE[@]}" "${CTRL[@]}"; do
  fasterq-dump -e 16 -p -O . "$r"
  gzip -f "${r}.fastq"
done

# Move to fastq dir
# Rename only selected FASTQs
mv SRR20315306.fastq.gz WT_Stroke_1.fastq.gz
mv SRR20315307.fastq.gz WT_Stroke_2.fastq.gz
mv SRR20315314.fastq.gz WT_Stroke_3.fastq.gz
mv SRR20315326.fastq.gz WT_Stroke_4.fastq.gz

mv SRR20315316.fastq.gz WT_Ctrl_1.fastq.gz
mv SRR20315322.fastq.gz WT_Ctrl_2.fastq.gz
mv SRR20315324.fastq.gz WT_Ctrl_3.fastq.gz

# Move only renamed WT files to fastq/
mv WT_*.fastq.gz ../fastq/

# -------------------- QC --------------------

cd ~/stroke_mice/fastq
fastqc *.fastq.gz -o ../qc --threads 16


# -------------------- Alignment --------------------

mkdir -p ~/stroke_mice/hisat2_index
cd ~/stroke_mice/hisat2_index

# Download and unpack HISAT2 mm10 index
curl -O ftp://ftp.ccb.jhu.edu/pub/infphilo/hisat2/data/mm10.tar.gz

cd ~/stroke_mice/hisat2_index
tar -xzf mm10.tar.gz

ls mm10


cd ~/stroke_mice
SAMPLES=(WT_Stroke_1 WT_Stroke_2 WT_Stroke_3 WT_Stroke_4 WT_Ctrl_1 WT_Ctrl_2 WT_Ctrl_3)

# Align each sample
for sample in "${SAMPLES[@]}"
do
  echo "Aligning ${sample}..."
  hisat2 -p 4 \
    -x hisat2_index/mm10/genome \
    -U fastq/${sample}.fastq.gz \
    2> logs/${sample}_hisat2.log | \
    samtools sort -@ 4 -o aligned/${sample}.bam
  samtools index aligned/${sample}.bam
  echo "${sample} alignment done."
done



# -------------------- Post-alignment QC --------------------

echo -e "\n Alignment summary (samtools flagstat):"
for sample in "${SAMPLES[@]}"
do
  echo -e "\n--- ${sample} ---"
  samtools flagstat aligned/${sample}.bam | head
done

echo -e "\n Last 5 lines of HISAT2 logs:"
for sample in "${SAMPLES[@]}"
do
  echo -e "\n--- ${sample}_hisat2.log ---"
  tail -n 5 logs/${sample}_hisat2.log
done


# -------------------- Quantification --------------------

cd ~/stroke_mice

# Download GENCODE mouse annotation
curl -O https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M32/gencode.vM32.annotation.gtf.gz
gunzip -f gencode.vM32.annotation.gtf.gz

# Count reads per gene
featureCounts -T 16 -t exon -g gene_name \
  -a gencode.vM32.annotation.gtf \
  -o counts/raw_counts_gene_sym.txt aligned/*.bam \
  &> logs/featureCounts_gene_sym.log

# Build clean matrix (GeneSymbol + samples)
{ printf "GeneSymbol\t"; head -n 2 counts/raw_counts_gene_sym.txt | tail -n 1 | cut -f7-; } > counts/final_counts_symbols.tsv
tail -n +3 counts/raw_counts_gene_sym.txt | \
  awk -v OFS="\t" '{ out=$1; for(i=7;i<=NF;i++) out=out OFS $i; print out }' >> counts/final_counts_symbols.tsv
sed -i '' '1 s|aligned/||g; 1 s|\.bam||g' counts/final_counts_symbols.tsv

# Preview
head counts/final_counts_symbols.tsv | column -t

echo "✅ RNA-seq pipeline complete: WT Stroke vs WT Control."
