#!/usr/bin/env bash

# Bulk RNA-Seq analysis - GEO:GSE52778
# Biological system: Primary human airway smooth muscle cells
# Experimental groups: Untreated and dexamethasone-treated
# Sample subset: 8 paired-end samples selected from the original dataset

set -euo pipefail

# Project directory 

cd ~/bulk_RNA_dexamethasone

# Sample IDs 

samples=(
    SRR1039508
    SRR1039509
    SRR1039512
    SRR1039513
    SRR1039516
    SRR1039517
    SRR1039520
    SRR1039521
)

# Step 1: Create project directories 

mkdir -p raw
mkdir -p fastqc/before_trimming
mkdir -p fastqc/after_trimming
mkdir -p trimmed
mkdir -p fastp_reports
mkdir -p multiqc/before_trimming
mkdir -p multiqc/after_trimming
mkdir -p reference/salmon_index
mkdir -p quant
mkdir -p figures
mkdir -p results
mkdir -p scripts

# Step 2: FASTQ quality control before trimming 

fastqc \
    raw/*.fastq.gz \
    -o fastqc/before_trimming

# Step 3: MultiQC before trimming 

multiqc \
    fastqc/before_trimming \
    -o multiqc/before_trimming

# Step 4: Adapter trimming and quality filtering with fastp 

for sample in "${samples[@]}"
do
    fastp \
        -w 8 \
        -i raw/${sample}_1.fastq.gz \
        -I raw/${sample}_2.fastq.gz \
        -o trimmed/${sample}_trimmed_1.fastq.gz \
        -O trimmed/${sample}_trimmed_2.fastq.gz \
        --detect_adapter_for_pe \
        -h fastp_reports/${sample}_fastp.html \
        -j fastp_reports/${sample}_fastp.json
done

# Step 5: FASTQ quality control after trimming 

fastqc \
    trimmed/*.fastq.gz \
    -o fastqc/after_trimming

# Step 6: MultiQC after trimming 

multiqc \
    fastqc/after_trimming \
    -o multiqc/after_trimming

# Step 7: Build Salmon transcriptome index 

if [[ ! -f reference/salmon_index/versionInfo.json ]]
then
    salmon index \
        -t reference/Homo_sapiens.GRCh38.cdna.all.fa \
        -i reference/salmon_index \
        -k 31
fi

# Step 8: Salmon transcript quantification 

for sample in "${samples[@]}"
do
    salmon quant \
        -i reference/salmon_index \
        -l A \
        -1 trimmed/${sample}_trimmed_1.fastq.gz \
        -2 trimmed/${sample}_trimmed_2.fastq.gz \
        -p 8 \
        -o quant/${sample}_quant
done

# Step 9: Continue downstream analysis in R 
# tximport, DESeq2, differential expression, pathway analysis, visualization, 
# and literature validation are performed using the R script in scripts/.