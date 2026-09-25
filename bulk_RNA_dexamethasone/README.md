# Bulk RNA-seq Analysis of the Dexamethasone Response in Human Airway Smooth Muscle Cells

## Project overview

This repository contains an independent computational analysis of a publicly available bulk RNA-seq dataset generated from primary human airway smooth muscle (ASM) cells treated with dexamethasone.

The purpose of this project was to develop and document a complete bulk RNA-seq analysis workflow, starting from publicly available sequencing data and proceeding through quality control, transcript quantification, differential expression analysis, functional enrichment, visualization, and literature-based validation.

The analysis was performed independently using R, Bioconductor, and command-line bioinformatics tools.

## Data and experimental ownership

**I did not perform the original biological experiments and I do not own the experimental data, sequencing data, or published findings associated with this study.**

The biological samples, dexamethasone treatment, RNA extraction, library preparation, sequencing, and original experimental work were performed by the authors of the original study.

The sequencing data analyzed in this repository are publicly available through the relevant public repositories. The original study and its authors retain ownership and credit for the experimental work and associated scientific findings.

This repository contains **my independent computational analysis of those publicly available data**. The code, analysis workflow, visualizations, and interpretations presented here were generated as part of this portfolio project.

The original study should be cited when referring to the experimental design, biological samples, sequencing data, or previously reported findings.

## Dataset

- **GEO accession:** GSE52778
- **ENA/SRA project:** PRJNA229998
- **Organism:** *Homo sapiens*
- **Biological system:** Primary human airway smooth muscle cells
- **Treatment:** Dexamethasone
- **Control:** Untreated
- **Analysis subset:** 8 paired-end RNA-seq samples from 4 donors
- **Experimental structure:** Paired donor design
- **Paper** Himes BE, Jiang X, Wagner P, Hu R et al. RNA-Seq transcriptome profiling identifies CRISPLD2 as a glucocorticoid responsive gene that modulates cytokine function in airway smooth muscle cells. PLoS One 2014;9(6):e99625. PMID: 24926665

The analysis uses publicly available FASTQ data associated with the original study. The raw sequencing files are **not included in this GitHub repository** because of their size and because they are third-party publicly available data.

## Analytical workflow

```text
Publicly available FASTQ data
            ↓
        FastQC
            ↓
       MultiQC
            ↓
     Adapter trimming
        and QC
            ↓
        Salmon
   transcript quantification
            ↓
        tximport
            ↓
        DESeq2
   differential expression
            ↓
      DEG filtering
            ↓
      Visualization
            ↓
   DAVID functional
      enrichment
            ↓
 Literature-based validation
```
## License

The source code and original materials created for this repository are distributed under the MIT License. See LICENSE for details.