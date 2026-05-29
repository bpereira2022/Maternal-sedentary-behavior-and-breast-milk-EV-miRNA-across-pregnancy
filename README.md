# Maternal Sedentary Behavior and Breast Milk EV miRNA Across Pregnancy

## Overview

This repository contains analysis code and data for:
**"Maternal sedentary behavior and breast milk EV miRNA across pregnancy: evidence for early-mid gestation sensitivity"**.

The project evaluates how sedentary behavior (sed60num) and moderate-to-vigorous physical activity (MVPA) relate to breast milk extracellular vesicle (EV) miRNA expression across pregnancy and postpartum.

## Repository Structure

```
.
├── Analysis_Plots_Updated.Rmd    # Main analysis and figure-generation workflow
├── ColData.csv                   # Sample metadata and activity variables
├── Merged_Count.csv              # miRNA count matrix
├── Pathway Analysis Files/
│   ├── T2_sednum60_pathways.csv
│   └── background_gene_list_panther.txt
├── Raw Data Processing           # Bash pipeline for raw small RNA-seq processing
└── README.md
```

## Study Timepoints

- Trimester 1 (T1)
- Trimester 2 (T2)
- Trimester 3 (T3)
- Postpartum (PP)

## Main Analysis Workflow

`Analysis_Plots_Updated.Rmd` performs:

1. Data loading and alignment of sample IDs (`Merged_Count.csv`, `ColData.csv`)
2. Low-count filtering and VST normalization (DESeq2)
3. PCA and variance summaries
4. Differential expression testing by sedentary-behavior and MVPA quartiles
5. Figure generation (TIFF outputs)
6. Export of DEG tables and auxiliary outputs

### Key generated outputs

- `DIANA_background_miRNAs.txt`
- `PCA_percent_variance.csv`
- `PCA_sample_scores.csv`
- `PCA_feature_loadings_topGenes.csv`
- `results_activity_T1_DEGs_NAomit.csv`
- `results_activity_T2_DEGs_NAomit.csv`
- `results_activity_T3_DEGs_NAomit.csv`
- `results_activity_PP_DEGs_NAomit.csv`
- `results_mvpa_T1_DEGs_NAomit.csv`
- `results_mvpa_T2_DEGs_NAomit.csv`
- `results_mvpa_T3_DEGs_NAomit.csv`
- `results_mvpa_PP_DEGs_NAomit.csv`
- `Figure_5B.tiff`, `Figure_5C.tiff`, `Figure_5D.tiff`
- `Figure_6_Maternal_Activity_ISEV.tiff`
- `Figure_7A.tiff`, `Figure_7B.tiff`, `Figure_7C.tiff`
- `Figure_8A.tiff`

## Raw Data Processing

`Raw Data Processing` is a bash script for small RNA-seq preprocessing:

- FastQC quality control
- Adapter trimming (bbduk.sh)
- Length filtering (cutadapt; 18–30 nt)
- Alignment to miRBase (bowtie)
- BAM processing and indexing (samtools)
- miRNA count extraction (`samtools idxstats`)

> Update `INPUTDIR`, `OUTPUTDIR`, and `HUMANMIRBASEINDEX` in the script before running.

## Requirements

### R/Bioconductor and CRAN packages

- BiocManager
- DESeq2
- vsn
- pheatmap
- RColorBrewer
- ggplot2
- matrixStats
- tidyverse
- patchwork
- EnhancedVolcano
- readr
- dplyr
- ggrepel (used for PCA point labels)

### Command-line tools (raw processing)

- FastQC
- bbduk.sh (BBTools)
- cutadapt
- bowtie
- samtools

## Usage

### Run the R analysis

From repository root:

```r
rmarkdown::render("Analysis_Plots_Updated.Rmd")
```

### Run raw data processing

```bash
bash "Raw Data Processing"
```

## Citation

If you use this repository, please cite:

*Elliott, B.B. (2026). Maternal sedentary behavior and breast milk EV miRNA across pregnancy: evidence for early-mid gestation sensitivity.*

## Contact

- Author: Beatriz B. Elliott
- Repository: <https://github.com/bpereira2022/Maternal-sedentary-behavior-and-breast-milk-EV-miRNA-across-pregnancy>
