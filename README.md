# Maternal Sedentary Behavior and Breast Milk EV miRNA Across Pregnancy

## Overview

This repository contains the complete analysis pipeline for the research study: **"Maternal sedentary behavior and breast milk EV miRNA across pregnancy: evidence for early-mid gestation sensitivity"** by Beatriz B. Elliott.

This study investigates the relationship between maternal sedentary behavior and extracellular vesicle (EV) microRNA (miRNA) expression in breast milk across different pregnancy trimesters and postpartum. The analysis examines how prolonged sedentary bouts (≥60 minutes) and moderate-to-vigorous physical activity (MVPA) correlate with miRNA expression patterns.

## Repository Structure

```
.
├── Anaysis&Plots.Rmd          # Main R Markdown analysis file
├── ColData.csv                # Sample metadata with activity measurements
├── Merged_Count.csv           # miRNA count data matrix
├── Pathway Analysis Files/    # Pathway enrichment analysis results
│   ├── T2_sednum60_pathways.csv
│   └── background_gene_list_panther.txt
├── Raw Data Processing        # Bash script for processing raw sequencing data
└── README.md                  # This file
```

## Study Design

The study analyzes breast milk samples collected from pregnant participants across multiple time points:
- **Trimester 1 (T1)**: Early pregnancy
- **Trimester 2 (T2)**: Mid-pregnancy
- **Trimester 3 (T3)**: Late pregnancy
- **Postpartum (PP)**: After delivery

### Activity Measurements

Two main physical activity metrics were assessed:
1. **Sedentary Behavior (sed60num)**: Minutes per day spent in sedentary bouts ≥60 minutes
2. **MVPA (mvpa)**: Minutes per day spent in moderate-to-vigorous physical activity

Participants were categorized into quartiles:
- **Low**: ≤ Q1 (quartile 1)
- **Interquartile range**: Q1-Q3
- **High**: ≥ Q3 (quartile 3)

## Data Files

### ColData.csv
Sample metadata file containing:
- Participant IDs
- MVPA measurements and quartile classifications for each time point
- Sedentary behavior measurements (sed60num) and quartile classifications for each time point

### Merged_Count.csv
miRNA expression count matrix:
- Rows: miRNA features (e.g., hsa-let-7a-5p, hsa-miR-21-5p)
- Columns: Sample IDs
- Values: Raw read counts for each miRNA

### Pathway Analysis Files
- **T2_sednum60_pathways.csv**: GO Biological Process enrichment results for Trimester 2 sedentary behavior analysis
- **background_gene_list_panther.txt**: Background gene list for pathway analysis

## Analysis Pipeline

The analysis consists of several key components:

### 1. Raw Data Processing
The `Raw Data Processing` script performs the following steps:
- **Quality Control**: FastQC analysis of raw sequencing reads
- **Adapter Trimming**: Using bbduk.sh to remove adapter sequences
- **Length Filtering**: Cutadapt to select reads between 18-30 nucleotides
- **Alignment**: Bowtie alignment to human miRBase reference
- **Quantification**: Counting aligned reads per miRNA using samtools

### 2. Quality Control & Normalization (`Anaysis&Plots.Rmd`)
- Filtering low-count miRNAs (≥5 counts in ≥3 samples)
- DESeq2 normalization
- Library size distribution analysis
- Variance Stabilizing Transformation (VST)

### 3. Principal Component Analysis (PCA)
- VST transformation of count data
- Selection of top 500 variable genes
- PCA calculation and visualization
- Correlation analysis between PCs and activity metrics
- Bootstrap confidence interval estimation

### 4. Differential Expression Analysis
Using DESeq2, the analysis compares:
- **Sedentary Behavior**: Low vs. High sedentary time (Q1 vs. Q3)
- **Active Behavior**: Low vs. High MVPA (Q1 vs. Q3)

Analyses are performed separately for:
- Trimester 1
- Trimester 2
- Trimester 3
- Postpartum

Differentially expressed miRNAs are identified using an adjusted p-value threshold of 0.10.

### 5. Visualization
- **Volcano Plots**: Display log2 fold change vs. adjusted p-value
- **Bubble Plots**: Show enriched GO pathways with fold enrichment and significance
- **PCA Plots**: Visualize sample clustering by activity quartiles
- **Regression Plots**: Illustrate relationships between PCs and activity metrics

## Key Findings

The analysis generates several output files:
- `results_activity_T1_DEGs_NAomit.csv`: DEGs for sedentary behavior in T1
- `results_activity_T2_DEGs_NAomit.csv`: DEGs for sedentary behavior in T2
- `results_activity_T3_DEGs_NAomit.csv`: DEGs for sedentary behavior in T3
- `results_activity_PP_DEGs_NAomit.csv`: DEGs for sedentary behavior postpartum
- `results_mvpa_T1_DEGs_NAomit.csv`: DEGs for MVPA in T1
- `DIANA_background_miRNAs.txt`: List of miRNAs passing QC filters for pathway analysis

## Requirements

### R Dependencies
```r
# Bioconductor packages
BiocManager::install("DESeq2")
BiocManager::install("vsn")

# CRAN packages
install.packages(c(
  "pheatmap",
  "RColorBrewer",
  "ggplot2",
  "tidyverse",
  "dplyr",
  "matrixStats",
  "EnhancedVolcano",
  "readr"
))
```

### Command-line Tools (for raw data processing)
- FastQC
- bbduk.sh (from BBTools suite)
- Cutadapt
- Bowtie
- Samtools
- Human miRBase reference database

## Usage

### Running the Main Analysis

1. Ensure all required R packages are installed
2. Place data files (`ColData.csv` and `Merged_Count.csv`) in the working directory
3. Open and run `Anaysis&Plots.Rmd` in RStudio or with:

```r
rmarkdown::render("Anaysis&Plots.Rmd")
```

### Processing Raw Sequencing Data

1. Update the directory paths in the `Raw Data Processing` script:
   - `INPUTDIR`: Directory containing raw FASTQ files
   - `OUTPUTDIR`: Directory for processed output
   - `HUMANMIRBASEINDEX`: Path to miRBase reference index

2. Run the bash script:
```bash
bash "Raw Data Processing"
```

## Statistical Methods

- **Differential Expression**: DESeq2 with Wald test
- **Multiple Testing Correction**: Benjamini-Hochberg (FDR)
- **Significance Threshold**: Adjusted p-value < 0.10
- **Correlation Analysis**: Pearson and Spearman with bootstrap confidence intervals (B=5000)
- **Pathway Enrichment**: GO Biological Process with Bonferroni correction

## Results Interpretation

### Volcano Plots
- X-axis: log2 fold change (effect size)
- Y-axis: -log10(adjusted p-value) (significance)
- Points colored by significance and fold change thresholds

### Bubble Plots
- Size: Fold enrichment of pathway
- Color: -log10(Bonferroni p-value)
- Y-axis: Enriched GO terms

### PCA Plots
- Convex hulls show grouping by activity quartiles
- PC2 appears to be sensitive to sedentary behavior in Trimester 2

## Citation

If you use this code or data, please cite:

*Elliott, B.B. (2026). Maternal sedentary behavior and breast milk EV miRNA across pregnancy: evidence for early-mid gestation sensitivity.*

## Contact

For questions or additional information, please contact:
- **Author**: Beatriz B. Elliott
- **Repository**: [github.com/bpereira2022/Maternal-sedentary-behavior-and-breast-milk-EV-miRNA-across-pregnancy](https://github.com/bpereira2022/Maternal-sedentary-behavior-and-breast-milk-EV-miRNA-across-pregnancy)

## License

Please contact the author regarding data usage and licensing.

---

*Last updated: January 2026*