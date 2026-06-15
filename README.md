# Maternal Sedentary Behavior and Breast Milk EV miRNA Across Pregnancy: Evidence for Early-Mid Gestation Sensitivity

**Author:** Beatriz B. Elliott  
**Last updated:** 2026-01-21

---

## Project Overview

This project investigates the relationship between maternal physical activity patterns, specifically prolonged sedentary behavior (SED60+) and moderate-to-vigorous physical activity (MVPA), and small extracellular vesicle (EV) miRNA profiles in breast milk across pregnancy and the postpartum period. miRNA sequencing was performed on 11 breast milk samples and analyzed for differential expression and principal component associations with activity measures collected across four timepoints: Trimester 1, Trimester 2, Trimester 3, and Postnatal.

---

## Repository Contents

| File | Description |
|------|-------------|
| `Raw_Data_Processing.txt` | Bash pipeline for raw FASTQ processing through to mature miRNA counts |
| `Merged_Count.csv` | Raw miRNA count matrix (2,655 miRNAs × 11 samples) |
| `ColData.csv` | Sample metadata including MVPA and SED60+ values and quartile classifications across all four timepoints |
| `Analysis_Plots_NOU_Activity.Rmd` | R Markdown file containing all QC, PCA, correlation, DESeq2, and figure generation code |
| `DIANA_background_miRNAs.txt` | Background miRNA list for pathway enrichment (generated on run) |
| `T2_sednum60_pathways.csv` | GO pathway enrichment results for Trimester 2 SED60+ upregulated miRNA targets (input for Figure 8A) |

---

## Raw Data Processing Pipeline (`Raw_Data_Processing.txt`)

Run in bash. Requires: `fastqc`, `bbduk.sh`, `cutadapt`, `bowtie`, `samtools`, and a human miRBase Bowtie index.

**Steps per sample:**

1. **FastQC** — quality control on raw FASTQ files
2. **bbduk.sh** — adapter trimming (Illumina TruSeq small RNA adapter)
3. **Cutadapt** — read length filtering (18–30 bp)
4. **Bowtie** — alignment to human mature miRNA sequences from miRBase (`-n 0 -l 32 --norc --best --strata -m 1`; unique alignments only)
5. **samtools** — SAM → BAM conversion, sorting, and indexing
6. **samtools idxstats** — per-miRNA read count extraction

Output per sample: a count file in `$OUTPUTDIR/maturemiRNAcounts/`. These are merged across samples to produce `Merged_Count.csv`.

**Before running**, update the following paths at the top of the script:
```bash
INPUTDIR="/path/to/your/fastq/files"
OUTPUTDIR="/path/to/your/output/directory"
HUMANMIRBASEINDEX="/path/to/miRBase/bowtie/index"
```

---

## R Analysis (`Analysis_Plots_NOU_Activity.Rmd`)

#### `load_data`
- Reads `Merged_Count.csv` and `ColData.csv`
- Aligns samples between count matrix and metadata
- Builds a `DESeqDataSet` and filters low-count miRNAs (≥5 counts in ≥3 samples)
- Saves background miRNA list to `DIANA_background_miRNAs.txt`
- Applies variance stabilizing transformation (VST, blind) → `count_mat`
- Defines shared `jn_theme` for ggplot figures

#### Figure 5A — Qubit RNA Concentrations
Bar plot of RNA concentrations (ng/µL) per sample as measured by Qubit. Saved as `Figure_4A.tiff`.

#### Figure 5B — Raw vs. Normalized Library Sizes
Side-by-side barplot comparing raw and DESeq2 size-factor-normalized total counts per sample. Saved as `Figure_4B.tiff`.

#### Figure 5C — PCA Plot
PCA on VST-normalized counts across all miRNAs. Samples labeled by ID. Saved as `Figure_4C.tiff`.

#### Figure 5D — Top 20 Most Expressed miRNAs
Horizontal bar plot of the 20 miRNAs with the highest mean VST-normalized expression across all samples. Saved as `Figure_4D.tiff`.

#### `fig6` — Correlation Analysis + Regression Plots

**Correlation analysis** (printed to console, not saved as a figure):
- PCA is re-run using the **top 500 most variable miRNAs** (`rowVars`), with `scale. = TRUE`
- Pearson and Spearman correlations with bootstrapped 95% CIs (B = 5,000) are computed for PC1 and PC2 against each activity predictor across all four timepoints:
  - SED60+: `sed60num_trimester1/2/3/postpartum`
  - MVPA: `mvpa30num_trimester1/2/3/postpartum`
- To change the predictor being examined, update the `predictor` variable

**Regression plots:**
- 4-panel figure: PC2 vs. SED60+ minutes/day, one panel per timepoint, colored by quartile (Q1 = red, IQR = green, Q3 = blue), with inline Pearson and Spearman annotation
- Saved as `Figure_6_SED60_Maternal_Activity_ISEV.tiff`
- A matching 4-panel MVPA figure is saved as `Figure_6_MVPA30_Maternal_Activity_ISEV.tiff`

#### `deseq2` — Differential Expression
DESeq2 analysis comparing low (Q1) vs. high (Q3) groups, excluding the interquartile range, for:
- SED60+ at T1, T2, T3, Postpartum → `results_activity_T1/T2/T3/PP_DEGs_NAomit.csv`

DEGs are defined as padj < 0.10. Results objects are retained in the environment for volcano plots.

#### `fig7` — Volcano Plots
EnhancedVolcano plots for SED60+ DESeq2 results at T1, T2, and T3. Thresholds: padj < 0.054, |log2FC| > 1. Colorblind-friendly palette. Saved individually as `Figure_7A/B/C.tiff` for assembly in PowerPoint.

#### `fig8A` — Bubble Plot
GO Biological Process enrichment bubble plot for targets of upregulated miRNAs in the Trimester 2 SED60+ comparison. Requires `T2_sednum60_pathways.csv`. Bubble size = fold enrichment; color = −log10(Bonferroni). Saved as `Figure_8A.tiff`.

---

## Activity Variables

| Variable | Description |
|----------|-------------|
| `sed60num_trimester1/2/3/postpartum` | Average minutes/day in sedentary bouts ≥60 min at each timepoint |
| `sed60num_quartile_trimester1/2/3/postpartum` | Quartile classification for SED60+ ("Lower than or equal to X", "Interquartile range", "Greater than or equal to X") |
| `mvpa30num_trimester1/2/3/postpartum` | Average minutes/day in MVPA bouts ≥30 min at each timepoint |
| `mvpa_quartile_trimester1/2/3/postpartum` | Quartile classification for MVPA |

---
