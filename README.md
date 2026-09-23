# Maternal Sedentary Behavior and Breast Milk EV miRNA Across Pregnancy: A Pilot Study

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
| `ColData.csv` | Sample metadata including MVPA and SED60+ values |
| `Analysis_Plots_NOU_Activity.Rmd` | R Markdown file containing all QC, PCA, correlation, DESeq2, and figure generation code |
| `QC` | Sample and library quality-control figures
| `PCA_and_Correlation` | PCA, PC-vs-activity correlations, loading diagnostics
| `Differential_Expression` | DESeq2 results (quartile model), volcano plots, per-sample counts
| `Sensitivity_Analyses` | Robustness checks for the differential expression and PCA results
| `Pathway Analysis Files` | Background miRNA list and GO pathway enrichment inputs/outputs

---

## Step 1: Raw data processing (FASTQ to miRNA counts)

**Script:** `Raw_Data_Processing.sh` (Bash)

For each sample, the script runs the following steps:

| Step | Tool | What it does |
|---|---|---|
| 1 | FastQC | Quality-control report for each FASTQ (`QCReports/`) |
| 2 | `bbduk.sh` (BBTools) | Trims the 3' Illumina adapter (`ktrim=r k=21 mink=11 hdist=2`) |
| 3 | Cutadapt | Keeps reads 18-30 nt long (typical miRNA size range) |
| 4 | Bowtie 1 | Aligns reads to a human-only mature miRNA index built from miRBase (`-n 0 -l 32 --norc --best --strata -m 1`: zero mismatches, forward strand only, uniquely-best alignments only) |
| 5 | SAMtools | Converts SAM to sorted, indexed BAM; removes the intermediate SAM |
| 6 | `samtools idxstats` | Extracts per-miRNA read counts into `<sample>-maturemiRNA-counts.txt` |

**Outputs (per sample):** trimmed FASTQ, length-filtered FASTQ, aligned/unaligned reads, sorted BAM + index, Bowtie and Cutadapt logs, and a two-column count table (`miRNA`, `<sample>-miRNAcount`).

**Before running:**

- Edit `INPUTDIR`, `OUTPUTDIR`, and `HUMANMIRBASEINDEX` at the top of the script to match your machine. The Bowtie index (`only_human_mirna_DNA`) must be built from miRBase human mature miRNA sequences in advance.
- Update the sample list in the `for SAMPLE in ...` line to include every sample.
- The per-sample count tables are merged into `Merged_Count.csv` (one column per participant, renamed to `NOU-XXXX` IDs) before the R analysis.

**Requires:** FastQC, BBTools (`bbduk.sh`), Cutadapt, Bowtie 1, SAMtools.
```

---

## Step 2: Statistical analysis (`Analysis_Plots_NOU_Activity.Rmd`)

Run top to bottom in RStudio (or knit to HTML). A fixed seed (`set.seed(123)`) is used for bootstrap resampling.

### 2.1 Setup and preprocessing

- Loads `Merged_Count.csv` and `ColData.csv` and aligns samples between them.
- **Quartile groups.** For each of the 8 exposure variables (SED60+ and MVPA × 4 timepoints), participants are labeled by within-timepoint quartile cut-points (type 7 quantiles):
  - **Low:** value ≤ Q1
  - **High:** value ≥ Q3
  - **Interquartile range:** everything in between
- **Filtering.** miRNAs are retained if they have at least 5 counts in at least 3 samples.
- **Normalization.** Variance-stabilizing transformation (`varianceStabilizingTransformation`, `blind = TRUE`) for PCA and expression plots.
- **Background list.** The filtered miRNA set is written to `DIANA_background_miRNAs.txt` for use as the background universe in pathway enrichment.

### 2.2 QC and descriptive figures (Figure 5)

- **5A:** RNA concentration (Qubit) per sample. Values are entered manually in the code.
- **5B:** Raw vs. DESeq2-normalized total counts per sample.
- **5C:** PCA of all filtered miRNAs (PC1 vs. PC2, labeled by sample).
- **5D:** Top 20 most highly expressed miRNAs (mean VST expression).

### 2.3 PCA and correlation with activity (Figure 6, Table S3)

- **PCA:** top 500 most variable miRNAs (VST), centered and scaled (`prcomp`).
- **Correlation grid:** PC1 and PC2 scores vs. each exposure at each timepoint, giving **16 tests** (2 activities × 4 timepoints × 2 PCs). Each test reports:
  - Pearson r and Spearman ρ with p-values
  - 5,000-iteration bootstrap 95% confidence intervals
  - Benjamini-Hochberg (BH) adjusted p-values, applied across all 16 tests separately for Pearson and Spearman
- **Regression plots (Figure 6A-D):** PC2 vs. SED60+ at each timepoint, points colored by quartile group, annotated with raw and BH-adjusted p-values. A combined 4-panel version is also saved.
- **Loading contribution diagnostic:** quantifies how much hsa-miR-148a-3p contributes to PC1 and PC2 relative to a uniform expectation, and plots the top 20 contributors to each component.

### 2.4 Differential expression with DESeq2 (Table 2, Table S5, Figure 7)

- **Design:** `~ Group` comparing **Low vs. High** quartile participants. Interquartile-range participants are excluded from the fit.
- **Per-comparison filtering:** normalized count ≥ 5 in at least 3 samples.
- **Significance:** miRNAs with `padj < 0.10` are written to the `*_DEGs_NAomit.csv` files.
- **Coverage:** 8 comparisons total (SED60+ and MVPA × Trimester 1 / 2 / 3 / postpartum).
- **Volcano plots (Figure 7A-C):** SED60+ in Trimesters 1-3 (EnhancedVolcano; labeled at `padj < 0.05`, |log2FC| > 1).

### 2.5 Sensitivity analyses

| Analysis | Question addressed | Key outputs |
|---|---|---|
| **PCA: quartile vs. continuous concordance** | Do the extreme-quartile PC score differences (Wilcoxon) agree in direction with the full-sample continuous correlations? | `PCA_Quartile_vs_Continuous_Concordance.csv` |
| **DE: continuous predictor (full n = 11)** | Do quartile-based DEGs replicate when every participant is used and the exposure is modeled continuously (z-scored, log2FC per 1 SD)? | `results_*_continuous_DEGs_NAomit.csv`, `Quartile_vs_Continuous_Concordance.csv` |
| **DE: full-dispersion 3-level group model** | Do the DEGs hold when dispersion and size factors are estimated on all 11 participants (Low / Mid / High), with the contrast still Low vs. High? | `results_*_fulldispersion_DEGs_NAomit.csv`, `Quartile_vs_FullDispersion_Comparison.csv` |
| **DE: exclude hsa-miR-148a-3p** | Are results robust to removing miR-148a-3p, which is highly abundant and may dominate variance? Evaluated at `padj < 0.05`. | `results_*_excl148a_DEGs_NAomit.csv`, `Sensitivity_Exclude_miR148a_AllTimepoints.csv` |
| **PCA: exclude hsa-miR-148a-3p** | Do the 16 PC-activity correlations keep the same direction when the PCA is rebuilt without miR-148a-3p? | `PCA_Sensitivity_Exclude_miR148a.csv` |

### 2.6 Supplementary analyses

- **Per-participant normalized counts** for every significant miRNA in Trimester 1 and Trimester 2 (long-format CSVs, a combined master table, and faceted dot plots).
- **Dispersion diagnostics:** `plotDispEsts` for the Trimester 1 and Trimester 2 models.
- **PC2 loading consistency:** for each significant miRNA, checks whether the direction of differential expression agrees with the direction implied by its PC2 loading and the sign of the PC2-SED60+ correlation. An exact binomial test (vs. 50%) is run per timepoint and BH-corrected across the two timepoints.

### 2.7 Pathway enrichment (Figure 8A)

Bubble plot of the top GO Biological Process terms enriched among targets of Trimester 2 SED60+-associated miRNAs (point size = fold enrichment; color = −log10 Bonferroni-adjusted p-value). Enrichment is run outside R using `DIANA_background_miRNAs.txt` as the background; results are read from `T2_sednum60_pathways.csv`.

### 2.8 Final outputs

The last chunk copies all deliverables into `Final_Outputs/<category>/` and writes `Final_Outputs/_manifest.csv` listing every expected file and whether it was found. Any missing files are printed to the console.

---

## Figures and tables: where to find them

| Item | Description | Generated file(s) | Folder |
|---|---|---|---|
| Figure 5A | Qubit RNA concentrations | `Figure_4A.tiff` | `QC/` |
| Figure 5B | Raw vs. normalized library sizes | `Figure_4B.tiff` | `QC/` |
| Figure 5C | PCA (PC1 vs. PC2) | `Figure_4C.tiff` | `PCA_and_Correlation/` |
| Figure 5D | Top 20 most expressed miRNAs | `Figure_4D.tiff` | `PCA_and_Correlation/` |
| Figure 6A-D | PC2 vs. SED60+ (T1, T2, T3, postpartum) | `Figure_6A.tiff` ... `Figure_6D.tiff`, `Figure_6_Combined.tiff` | `PCA_and_Correlation/` |
| Table S3 | PC1/PC2 vs. SED60+ and MVPA, all timepoints | `TableS3_PC_Correlations_Full16.csv`, `PC_Correlation_Summary_Full16.csv` | `PCA_and_Correlation/` |
| Table 2 | DE miRNAs, SED60+ (T1, T2, T3, postpartum) | `results_sedentary_*_DEGs_NAomit.csv` | `Differential_Expression/` |
| Table S5 | DE miRNAs, MVPA (T1, T2, T3, postpartum) | `results_mvpa_*_DEGs_NAomit.csv` | `Differential_Expression/` |
| Figure 7A-C | Volcano plots, SED60+ T1-T3 | `Figure_7A.tiff` ... `Figure_7C.tiff` | `Differential_Expression/` |
| Supplementary | Per-participant normalized counts (T1, T2) | `Supplementary_Table_Individual_Normalized_Counts_ALL.csv`, `Supplementary_IndividualCounts_T*.csv`, `Supp_Figure_IndividualCounts_T*.tiff` | `Differential_Expression/` |
| Supplementary | Dispersion estimate plots | `DispersionEstimates_Trimester1.tiff`, `DispersionEstimates_Trimester2.tiff` | `Differential_Expression/` |
| Supplementary | PC2 loading consistency | `PC2_Consistency_T1.csv`, `PC2_Consistency_T2.csv` | `PCA_and_Correlation/` |
| Supplementary | Loading contributions | `PCA_Loading_Contributions_AllGenes.csv`, `Figure_PCA_Contribution_PC1.tiff`, `Figure_PCA_Contribution_PC2.tiff` | `PCA_and_Correlation/` |
| Sensitivity | All sensitivity results | see [Section 2.5](#25-sensitivity-analyses) | `Sensitivity_Analyses/` |
| Figure 8A | GO enrichment bubble plot (T2 SED60+) | `Figure_8A.tiff` | `Pathway Analysis Files/` |
| Background list | miRNA universe for enrichment | `DIANA_background_miRNAs.txt` | `Pathway Analysis Files/` |

> **Note:** Panels 5A-5D are saved with `Figure_4*` filenames in the code. Panels B-D of the QC figure were assembled and lettered in PowerPoint for the manuscript.

---

## Sign conventions (read before interpreting log2FC)

Contrast direction differs between the SED60+ and MVPA models:

| Model | Reference | Comparison | Positive log2FC means... |
|---|---|---|---|
| **SED60+** (quartile) | High (Q3, more sedentary) | Low (Q1, less sedentary) | **Higher** in the **less** sedentary group |
| **MVPA** (quartile) | Low (Q1, less active) | High (Q3, more active) | **Higher** in the **more** active group |
| **Continuous** (either) | n/a | per +1 SD of the predictor | Higher with **increasing** predictor value |

Because the SED60+ quartile contrast is "Low vs. High", its sign is opposite to the continuous model's. The concordance code inverts the continuous log2FC for SED60+ so directions are compared in real-world terms. Each DEG file also includes an explicit `Direction` column stating which group the miRNA is higher in.

---

## Software and dependencies

**Command-line (Step 1):** FastQC, BBTools, Cutadapt, Bowtie 1, SAMtools.

**R (Step 2):** The `.Rmd` installs any missing packages on first run.

- Bioconductor: `DESeq2`, `vsn`, `EnhancedVolcano`
- CRAN: `tidyverse` (includes `ggplot2`, `dplyr`, `readr`), `pheatmap`, `RColorBrewer`, `matrixStats`, `patchwork`, `ggrepel`

Record your exact versions with `sessionInfo()` after running.

---

## How to reproduce

1. **Clone the repository**
   ```bash
   git clone https://github.com/bpereira2022/Maternal-sedentary-behavior-and-breast-milk-EV-miRNA-across-pregnancy.git
   cd Maternal-sedentary-behavior-and-breast-milk-EV-miRNA-across-pregnancy
   ```
2. **(Optional) Regenerate counts from raw FASTQ.** Edit the paths and sample list in `Raw_Data_Processing.sh`, run it, and merge the per-sample count tables into `Merged_Count.csv`. This step requires the raw FASTQ files, which are not stored in this repository.
3. **Run the analysis.** Open `Analysis_Plots_NOU_Activity.Rmd` in RStudio, set the working directory to the folder containing `Merged_Count.csv` and `ColData.csv`, and run all chunks (or knit).
4. **Pathway enrichment.** Use `DIANA_background_miRNAs.txt` as the background and place the enrichment results at `T2_sednum60_pathways.csv` before running the Figure 8A chunk.
5. **Collect outputs.** The final chunk organizes everything into `Final_Outputs/`, with a manifest of found and missing files.

---

## Citation and contact

If you use this code or data, please cite the associated manuscript:

> Elliott BB, *et al.* Maternal Sedentary Behavior and Breast Milk Extracellular Vesicle MicroRNA Across Pregnancy: A Pilot Study.

**Contact:** Beatriz B. Elliott, GitHub: [@bpereira2022](https://github.com/bpereira2022)
