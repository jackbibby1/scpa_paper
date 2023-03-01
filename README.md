# Systematic Single Cell Pathway Analysis (SCPA) to characterize early T cell activation
Repo contains all the R scripts to replicate figures of the paper, and bulk RNA-seq analysis of Ifnar1-/- mice.

Raw and processed sequencing data for both bulk and scRNA-seq can be found at [GSE212270](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE212270) and [SRA PRJNA874734](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA874734&o=acc_s%3Aa).

GSE212270 data contain:

- scRNA-seq data:
  - GSE212270_RAW.tar containing matrix, features, and barcode files for each sample
  - Seurat object .rds file for each sorted T cell population, integrated across the three time points (i.e. Naive CD4+ T cells at 0, 12, and 24 hours integrated, then Memory CD4+ T cells at 0, 12, and 24 hours integrated...)

- Bulk RNA-seq data
  - GSE212270_counts.txt.gz containing featureCounts output of Ifnar1-/- versus WT bulk RNA-sequencing
  - GSE212270_tmm_cpm_normalised.csv.gz containing TMM normalised counts per million values of edgeR processed Ifnar1-/- versus WT bulk RNA-sequencing
