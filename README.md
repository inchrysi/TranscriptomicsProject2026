# TranscriptomicsProject2026

Course project for **Transcriptomics** (MSc in Bioinformatics for Computational Genomics, University of Milan and Politecnico di Milano). It contains two R Markdown analyses:

1. **Bulk RNA-seq differential expression** across three human tissues (brain, blood, lung) with edgeR.
2. **Single-cell RNA-seq analysis** of a colon dataset with Seurat, from quality control to cluster annotation.

## Repository content

| File | Description |
|---|---|
| `Transcriptomics_WFforDE.rmd` | Bulk RNA-seq workflow: sample QC, normalization, differential expression, tissue-specific genes |
| `ScsHW2.rmd` | scRNA-seq workflow: QC, normalization, dimensionality reduction, clustering, marker genes |

## 1. Bulk RNA-seq differential expression

**Data:** GTEx samples for brain, blood and lung, retrieved through [recount3](https://bioconductor.org/packages/recount3/) (3 replicates per tissue).

**Workflow**

- **Sample QC and replicate selection:** RIN score (≥ 6), percentage of uniquely mapped reads, fraction of rRNA reads. Replicates failing QC were replaced with new ones.
- **Filtering and normalization:** low-expression genes removed with `filterByExpr`, TMM normalization, log-CPM values compared before and after normalization.
- **Exploratory analysis:** MDS plots coloured by tissue and by metadata (rRNA fraction, mitochondrial reads, age, RIN) to check for technical effects.
- **Differential expression:** dispersion estimation, quasi-likelihood GLM fit (`glmQLFit`, `glmQLFTest`), pairwise contrasts (blood vs brain, blood vs lung, brain vs lung), Benjamini-Hochberg FDR < 0.05, with and without a log2 fold-change threshold.
- **Tissue-specific genes:** one-vs-others contrasts and intersection of pairwise results to select the top 100 genes per tissue, to be used for enrichment analysis (Enrichr).
- **Validation:** expression of selected genes checked across the full GTEx tissue set using TPM values and Wilcoxon tests.

**R packages:** `recount3`, `recount`, `edgeR`, `limma`, `ggplot2`, `ggrepel`, `ggforce`, `dplyr`, `tidyr`, `tibble`

## 2. Single-cell RNA-seq analysis

**Data:** sparse count matrix of a colon scRNA-seq dataset.

**Workflow**

- **QC and filtering:** number of genes and counts per cell, mitochondrial and ribosomal gene percentages; cells kept with 200-6000 detected genes and < 20% mitochondrial reads.
- **Normalization:** `LogNormalize`, cell cycle scoring, 2,000 highly variable genes (`vst`), scaling.
- **Dimensionality reduction:** PCA, choice of the number of components with elbow plot and explained variance, t-SNE and UMAP.
- **Clustering:** shared nearest neighbour graph and Louvain clustering, comparing different numbers of PCs and resolutions.
- **Marker genes:** Wilcoxon test with `FindMarkers` / `FindAllMarkers` (one cluster vs the rest, and cluster vs cluster comparisons).
- **Annotation:** manual annotation of clusters from marker genes (e.g. fibroblast subpopulations, endothelial cells, enteric glia, smooth muscle).

**R packages:** `Seurat`, `dplyr`, `patchwork`, `ggplot2`

## How to run

1. Install R and RStudio, then install the packages listed above.
2. Put the input data in the working directory (the data are not included in this repository):
   - bulk: `rse_brain.RDS`, `rse_blood.RDS`, `rse_lung.RDS`
   - single-cell: the sparse count matrix loaded as `sm`
3. Open the `.rmd` file in RStudio and knit it, or run the chunks one by one.
