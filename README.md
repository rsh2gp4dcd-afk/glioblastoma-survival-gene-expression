# GBM Survival–Pathway Analysis

## Overview

This project investigates associations between gene expression and survival duration in glioblastoma (GBM) samples using gene-level correlation analysis and pathway-level enrichment analysis.

The analysis uses data from the Gene Expression Omnibus (GEO) dataset GSE149009.

## Research Question

Which gene-expression patterns and biological pathways are associated with survival duration among the available GBM samples?

## Methods

1. Matched gene-expression data to clinical metadata using sample identifiers.
2. Excluded samples with missing survival information from gene–survival correlation analysis.
3. Calculated Spearman rank correlations between individual gene-expression values and survival days.
4. Applied Benjamini–Hochberg false-discovery-rate correction.
5. Ranked genes by signed Spearman correlation.
6. Mapped Ensembl gene identifiers to human gene symbols and constructed a unique-symbol ranking.
7. Conducted preranked Gene Set Enrichment Analysis (GSEA) using GO Biological Process gene sets.

The documented analysis settings use GSEApy prerank, 1,000 permutations, random seed 42, and a significance threshold of FDR q-value < 0.05. Package versions will be recorded from the final analysis environment.

## Results

- 68 samples were present in the source cohort.
- 59 samples had available survival days and were included in gene–survival correlations.
- 30,617 nonconstant genes were tested.
- No individual gene passed FDR < 0.05 in the reported gene-level analysis.
- 2,710 GO Biological Process pathways were tested.
- 66 pathways passed FDR q-value < 0.05.
- Of these significant pathways, 59 had positive normalized enrichment scores (NES) and 7 had negative NES.

Leading enriched terms included mitochondrial translation, mitochondrial gene expression, chromatin remodeling, and RNA processing and splicing.

## Interpretation

Positive NES indicates enrichment among genes positively correlated with survival; negative NES indicates enrichment among genes negatively correlated with survival, based on the ranking convention used in this analysis.

These results identify associations and hypotheses for further investigation. They do not establish that any gene or pathway causes longer or shorter survival.

## Limitations

- The gene–survival analysis used a relatively small sample of 59 cases.
- Nine samples lacked survival data and were excluded from the correlation analysis.
- The analysis does not establish causality.
- Potential clinical confounders were not fully adjusted for.
- Correlation with recorded survival days does not explicitly model censoring as a dedicated survival-analysis method would.
- Pathway findings require validation using independent data.

## Repository Contents

- `notebooks/`: analysis notebook.
- `results/gene_survival/`: gene-level association results and ranked gene lists.
- `results/pathway_analysis/`: complete GSEA results and significant pathway tables.
- `figures/`: pathway visualisations.
- `docs/`: detailed research documentation.

## Reproducibility

The analysis notebook should document the data source, preprocessing, package versions, random seed, gene-set collection, permutation count, and output filenames.

## Data Source

Gene Expression Omnibus (GEO), accession GSE149009. Add the verified original dataset citation before public release.
