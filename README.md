# GBM Gene Expression and Survival-Associated Pathway Analysis

## Overview

This project investigates associations between gene-expression profiles and recorded survival duration in glioblastoma (GBM), using publicly available data from the Gene Expression Omnibus (GEO) dataset GSE149009.

The analysis combines exploratory data analysis, gene-level correlation analysis, and Gene Ontology (GO) Biological Process pathway enrichment to investigate molecular processes associated with variation in recorded survival duration.

## Research Objectives

- Examine and verify the mapping between clinical metadata and gene-expression samples.
- Explore the gene-expression dataset and its clinical variables.
- Rank genes by their association with recorded survival duration using Spearman correlation.
- Perform preranked gene-set enrichment analysis to identify biological processes associated with the ranked gene list.
- Examine pathway-enrichment results under a protein-coding gene sensitivity analysis, where applicable.

## Data Source

**Dataset:** GSE149009  
**Repository:** [NCBI Gene Expression Omnibus](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149009)

The original dataset contains gene-expression profiles from glioblastoma samples alongside associated clinical information. The original study should be cited when using these data.

## Repository Structure

```text
gbm-survival-pathway-analysis/
├── README.md
├── 01_verify_sample_mapping.py
├── 02_expression_eda.py
├── GSE149009_Exploratory_Analysis.ipynb
├── GSE149009_Survival_Analysis.ipynb
├── GSE14909_Pathway_Analysis.ipynb
├── 01_full_gene_survival_rankings_30617.csv
├── 02_prerank_all_Ensembl_IDs.csv
├── 05_GO_Biological_Process_prerank_REBUILT.csv
├── 06_significant_GO_pathways_FDR05.csv
├── top_negative_pathways.jpeg
└── top_postive_pathways.jpeg
```

## Methods

### 1. Sample Mapping and Exploratory Analysis

Sample-mapping checks and exploratory analyses were used to examine the correspondence between clinical metadata and expression samples, as well as the structure and characteristics of the dataset.

### 2. Gene-Level Association Analysis

Spearman rank correlation was used to assess associations between gene-expression measurements and recorded survival duration. Multiple-testing correction was applied to account for the large number of genes tested.

### 3. Pathway Enrichment Analysis

Genes were ranked by their association statistics and analysed using preranked gene-set enrichment analysis with the GO Biological Process 2023 gene-set collection and GSEApy.

Pathway-level enrichment was evaluated using normalized enrichment scores (NES) and false discovery rate (FDR)-adjusted values.

### 4. Sensitivity Analysis

A protein-coding gene analysis was conducted to examine whether pathway-level patterns were broadly consistent when the gene ranking was restricted to protein-coding genes. The final manuscript should specify the exact filtering and identifier-mapping procedures used.

## Main Analysis Outputs

- `01_full_gene_survival_rankings_30617.csv`: gene-level association rankings.
- `02_prerank_all_Ensembl_IDs.csv`: ranked gene list used for pathway analysis.
- `05_GO_Biological_Process_prerank_REBUILT.csv`: rebuilt GO Biological Process enrichment results.
- `06_significant_GO_pathways_FDR05.csv`: pathway results meeting the specified FDR threshold.
- `top_negative_pathways.jpeg`: visualization of negatively enriched pathways.
- `top_postive_pathways.jpeg`: visualization of positively enriched pathways.

The exact relationship between the rebuilt enrichment results, the significant-pathway table, and the figures should be checked against the final manuscript before publication.

## Limitations

This analysis examines associations with recorded survival duration and should not be interpreted as demonstrating causation or establishing a validated prognostic biomarker.

The available clinical table contains missing survival-duration values and does not include a verified event/censoring indicator. Therefore, the analysis should not be described as a formal censored-survival analysis unless the required outcome information is independently verified.

Other limitations include potential effects of sample mapping, expression preprocessing, gene identifier mapping, duplicate gene symbols, and tied ranking statistics. Independent validation is required before drawing clinical conclusions.

## Reproducibility

Analysis scripts, notebooks, and derived result tables are provided in this repository to document the computational workflow.

Software versions, preprocessing choices, input-file requirements, and execution instructions should be recorded in the notebooks or accompanying documentation before publication.

The original expression data are available from GEO using the accession above.

## Citation

The original dataset is associated with the following publication:

Yanovich-Arad, G., et al. (2021). *Proteogenomics of glioblastoma associates molecular patterns with survival*. Cell Reports, 34(9), 108787. https://doi.org/10.1016/j.celrep.2021.108787

If this repository is used in future work, please cite the original study and the associated manuscript when a stable citation becomes available.
