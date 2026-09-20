import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler

# ============================================================
# GSE149009 — Gene Expression Exploratory Data Analysis
# ============================================================

# File locations
CLINICAL_FILE = Path("../data/GSE149009_clean_clinical_metadata.csv")
EXPRESSION_FILE = Path("../data/GSE149009_gbm_rna_deseq2.txt")

# Output folders
FIGURES = Path("../figures")
RESULTS = Path("../results/expression_eda")

FIGURES.mkdir(parents=True, exist_ok=True)
RESULTS.mkdir(parents=True, exist_ok=True)


# ============================================================
# 1. Load data
# ============================================================

print("Loading clinical metadata...")

clinical = pd.read_csv(CLINICAL_FILE)

print(f"Clinical samples: {len(clinical)}")

print("Loading expression matrix...")

expression = pd.read_csv(
    EXPRESSION_FILE,
    sep="\t",
    index_col=0
)

print(f"Genes: {expression.shape[0]}")
print(f"Samples: {expression.shape[1]}")


# ============================================================
# 2. Align expression samples with clinical metadata
# ============================================================

samples = clinical["expr_column"].tolist()

missing_samples = set(samples) - set(expression.columns)

if missing_samples:
    raise ValueError(
        f"These clinical samples are missing from the expression matrix: "
        f"{sorted(missing_samples)}"
    )

# Put expression samples in the exact same order as the clinical data
expression = expression[samples]

print("Sample alignment verified.")
print(f"Aligned samples: {expression.shape[1]}")


# ============================================================
# 3. Basic expression-data quality control
# ============================================================

print("\nRunning basic QC...")

total_missing = int(expression.isna().sum().sum())
total_zero = int((expression == 0).sum().sum())

qc_summary = pd.DataFrame({
    "number_of_genes": [expression.shape[0]],
    "number_of_samples": [expression.shape[1]],
    "total_missing_values": [total_missing],
    "total_zero_values": [total_zero]
})

qc_summary.to_csv(
    RESULTS / "expression_qc_summary.csv",
    index=False
)

print(f"Missing values: {total_missing}")
print(f"Zero values: {total_zero}")


# ============================================================
# 4. Expression distribution
# ============================================================

print("\nCreating expression distribution plot...")

expression_values = expression.to_numpy().flatten()

plt.figure(figsize=(8, 5))

plt.hist(
    expression_values,
    bins=100
)

plt.xlabel("Expression value")
plt.ylabel("Frequency")
plt.title("Distribution of Gene-Expression Values")

plt.tight_layout()

plt.savefig(
    FIGURES / "expression_distribution.png",
    dpi=300
)

plt.close()


# ============================================================
# 5. Mean expression by sample
# ============================================================

print("Calculating mean expression per sample...")

sample_means = expression.mean(axis=0)

sample_summary = pd.DataFrame({
    "sample": sample_means.index,
    "mean_expression": sample_means.values
})

sample_summary.to_csv(
    RESULTS / "sample_expression_summary.csv",
    index=False
)

plt.figure(figsize=(10, 5))

plt.bar(
    range(len(sample_means)),
    sample_means.values
)

plt.xlabel("Sample")
plt.ylabel("Mean expression")
plt.title("Mean Gene Expression by Sample")

plt.tight_layout()

plt.savefig(
    FIGURES / "mean_expression_by_sample.png",
    dpi=300
)

plt.close()


# ============================================================
# 6. Identify highly variable genes
# ============================================================

print("Identifying highly variable genes...")

gene_variance = expression.var(axis=1)

gene_variance = gene_variance.sort_values(
    ascending=False
)

top_1000 = gene_variance.head(1000)

top_1000.to_csv(
    RESULTS / "top_1000_variable_genes.csv",
    header=["variance"]
)

print("Top 1,000 variable genes saved.")


# ============================================================
# 7. PCA
# ============================================================

print("Running PCA...")

# Use the top 1,000 most variable genes
pca_data = expression.loc[top_1000.index].T

# Standardize genes before PCA
scaler = StandardScaler()

scaled_data = scaler.fit_transform(
    pca_data
)

# Calculate the first two principal components
pca = PCA(
    n_components=2
)

principal_components = pca.fit_transform(
    scaled_data
)


# ============================================================
# 8. Save PCA coordinates
# ============================================================

pca_results = pd.DataFrame({
    "sample": samples,
    "PC1": principal_components[:, 0],
    "PC2": principal_components[:, 1]
})

pca_results.to_csv(
    RESULTS / "pca_coordinates.csv",
    index=False
)


# ============================================================
# 9. PCA plot
# ============================================================

plt.figure(figsize=(7, 6))

plt.scatter(
    pca_results["PC1"],
    pca_results["PC2"]
)

pc1_percent = pca.explained_variance_ratio_[0] * 100
pc2_percent = pca.explained_variance_ratio_[1] * 100

plt.xlabel(f"PC1 ({pc1_percent:.1f}% variance)")
plt.ylabel(f"PC2 ({pc2_percent:.1f}% variance)")

plt.title(
    "PCA of the 1,000 Most Variable Genes"
)

plt.tight_layout()

plt.savefig(
    FIGURES / "PCA.png",
    dpi=300
)

plt.close()


# ============================================================
# 10. Finish
# ============================================================

print("\n======================================")
print("Expression EDA completed successfully")
print("======================================")

print(f"Genes analyzed: {expression.shape[0]}")
print(f"Samples analyzed: {expression.shape[1]}")

print("\nResults saved to:")
print("results/expression_eda/")

print("\nFigures saved to:")
print("figures/")
