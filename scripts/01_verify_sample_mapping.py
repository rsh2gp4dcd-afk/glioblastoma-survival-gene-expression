import pandas as pd
from pathlib import Path

# --------------------------------------------------
# GSE149009 sample mapping verification
# --------------------------------------------------

CLINICAL_FILE = Path("../data/GSE149009_clean_clinical_metadata.csv")
EXPRESSION_FILE = Path("../data/GSE149009_gbm_rna_deseq2.txt")

# Load clinical metadata
clinical = pd.read_csv(CLINICAL_FILE)

print("Clinical cohort:")
print(f"Number of patients: {len(clinical)}")

# Check that each patient has one GSM ID
assert clinical["gsm"].notna().all(), "Missing GSM IDs found."
assert clinical["gsm"].is_unique, "Duplicate GSM IDs found."

print("\nGSM IDs are present and unique.")

# --------------------------------------------------
# Load expression matrix
# --------------------------------------------------

expression = pd.read_csv(
    EXPRESSION_FILE,
    sep="\t",
    index_col=0
)

print("\nExpression matrix:")
print(f"Genes: {expression.shape[0]}")
print(f"Samples: {expression.shape[1]}")

expression_samples = set(expression.columns)
clinical_samples = set(clinical["expr_column"])

# --------------------------------------------------
# Verify clinical -> expression mapping
# --------------------------------------------------

missing_from_expression = clinical_samples - expression_samples
missing_from_clinical = expression_samples - clinical_samples

print("\nMapping check:")

if missing_from_expression:
    print("Clinical samples missing from expression matrix:")
    print(sorted(missing_from_expression))
else:
    print("All clinical expression-column IDs were found.")

if missing_from_clinical:
    print("\nExpression samples not represented in clinical cohort:")
    print(sorted(missing_from_clinical))
else:
    print("All expression samples are represented in the clinical cohort.")

# --------------------------------------------------
# Verify one-to-one mapping
# --------------------------------------------------

assert len(missing_from_expression) == 0, (
    "Some clinical samples are missing from the expression matrix."
)

# Reorder expression matrix to exactly match clinical cohort
ordered_samples = clinical["expr_column"].tolist()

expression_aligned = expression[ordered_samples]

assert list(expression_aligned.columns) == ordered_samples

print("\nSUCCESS:")
print("Expression matrix has been aligned to the clinical cohort.")
print(f"Aligned samples: {expression_aligned.shape[1]}")
print(f"Genes: {expression_aligned.shape[0]}")

# --------------------------------------------------
# Save mapping table
# --------------------------------------------------

mapping = clinical[
    ["index", "patient_id", "gsm", "expr_column"]
].copy()

output = Path("../results/mapping")
output.mkdir(parents=True, exist_ok=True)

mapping.to_csv(
    output / "verified_sample_mapping.csv",
    index=False
)

print("\nSaved:")
print("results/mapping/verified_sample_mapping.csv")
