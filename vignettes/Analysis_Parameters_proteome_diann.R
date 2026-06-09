# =============================================================================
# Analysis_Parameters.R
# Parameters for the proteomics (DIA-NN) workflow vignette.
# Adjust all paths and settings below before running the analysis.
# =============================================================================

# -----------------------------------------------------------------------------
# Input files
# -----------------------------------------------------------------------------

# Path to the DIA-NN Protein Group matrix .tsv report
PATH_TO_INPUT_Proteome <- "path/to/DIANN_protein_report.tsv"

# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------

# Directory where the Excel results file will be written
PATH_TO_OUTPUT <- "path/to/output/"

# -----------------------------------------------------------------------------
# Contaminant filtering
# -----------------------------------------------------------------------------

# Column in the DIA-NN Protein Group matrix that contains contaminant information
CONTAMINANT_COLUMN <- "Protein.Group"

# String pattern that identifies contaminant entries in CONTAMINANT_COLUMN
CONTAMINANT_FLAG <- "cRAP-"

# -----------------------------------------------------------------------------
# Missing value filtering
# -----------------------------------------------------------------------------

# Minimum number of measured (non-missing) values required per feature
# per condition to retain that feature
min_measured_val <- 2
