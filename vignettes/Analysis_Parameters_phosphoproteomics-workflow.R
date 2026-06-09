# =============================================================================
# Analysis_Parameters.R
# Parameters for the phosphoproteomics workflow vignette.
# Adjust all paths and settings below before running the analysis.
# =============================================================================

# -----------------------------------------------------------------------------
# Input files
# -----------------------------------------------------------------------------

# Path to the Spectronaut Protein Group (MLP) .tsv report
PATH_TO_INPUT_Proteome <- "path/to/spectronaut_protein_report.tsv"

# Path to the Spectronaut PTM Site BGS .tsv report
PATH_TO_INPUT_PTM <- "path/to/spectronaut_ptm_report.tsv"

# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------

# Directory where the Excel results file will be written
PATH_TO_OUTPUT <- "path/to/output/"

# -----------------------------------------------------------------------------
# Contaminant filtering
# -----------------------------------------------------------------------------

# Column in the Spectronaut report that contains contaminant information
CONTAMINANT_COLUMN <- "PG.FastaHeaders"

# String pattern that identifies contaminant entries in CONTAMINANT_COLUMN
CONTAMINANT_FLAG <- "Cont_"

# -----------------------------------------------------------------------------
# Missing value filtering
# -----------------------------------------------------------------------------

# Minimum number of measured (non-missing) values required per feature
# per condition to retain that feature
min_measured_val <- 2
