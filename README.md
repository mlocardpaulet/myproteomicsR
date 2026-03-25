# myproteomicsR

R package for statistical analysis of DIA phosphoproteomics data from Spectronaut. Covers data loading, quality control, missing value imputation, limma-based differential abundance analysis, and Excel export.

## Installation

```r
# Install devtools if needed
install.packages("devtools")

# Install myproteomicsR from a local path
devtools::install("path/to/myproteomicsR")
```

If the source is on GitHub:

```r
# Install myproteomicsR
pak::pak("mlocardpaulet/myproteomicsR")
```

## Usage

See the vignette for a full worked example:

```r
library(myproteomicsR)
vignette("phosphoproteomics-workflow", package = "myproteomicsR")
```

## Dependencies

Imports: `data.table`, `dplyr`, `ggplot2`, `ggrepel`, `knitr`, `limma`, `magrittr`, `openxlsx`, `tibble`, `tidyr`
