#' Build a limma design matrix and contrasts for differential abundance analysis
#'
#' Constructs a model matrix from the variables specified in the \code{formula} and builds
#' a contrast matrix for the specified comparisons using \code{\link[limma]{makeContrasts}}.
#' This function is designed for experimental designs with categorical conditions (e.g., "BFG", "ASC").
#'
#' @param input_data A long-format data frame containing at least the columns:
#'   \itemize{
#'     \item{\code{filename}: Unique identifiers for each run/sample.}
#'     \item{\code{condition}: Experimental condition (e.g., "BFG", "ASC").}
#'     \item{Other variables included in \code{formula} (if any).}
#'   }
#'   One row per feature-run combination is expected.
#'
#' @param formula A formula specifying the design matrix. Default: \code{~ condition},
#'   which creates a design matrix for the \code{condition} variable only.
#'   Use standard R formula syntax to include additional terms or interactions. Examples:
#'   \itemize{
#'     \item{\code{~ condition}: Main effect of condition (default).}
#'     \item{\code{~ condition + BMI}: Main effects of condition and BMI.}
#'     \item{\code{~ condition * BMI}: Main effects + interaction between condition and BMI.}
#'   }
#'
#' @param all_cond_pairs A character vector of contrasts to pass to
#'   \code{\link[limma]{makeContrasts}}, in the form \code{"term1-term2"} or \code{"term"}.
#'   If \code{NULL}, the function will stop and suggest all possible pairwise contrasts
#'   based on the design matrix. Examples:
#'   \itemize{
#'     \item{\code{"BFG-ASC"}: Difference between BFG and ASC.}
#'     \item{\code{"BMI"}: Main effect of BMI (if included in the formula).}
#'     \item{\code{"conditionBFG:BMI"}: Interaction effect (if interaction is included in the formula).}
#'   }
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{design}{A numeric \code{matrix} — the model (design) matrix with filenames as row names
#'       and terms from the formula as column names. Column names are cleaned to be valid R names.}
#'     \item{contrasts}{A contrast matrix produced by \code{\link[limma]{makeContrasts}} for the
#'       specified \code{all_cond_pairs}.}
#'   }
#'
#' @section Design Matrix and Contrasts:
#' The design matrix is built using the provided \code{formula}, which allows for flexible
#' modeling of main effects and interactions. The function automatically:
#' \itemize{
#'   \item{Converts non-factor columns (except \code{filename}) to factors.}
#'   \item{Cleans column names to remove \code{fcon$} prefixes and ensure validity for \code{limma}.}
#'   \item{Generates contrasts for the specified comparisons.}
#' }
#' If \code{all_cond_pairs} is \code{NULL}, the function will suggest all possible pairwise contrasts
#' based on the design matrix.
#'
#' @section Examples:
#' \preformatted{
#' # Example 1: Default usage (condition only)
#' design_result <- create_limma_design(
#'   input_data = my_data,
#'   all_cond_pairs = c("BFG-ASC")
#' )
#'
#' # Example 2: Include BMI as a covariate
#' design_result <- create_limma_design(
#'   input_data = my_data,
#'   formula = ~ condition + BMI,
#'   all_cond_pairs = c("BFG-ASC", "BMI")
#' )
#'
#' # Example 3: With interaction between condition and BMI
#' design_result <- create_limma_design(
#'   input_data = my_data,
#'   formula = ~ condition * BMI,
#'   all_cond_pairs = c("BFG-ASC", "BMI.conditionBFG")
#' )
#'
#' # Example 4: Let the function suggest contrasts
#' tryCatch(
#'   design_result <- create_limma_design(
#'     input_data = my_data,
#'     formula = ~ condition + BMI
#'   ),
#'   error = function(e) print(e$message)
#' )
#' }
#'
#' @seealso
#' \code{\link[limma]{makeContrasts}} for details on contrast specifications.
#' \code{\link[limma]{model.matrix}} for details on design matrices.
#'
#' @importFrom dplyr select distinct
#' @importFrom limma makeContrasts
#' @export
create_limma_design <- function(
    input_data,
    formula = ~ condition,  # Default formula (no interaction)
    all_cond_pairs = NULL
) {
  
  # Extract unique combinations of filename and factors in the formula
  fcon <- input_data %>%
    select(c("filename", all.vars(formula))) %>%
    distinct()
  
  # Convert non-factor columns to factors (except filename)
  fcon[all.vars(formula)] <- lapply(
    fcon[all.vars(formula)],
    function(x) if (!is.factor(x) && !is.numeric(x)) factor(x) else x
  )
  
  # Build the model matrix using the provided formula
  mat <- model.matrix(formula, data = fcon)
  rownames(mat) <- fcon$filename
  
  # Clean column names: remove "fcon$" and make them valid R names
  colnames(mat) <- gsub("fcon\\$", "", colnames(mat))
  colnames(mat) <- make.names(colnames(mat), unique = TRUE)
  
  # Check if all_cond_pairs is NULL
  if (is.null(all_cond_pairs) || length(all_cond_pairs) == 0) {
    # Generate all possible pairwise contrasts from the design matrix
    possible_contrasts <- combn(colnames(mat), 2, simplify = FALSE) %>%
      sapply(paste, collapse = "-")
    stop(
      "Error: 'all_cond_pairs' cannot be NULL or empty. ",
      "It should be a character vector of contrasts in the form 'condA-condB'. ",
      "Possible options based on your design matrix are:\n",
      paste("-", possible_contrasts, collapse = "\n")
    )
  }
  
  # Create contrasts
  contrasts_mat <- makeContrasts(contrasts = all_cond_pairs, levels = colnames(mat))
  
  list(design = mat, contrasts = contrasts_mat)
}


#' Run limma differential abundance analysis
#'
#' Fits a linear model with empirical Bayes moderation to quantification data,
#' applies the supplied contrasts, and returns tidy top-table results for every
#' contrast. Results are augmented with per-condition measured-value counts.
#'
#' @param input_data A long-format data frame containing at least the columns:
#'   \itemize{
#'     \item{\code{filename}: Unique identifiers for each run/sample.}
#'     \item{\code{log2_quan}: Log2-transformed quantification values.}
#'     \item{A column specified by \code{feature_name} (e.g., \code{protein_group_accessions}).}
#'   }
#'   One row per feature-run combination is expected.
#'
#' @param design A numeric design matrix as returned by \code{\link{create_limma_design}}.
#'   This matrix should have filenames as row names and design terms as column names.
#'
#' @param contrasts A contrast matrix as returned by \code{\link{create_limma_design}}.
#'   This matrix defines the comparisons of interest (e.g., "BFG-ASC").
#'
#' @param feature_name Character. Name of the column in \code{input_data} that identifies
#'   features (e.g., \code{"protein_group_accessions"} or \code{"site_ID"}). Must not be \code{NULL}.
#'
#' @return A \code{data.table} with one row per feature-contrast combination containing:
#'   \describe{
#'     \item{Standard limma topTable columns:}{ \code{logFC}, \code{AveExpr}, \code{t},
#'       \code{P.Value}, \code{adj.P.Val}, \code{B}.}
#'     \item{contrast:}{The contrast name (e.g., "BFG-ASC").}
#'     \item{feature:}{The feature identifier (e.g., protein accession).}
#'     \item{cond1, cond2:}{For pairwise contrasts, the two conditions being compared.}
#'     \item{count_val_cond1, count_val_cond2:}{For pairwise contrasts, the number of non-NA
#'       values for each condition.}
#'     \item{cond, count_val_cond:}{For single-term contrasts, the condition and count of non-NA values.}
#'   }
#'
#' @section Details:
#' This function performs the following steps:
#' \itemize{
#'   \item{Reshapes \code{input_data} into a wide-format matrix suitable for \code{limma}.}
#'   \item{Fits a linear model using \code{\link[limma]{lmFit}}.}
#'   \item{Applies contrasts using \code{\link[limma]{contrasts.fit}}.}
#'   \item{Prints an MA plot for diagnostic purposes using \code{\link[limma]{plotMA}}.}
#'   \item{Applies empirical Bayes moderation using \code{\link[limma]{eBayes}} with \code{robust = TRUE}.}
#'   \item{Extracts top tables for each contrast and augments them with per-condition counts.}
#' }
#'
#' The function handles both pairwise contrasts (e.g., "BFG-ASC") and single-term contrasts
#' (e.g., "BMI"). For pairwise contrasts, it adds columns for the two conditions and their
#' respective counts of non-NA values. For single-term contrasts, it adds a column for the
#' condition and its count of non-NA values.
#'
#' @section Examples:
#' \preformatted{
#' # Example 1: Basic usage with pairwise contrast
#' design_result <- create_limma_design(
#'   input_data = my_data,
#'   all_cond_pairs = c("BFG-ASC")
#' )
#' limma_results <- run_limma(
#'   input_data = my_data,
#'   design = design_result$design,
#'   contrasts = design_result$contrasts,
#'   feature_name = "protein_group_accessions"
#' )
#'
#' # Example 2: With a single-term contrast (e.g., BMI)
#' design_result <- create_limma_design(
#'   input_data = my_data,
#'   formula = ~ condition + BMI,
#'   all_cond_pairs = c("BMI")
#' )
#' limma_results <- run_limma(
#'   input_data = my_data,
#'   design = design_result$design,
#'   contrasts = design_result$contrasts,
#'   feature_name = "protein_group_accessions"
#' )
#' }
#'
#' @seealso
#' \code{\link{create_limma_design}} for creating the design and contrast matrices.
#' \code{\link[limma]{lmFit}}, \code{\link[limma]{contrasts.fit}}, \code{\link[limma]{eBayes}},
#' \code{\link[limma]{topTable}} for details on the underlying limma functions.
#'
#' @importFrom dplyr mutate select group_by ungroup
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames
#' @importFrom limma lmFit contrasts.fit eBayes topTable plotMA
#' @importFrom data.table rbindlist
#' @export
run_limma <- function(input_data, design, contrasts, feature_name = NULL) {
  
  if (is.null(feature_name)) {
    stop("Please provide a feature name for the analysis: it should be the name of the corresponding column.\n")
  }
  
  input_data %>% 
    mutate(feature = get(feature_name)) %>% 
    select(feature, filename, log2_quan) %>%
    group_by(feature) %>%
    pivot_wider(names_from = filename, values_from = log2_quan) %>%
    ungroup() %>% 
    column_to_rownames("feature") -> limma_input
  
  fit <- lmFit(limma_input, design)
  fit <- contrasts.fit(fit, contrasts)
  print(plotMA(fit))
  fit <- eBayes(fit, robust = TRUE)
  
  lapply(
    colnames(contrasts), 
    function(x) {
      tt <- topTable(fit, coef=x, adjust="BH", number = nrow(limma_input))
      tt$contrast = x
      tt$feature = rownames(tt)
      
      if (!grepl("-", x)) {
        runs = rownames(design)[design[, x] > 0]
        
        tt %>% 
          mutate(
            cond = x,
            count_val_cond = rowSums(!is.na(limma_input[, runs]))[match(feature, rownames(limma_input))]
          )
      } else {
        ## count values per condtion
        cond1 = gsub("-.+", "", x)
        cond2 = gsub(".+-", "", x)
        runs1 = rownames(design)[design[, cond1] == 1]
        runs2 = rownames(design)[design[, cond2] == 1]
        
        tt %>% 
          mutate(
            cond1 = cond1,
            cond2 = cond2,
            count_val_cond1 = rowSums(!is.na(limma_input[, runs1]))[match(feature, rownames(limma_input))],
            count_val_cond2 = rowSums(!is.na(limma_input[, runs2]))[match(feature, rownames(limma_input))]
          )
      }
    }) %>% 
    rbindlist()
  
}
