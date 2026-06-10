#' Build a limma design matrix and contrasts
#'
#' Constructs a no-intercept model matrix from the conditions present in
#' \code{input_data} and builds a contrast matrix for the specified comparisons
#' using \code{\link[limma]{makeContrasts}}.
#'
#' @param input_data A long-format data frame containing at least the columns
#'   \code{condition}, \code{replicate}, and \code{filename} (one row per
#'   feature-run combination).
#' @param matrix_columns Character vector. Columns used to build the model
#'   matrix. Default: \code{c("condition", "replicate")}.
#' @param all_cond_pairs Character vector of contrasts to pass to
#'   \code{\link[limma]{makeContrasts}}, in the form \code{"condA-condB"}.
#'   Default: \code{"V1-NEG"}.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{design}{A numeric \code{matrix} — the model (design) matrix with
#'       filenames as row names and conditions as column names.}
#'     \item{contrasts}{A contrast matrix produced by
#'       \code{\link[limma]{makeContrasts}} covering all pairwise condition
#'       comparisons.}
#'   }
#'
#' @importFrom dplyr select distinct
#' @importFrom limma makeContrasts
#' @export
create_limma_design <- function(input_data, matrix_columns = c("condition", "replicate"), all_cond_pairs = NULL) {
  
  fcon <- input_data %>%
    select(c("filename", matrix_columns)) %>%
    distinct()
  
  fcon[matrix_columns] <- lapply(fcon[matrix_columns], function(x) if (!is.factor(x) && !is.numeric(x)) factor(x) else x)
  
  formula_terms <- paste(paste0("fcon$", matrix_columns), collapse = "+")
  mat <- model.matrix(as.formula(paste("~ 0 +", formula_terms)))
  rownames(mat) <- fcon$filename
  
  # Clean column names: remove all "fcon$" prefixes and any remaining artifacts
  colnames(mat) <- gsub("fcon\\$", "", colnames(mat))
  
  # Check if all_cond_pairs is NULL
  if (is.null(all_cond_pairs) || length(all_cond_pairs) == 0) {
    possible_contrasts <- combn(colnames(mat), 2, simplify = FALSE) %>%
      sapply(paste, collapse = "-")
    stop(
      "Error: 'all_cond_pairs' cannot be NULL or empty. ",
      "It should be a character vector of contrasts in the form 'condA-condB'. ",
      "Possible options based on your design matrix are:\n",
      paste("-", possible_contrasts, collapse = "\n")
    )
  }
  
  list(design = mat, contrasts = makeContrasts(contrasts = all_cond_pairs, levels = colnames(mat)))
}

#' Run a limma differential abundance analysis
#'
#' Fits a linear model with empirical Bayes moderation to the quantification
#' data, applies the supplied contrasts, and returns tidy top-table results for
#' every pairwise comparison augmented with per-condition measured-value counts.
#'
#' @param input_data A long-format data frame containing at least the columns
#'   \code{filename}, \code{log2_quan}, and the column referenced by
#'   \code{feature_name}.
#' @param design A numeric design matrix as returned by the
#'   \code{design} element of \code{\link{create_limma_design}}.
#' @param contrasts A contrast matrix as returned by the \code{contrasts}
#'   element of \code{\link{create_limma_design}}.
#' @param feature_name Character. Name of the column that identifies features
#'   (e.g. \code{"protein_group_accessions"} or \code{"site_ID"}).
#'   Must not be \code{NULL}.
#'
#' @return A \code{data.table} with one row per feature-contrast combination
#'   containing the standard limma \code{topTable} columns (\code{logFC},
#'   \code{AveExpr}, \code{t}, \code{P.Value}, \code{adj.P.Val}, \code{B})
#'   plus \code{contrast}, \code{feature}, \code{cond1}, \code{cond2},
#'   \code{count_val_cond1}, and \code{count_val_cond2}.
#'
#' @seealso \code{\link{create_limma_design}}
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
