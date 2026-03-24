#' Count measured values per feature and condition
#'
#' Summarises the number of non-NA log2 quantification values for each feature
#' in each condition, then returns a wide-format table with per-row maximum and
#' minimum measured value counts across conditions.
#'
#' @param input_table A data frame (or data.table) in long format containing at
#'   least the columns \code{condition}, \code{log2_quan}, and the column
#'   referenced by \code{feature_name}.
#' @param feature_name Character. Name of the column that identifies features
#'   (e.g. \code{"protein_group_accessions"} or \code{"site_ID"}).
#'   Must not be \code{NULL}.
#' @param title Character. Optional title string (currently unused, kept for
#'   forward compatibility). Default: \code{""}.
#'
#' @return A wide-format \code{tibble} with one row per feature. Columns are
#'   the feature identifier, one column per condition containing the count of
#'   measured values, plus \code{max_measured_val} and \code{min_measured_val}.
#'
#' @importFrom dplyr mutate summarise ungroup
#' @importFrom tidyr pivot_wider
#' @export
missing_values_count <- function(input_table, feature_name = NULL, title = "") {
  
  if (is.null(feature_name)) {
    stop("Please provide a feature name for the plot: it should be the name of the corresponding column.\n")
  }
  
  res <- input_table %>% 
    mutate(feature = get(feature_name)) %>%
    summarise(measured_values_count = sum(!is.na(log2_quan)), .by = c(condition, feature)) %>% 
    ungroup()
  return(
    res %>% 
      pivot_wider(names_from = condition, values_from = measured_values_count, values_fill = 0) %>% 
      mutate(max_measured_val = sapply(1:nrow(.), function(i) max(.[i, -1], na.rm = TRUE)),
             min_measured_val = sapply(1:nrow(.), function(i) min(.[i, -1], na.rm = TRUE))) 
  )
  
}

#' Filter features by minimum measured value count
#'
#' Retains only features that have at least \code{min_measured_val} non-NA
#' log2 quantification values in at least one condition (evaluated against the
#' per-row minimum across conditions).
#'
#' @param input_table A long-format data frame containing at least the columns
#'   \code{condition}, \code{log2_quan}, and the column referenced by
#'   \code{feature_name}.
#' @param feature_name Character. Name of the column that identifies features.
#'   Must not be \code{NULL}.
#' @param title Character. Plot title forwarded to
#'   \code{\link{missing_values_count}}. Default: \code{""}.
#' @param min_measured_val Integer. Minimum number of measured (non-NA) values
#'   required. Default: \code{2}.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{filtered_table}{The filtered long-format data frame.}
#'     \item{value_count}{The wide-format value-count table returned by
#'       \code{\link{missing_values_count}}.}
#'   }
#'
#' @seealso \code{\link{missing_values_count}}
#' @importFrom dplyr filter
#' @export
filter_value_count <- function(input_table, feature_name = NULL, min_measured_val = 0, title = "") {
  
  if (is.null(feature_name)) {
    stop("Please provide a feature name for the plot: it should be the name of the corresponding column.\n")
  }
  
  cat("keep only the features wit", min_measured_val, "or more measured values in at least one condition.\n")
  val_count <- missing_values_count(input_table, feature_name = feature_name, title = title)
  tokeep <- val_count$feature[val_count$max_measured_val >= min_measured_val] #TODO: fix it must be the sum of measured values across conditions, not the minimum
  cat("Keep", length(tokeep), "features on the", length(unique(unlist(input_table[,names(input_table) == feature_name, with = F]))), "rows of the input table.\\n")
  
  res <- input_table %>% 
    filter(get(feature_name) %in% tokeep)
  
  lres <- list(
    filtered_table = res,
    value_count = val_count
  )
  
  return(lres)
  
}

#' Compute per-feature coefficients of variation per condition
#'
#' Calculates the coefficient of variation (standard deviation divided by mean)
#' of log2 quantification values for each feature within each condition and
#' returns the result in wide format.
#'
#' @param input_table A long-format data frame containing at least the columns
#'   \code{condition}, \code{log2_quan}, and the column referenced by
#'   \code{feature_name}.
#' @param feature_name Character. Name of the column that identifies features.
#'   Must not be \code{NULL}.
#' @param title Character. Optional title string (currently unused, kept for
#'   forward compatibility). Default: \code{""}.
#'
#' @return A wide-format \code{tibble} with one row per feature and one column
#'   per condition containing the coefficient of variation.
#'
#' @importFrom dplyr mutate summarise ungroup
#' @importFrom tidyr pivot_wider
#' @export
feature_cvs <- function(input_table, feature_name = NULL, title = "") {
  
  if (is.null(feature_name)) {
    stop("Please provide a feature name for the plot: it should be the name of the corresponding column.\n")
  }
  
  input_table %>% 
    mutate(feature = get(feature_name)) %>%
    summarise(cv_per_feature = sd(log2_quan, na.rm = T)/mean(log2_quan, na.rm = T), .by = c(condition, feature)) %>% 
    ungroup() %>% 
    pivot_wider(names_from = condition, values_from = cv_per_feature, values_fill = 0) 
  
}
