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