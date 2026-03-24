#' Replace missing values for features absent in exactly one condition
#'
#' For each condition in which a feature has no measured values at all, this
#' function imputes the missing values using the per-run quantile of all
#' measured log2 quantities (controlled by \code{quantile_replacement_val}).
#' Features that are missing across all conditions are left unchanged.
#' Only designed for experiments with exactly two conditions.
#'
#' @param input_table A long-format data frame containing at least the columns
#'   \code{condition}, \code{filename}, \code{replicate}, \code{quantity},
#'   \code{log2_quan}, and the column referenced by \code{feature_name}.
#' @param value_count A wide-format measured-value count table as returned by
#'   \code{\link{missing_values_count}}, with a \code{feature} column and one
#'   column per condition.
#' @param quantile_replacement_val Numeric in \code{[0, 1]}. Quantile of the
#'   per-run measured log2 values used as the replacement value. For example,
#'   \code{0.01} uses the 1st percentile.
#' @param feature_name Character. Name of the column that identifies features
#'   (e.g. \code{"protein_group_accessions"} or \code{"site_ID"}).
#'   Must not be \code{NULL}.
#'
#' @return The input table augmented with imputed rows for features that were
#'   entirely missing in one condition.
#'
#' @seealso \code{\link{missing_values_count}}, \code{\link{wrap_replace_missing_values}}
#' @importFrom dplyr select distinct filter mutate pull
#' @importFrom knitr kable
#' @importFrom data.table rbindlist
#' @export
replace_missing_values <- function(input_table, 
                                   value_count,
                                   quantile_replacement_val, 
                                   feature_name = NULL) {
  
  if (is.null(feature_name)) {
    stop("Please provide a feature name for the plot: it should be the name of the corresponding column.\n")
  }
  
  # replacement value:
  toreplace <- sapply(unique(input_table$filename), function(x) {
    quantile(input_table$log2_quan[input_table$filename == x], probs = quantile_replacement_val, na.rm = TRUE)
  })
  names(toreplace) <- gsub("\\.\\d\\%$", "", names(toreplace))
  cat("Values used as replacement when missing: \n")
  print(toreplace)
  
  ucond <- names(value_count)[!(names(value_count) %in% c("feature", "max_measured_val", "min_measured_val"))]
  
  if (length(ucond) > 2) {
    break("More than 2 conditions found in the value count table. Please check the input table and the value count table.\n")
  }
  
  lfeat <- lapply(ucond, function(c) {
    value_count %>% 
      filter(.[[c]] == 0) %>% 
      pull(feature)
  })
  names(lfeat) <- ucond
  
  maptab <- input_table %>% 
    select(condition, filename, replicate) %>% 
    distinct()
  
  cat("Experimental plan:\n")
  print(kable(maptab))
  
  cat("Replace missing values for the following number of features:\n")
  print(sapply(lfeat, length))
  
  tab_to_add <- lapply(seq_along(lfeat), function(i) {
    f <- lfeat[[i]]
    lapply(f, function(s) {
      input_table[input_table[[feature_name]] == s & input_table$condition == names(lfeat)[setdiff(1:2, i)],] %>% 
        mutate(quantity = NA,
               filename = maptab$filename[match(paste(names(lfeat)[i], replicate), paste(maptab$condition, maptab$replicate))],
               condition = names(lfeat)[i],
               log2_quan = toreplace[match(filename, names(toreplace))])
    }) %>% 
      rbindlist()
  }) %>% 
    rbindlist()
  
  rbind(input_table, tab_to_add)
  
}

#' Apply missing-value replacement to proteome and phosphoproteome tables
#'
#' Iterates over the filtered proteome and phosphoproteome tables produced by
#' \code{\link{filter_value_count}} and calls \code{\link{replace_missing_values}}
#' on each, using the 1st percentile of measured values per run as the
#' imputation value.
#'
#' @param ltab A named list with elements \code{proteome} and \code{ptm},
#'   each being a list with \code{filtered_table} and \code{value_count}
#'   elements as returned by \code{\link{filter_value_count}}.
#' @param quantile_replacement_val Numeric in \code{[0, 1]}. Quantile of the
#'   per-run measured log2 values used as the replacement value.
#'   Default: \code{0.01} (1st percentile).
#'
#' @return A named list with the same names as \code{ltab}, where each element
#'   is the imputed long-format data frame returned by
#'   \code{\link{replace_missing_values}}.
#'
#' @seealso \code{\link{replace_missing_values}}, \code{\link{filter_value_count}}
#' @export
wrap_replace_missing_values <- function(ltab, quantile_replacement_val = 0.01) {
  
  lres <- lapply(seq_along(ltab), function(i) {
    cat("\n", names(ltab[i]), ":\n")
    replace_missing_values(ltab[[i]]$filtered_table, 
                           ltab[[i]]$value_count,
                           quantile_replacement_val, 
                           feature_name = ifelse(names(ltab)[i] == "proteome", "protein_group_accessions", "site_ID"))
    
  })
  
  names(lres) <- names(ltab)
  return(lres)
  
}
