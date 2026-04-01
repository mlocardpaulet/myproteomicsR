#' Convert a long-format quantification table to wide format
#'
#' Pivots \code{log2_quan} values from long format (one row per
#' feature-sample combination) to wide format (one row per feature, one column
#' per sample), preserving feature-level metadata columns.
#'
#' @param long_table A long-format data frame containing at least the columns
#'   named by \code{feature_name}, \code{filename}, and \code{log2_quan}.
#' @param feature_name Character. Name of the column that identifies features
#'   (e.g. \code{"protein_group_accessions"} or \code{"site_ID"}).
#'
#' @return A wide-format \code{tibble} with one row per feature. Sample
#'   \code{log2_quan} values appear as columns named after each
#'   \code{filename}, and feature-level metadata columns are appended.
#'
#' @importFrom dplyr select distinct left_join all_of
#' @importFrom tidyr pivot_wider
pivot_quan_to_wide <- function(long_table, feature_name) {
  sample_varying_cols <- c("filename", "condition", "replicate", "log2_quan",
                           "quantity", "nr_precursors_used_for_quantification")
  meta_cols <- setdiff(names(long_table), c(sample_varying_cols, feature_name))

  wide_table <- long_table %>%
    select(all_of(c(feature_name, "filename", "log2_quan"))) %>%
    pivot_wider(names_from = filename, values_from = log2_quan)

  if (length(meta_cols) > 0) {
    feature_meta <- long_table %>%
      select(all_of(c(feature_name, meta_cols))) %>%
      distinct()
    wide_table <- left_join(wide_table, feature_meta, by = feature_name)
  }

  wide_table
}

#' Write analysis results to a multi-sheet Excel workbook
#'
#' Assembles proteome and phosphoproteome analysis results into a named list
#' and writes each element as a separate sheet in a dated \code{.xlsx} file.
#'
#' @param unfiltered_table_proteome A data frame. Raw (unfiltered) proteome
#'   table.
#' @param filtered_table_proteome A data frame. Proteome table after quality
#'   filtering, used as input to the statistical analysis.
#' @param stat_results_proteome A data frame. Limma output for the proteome
#'   (as returned by \code{\link{run_limma}}).
#' @param mv_counts_proteome A data frame. Measured-value counts for the
#'   proteome (as returned by \code{\link{missing_values_count}}).
#' @param proteome_after_replacement A data frame. Proteome table after missing
#'   value imputation (as returned by \code{\link{replace_missing_values}}).
#' @param unfiltered_table_phosphoproteome A data frame. Raw (unfiltered)
#'   phosphoproteome table.
#' @param filtered_table_phosphoproteome A data frame. Phosphoproteome table
#'   after quality filtering.
#' @param stat_results_phosphoproteome A data frame. Limma output for the
#'   phosphoproteome (as returned by \code{\link{run_limma}}).
#' @param mv_counts_phosphoproteome A data frame. Measured-value counts for
#'   the phosphoproteome (as returned by \code{\link{missing_values_count}}).
#' @param phosphoproteome_after_replacement A data frame. Phosphoproteome table
#'   after missing value imputation (as returned by
#'   \code{\link{replace_missing_values}}).
#' @param file_name Character. Base name for the output file (without date
#'   prefix or \code{.xlsx} extension). Default: \code{"output"}.
#' @param path Character. Directory path where the file should be saved.
#'   Must end with \code{"/"}. Default: \code{"/"}.
#' @param wide_format Logical. If \code{TRUE}, the quantification tables (raw,
#'   filtered, and after missing-value replacement) are exported in wide format
#'   with one row per feature and one column per sample rather than in the
#'   default long format. The measured-value count tables and the limma output
#'   tables are not affected by this option. Default: \code{FALSE}.
#'
#' @return Invisibly returns \code{NULL}. Called for its side effect of writing
#'   an Excel file to disk.
#'
#' @importFrom openxlsx write.xlsx
#' @export
write_excel_output <- function(
    unfiltered_table_proteome,
    filtered_table_proteome,
    stat_results_proteome,
    mv_counts_proteome,
    proteome_after_replacement,
    unfiltered_table_phosphoproteome,
    filtered_table_phosphoproteome,
    stat_results_phosphoproteome,
    mv_counts_phosphoproteome,
    phosphoproteome_after_replacement,
    file_name = "output",
    path = "/",
    wide_format = FALSE) {

  if (wide_format) {
    unfiltered_table_proteome        <- pivot_quan_to_wide(unfiltered_table_proteome,        "protein_group_accessions")
    filtered_table_proteome          <- pivot_quan_to_wide(filtered_table_proteome,          "protein_group_accessions")
    proteome_after_replacement       <- pivot_quan_to_wide(proteome_after_replacement,       "protein_group_accessions")
    unfiltered_table_phosphoproteome <- pivot_quan_to_wide(unfiltered_table_phosphoproteome, "site_ID")
    filtered_table_phosphoproteome   <- pivot_quan_to_wide(filtered_table_phosphoproteome,   "site_ID")
    phosphoproteome_after_replacement <- pivot_quan_to_wide(phosphoproteome_after_replacement, "site_ID")
  }

  file.name <- paste0(path, gsub("-", "", Sys.Date()), "_", file_name, ".xlsx")
  cat("Write output in", file.name, "...\n")
  l.stat.results <- vector(mode = "list")
  l.stat.results$`Protein groups raw` <- unfiltered_table_proteome
  l.stat.results$`Protein groups before stat` <- filtered_table_proteome
  l.stat.results$`Protein groups MV` <- mv_counts_proteome
  l.stat.results$`Protein groups after MV replace` <- proteome_after_replacement
  l.stat.results$`Protein groups Limma ouptut` <- stat_results_proteome
  l.stat.results$`Psites raw` <- unfiltered_table_phosphoproteome
  l.stat.results$`Psites before stat` <- filtered_table_phosphoproteome
  l.stat.results$`Psites MV` <- mv_counts_phosphoproteome
  l.stat.results$`Psites after MV replace` <- phosphoproteome_after_replacement
  l.stat.results$`Psites Limma ouptut` <- stat_results_phosphoproteome
  write.xlsx(x = l.stat.results, file = file.name)

}
