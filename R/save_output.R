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
    path = "/") {
  
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
