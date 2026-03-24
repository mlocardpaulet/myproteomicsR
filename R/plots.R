#' Boxplot of log2 quantification values by run
#'
#' Produces a \code{ggplot2} boxplot showing the distribution of log2
#' quantification values for each run, coloured by condition.
#'
#' @param input_table A long-format data frame containing at least
#'   \code{filename}, \code{log2_quan}, and \code{condition} columns.
#' @param title Character. Plot title. Default: \code{""}.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom ggplot2 ggplot aes geom_boxplot theme_minimal scale_fill_viridis_d labs theme element_text
#' @export
plot_boxplot <- function(input_table, title = "") {
  
  p <- ggplot(input_table, aes(x = filename, y = log2_quan, fill = condition)) +
    geom_boxplot(alpha = 0.8) +
    theme_minimal() +
    scale_fill_viridis_d() +
    labs(y = "log2(Quantity)",
         x = "",
         title = title) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
  return(p)
  
}

#' Bar chart of missing values per run
#'
#' Computes the percentage of missing (\code{NA}) log2 quantification values
#' for each run and displays them as a bar chart coloured by condition.
#'
#' @param input_table A long-format data frame containing at least
#'   \code{filename}, \code{condition}, and \code{log2_quan} columns.
#' @param title Character. Plot title. Default: \code{""}.
#'
#' @return A \code{ggplot} object.
#'
#' @importFrom dplyr group_by summarise mutate
#' @importFrom ggplot2 ggplot aes geom_bar scale_fill_viridis_d labs theme_minimal theme element_text
#' @export
plot_missing_values_per_run <- function(input_table, title = "") {
  
  input_table %>% 
    group_by(filename, condition) %>% 
    summarise(proportion_missing_values = sum(is.na(log2_quan))/length(log2_quan)) %>% 
    mutate(perc_missing_values = round(proportion_missing_values*100, 2)) %>%
    ggplot(aes(
      x = filename,
      fill = condition,
      y = perc_missing_values
    )) +
    geom_bar(stat = "identity", col = "black") +
    scale_fill_viridis_d() +
    labs(y = "Percentage of missing values per run",
         x = "",
         title = title) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) 
  
}

#' PCA scatter plots (PC1/PC2 and PC3/PC4)
#'
#' Performs a principal component analysis on the complete-cases subset of the
#' feature matrix (rows without any missing values) and returns two scatter
#' plots: PC1 vs PC2 and PC3 vs PC4, with samples coloured by condition and
#' labelled by replicate.
#'
#' @param input_table A long-format data frame containing at least
#'   \code{filename}, \code{condition}, \code{replicate}, \code{log2_quan}, and
#'   the column referenced by \code{feature_name}.
#' @param title Character. Plot title. Default: \code{""}.
#' @param feature_name Character. Name of the column that identifies features.
#'   Must not be \code{NULL}.
#'
#' @return A list of two \code{ggplot} objects: \code{[[1]]} for PC1 vs PC2 and
#'   \code{[[2]]} for PC3 vs PC4.
#'
#' @importFrom tidyr pivot_wider
#' @importFrom tibble column_to_rownames rownames_to_column
#' @importFrom dplyr select left_join distinct
#' @importFrom ggplot2 ggplot aes geom_point scale_color_viridis_d labs theme_minimal
#' @importFrom ggrepel geom_text_repel
#' @export
plot_PCA <- function(input_table, title = "", feature_name = NULL) {
  
  if (is.null(feature_name)) {
    stop("Please provide a feature name for the plot: it should be the name of the corresponding column.\n")
  }
  
  names(input_table)[names(input_table) == feature_name] <- "feature"
  
  input_table %>% 
    select(feature, filename, log2_quan) %>% 
    pivot_wider(names_from = filename, values_from = log2_quan) %>% 
    column_to_rownames("feature") -> pcatab
  
  cat(sum(rowSums(is.na(pcatab)) == 0), "rows do not have missing values and will be used for PCA (on the total of", nrow(pcatab), "rows).\n")
  
  pcatab <- pcatab[rowSums(is.na(pcatab)) == 0, ]
  
  pres <- prcomp(t(pcatab), scale. = FALSE)
  cat("Perform PCA with centering but without scaling the data...\n")
  print(summary(pres))
  
  gtab <- as.data.frame(pres$x) %>% 
    rownames_to_column("filename") %>% 
    left_join(input_table %>% select(filename, condition, replicate) %>% distinct(), by = "filename") 
  
  g1 <- gtab %>% 
    ggplot(aes(x = PC1,
                   y = PC2)) +
    geom_point(aes(color = condition), size = 3, alpha = 0.8) +
    scale_color_viridis_d() +
    geom_text_repel(aes(label = replicate), size = 3, max.overlaps = 20) +
    labs(x = paste0("PC1 (", round(summary(pres)$importance[2,1]*100, 2), "%)"),
         y = paste0("PC2 (", round(summary(pres)$importance[2,2]*100, 2), "%)"),
         title = title) +
    theme_minimal() 
  
  g2 <- gtab %>% 
    ggplot(aes(x = PC3,
                   y = PC4)) +
    geom_point(aes(color = condition), size = 3, alpha = 0.8) +
    scale_color_viridis_d() +
    geom_text_repel(aes(label = replicate), size = 3, max.overlaps = 20) +
    labs(x = paste0("PC3 (", round(summary(pres)$importance[2,3]*100, 2), "%)"),
         y = paste0("PC4 (", round(summary(pres)$importance[2,4]*100, 2), "%)"),
         title = title) +
    theme_minimal() 
  
  return(list(g1, g2))
}

#' Volcano plot of limma results
#'
#' Plots log2 fold-change against -log10 p-value for each contrast, with
#' points coloured by significance status at the specified adjusted p-value
#' threshold.
#'
#' @param limma_results A data frame of limma results as returned by
#'   \code{\link{run_limma}}, containing at least \code{logFC},
#'   \code{P.Value}, \code{adj.P.Val}, and \code{contrast} columns.
#' @param title Character. Plot title. Default: \code{""}.
#' @param adj.p.value.threshold Numeric. Adjusted p-value threshold for
#'   colouring significant hits. Default: \code{0.05}.
#'
#' @return A \code{ggplot} object faceted by contrast.
#'
#' @seealso \code{\link{run_limma}}
#' @importFrom ggplot2 ggplot aes geom_hline geom_vline geom_point labs scale_color_manual facet_grid theme_minimal
#' @export
plot_volcano <- function(limma_results, title = "", adj.p.value.threshold = 0.05) {
  
  limma_results %>% 
    ggplot(aes(x = logFC, 
               y = -log10(P.Value),
               col = adj.P.Val <= adj.p.value.threshold)) +
    geom_hline(yintercept = 0) +
    geom_vline(xintercept = 0) +
    geom_point(alpha = 0.3) +
    labs(x = "log2(Fold Change)",
         y = "-log10(p-value)",
         title = title) +
    scale_color_manual(values = c("black", "green"), labels = c("Not significant", "Significant"), name = "") +
    facet_grid(contrast ~ .) +
    theme_minimal() 
  
}

#' Histogram of measured value counts per feature and condition
#'
#' Displays the distribution of non-NA measurement counts across features and
#' conditions as a grouped histogram. The \code{max_measured_val} and
#' \code{min_measured_val} summary columns are excluded from the plot.
#'
#' @param input_measured_val_count A wide-format data frame as returned by
#'   \code{\link{missing_values_count}}, containing a \code{feature} column
#'   and one column per condition with the corresponding measured value counts.
#' @param title Character. Plot title. Default: \code{""}.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{missing_values_count}}
#' @importFrom dplyr select
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_histogram scale_fill_viridis_d labs theme_minimal
#' @export
plot_measured_values <- function(input_measured_val_count, title = "") {
  input_measured_val_count %>% 
    select(-max_measured_val, -min_measured_val) %>%
    pivot_longer(cols = names(.)[names(.) != "feature"]) %>% 
    ggplot(aes(x = value, fill = name)) +
    geom_histogram(position= "dodge", alpha=0.8, bins = 30, col = "black") +
    scale_fill_viridis_d() +
    labs(x = "Number of measured values per feature and condition",
         title = title) +
    theme_minimal()
}

#' Histogram of per-feature coefficients of variation per condition
#'
#' Displays the distribution of CV values across features and conditions as a
#' grouped histogram.
#'
#' @param input_cvs A wide-format data frame as returned by
#'   \code{\link{feature_cvs}}, containing a \code{feature} column and one
#'   column per condition with the corresponding CV values.
#' @param title Character. Plot title. Default: \code{""}.
#'
#' @return A \code{ggplot} object.
#'
#' @seealso \code{\link{feature_cvs}}
#' @importFrom tidyr pivot_longer
#' @importFrom ggplot2 ggplot aes geom_histogram scale_fill_viridis_d labs theme_minimal
#' @export
plot_cvs <- function(input_cvs, title = "") {
  gtab <- input_cvs %>% 
    pivot_longer(cols = names(.)[names(.) != "feature"])
  ggplot(data = gtab, 
         aes(x = value, fill = name)) +
    geom_histogram(position= "dodge", alpha=0.8, bins = 30, col = "black") +
    scale_fill_viridis_d() +
    labs(x = "CVs per feature and condition",
         title = title) +
    theme_minimal()
}
