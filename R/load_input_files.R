#' Load a Spectronaut Protein Group MLP report
#'
#' Reads a tab-separated Spectronaut Protein Group report, verifies that all
#' required columns are present, selects and renames relevant columns, and
#' adds a log2-transformed quantity column.
#'
#' @param path_to_input_protein Character. Path to the Spectronaut \code{.tsv}
#'   Protein Group (MLP) report file.
#'
#' @return A \code{tibble} in long format with one row per protein group-run
#'   combination and the following columns:
#'   \describe{
#'     \item{condition}{Experimental condition (\code{R.Condition}).}
#'     \item{filename}{Raw file name (\code{R.FileName}).}
#'     \item{replicate}{Replicate label (\code{R.Replicate}).}
#'     \item{pg_informations}{FASTA headers (\code{PG.FastaHeaders}).}
#'     \item{protein_group_accessions}{Protein accessions
#'       (\code{PG.ProteinAccessions}).}
#'     \item{nr_precursors_used_for_quantification}{Number of precursors used
#'       for quantification (\code{PG.NrOfPrecursorsUsedForQuantification}).}
#'     \item{quantity}{Raw protein group quantity (\code{PG.Quantity}).}
#'     \item{log2_quan}{Log2-transformed quantity.}
#'   }
#'
#' @importFrom data.table fread
#' @importFrom dplyr select rename mutate distinct all_of
#' @export
load_Spectronaut_protein_input <- function(path_to_input_protein) {
  
  # Load the proteome data from Spectronaut report
  proteome_data <- fread(path_to_input_protein)
  
  # Check if the necessary columns are present
  required_columns <- c("R.Condition",
                        "R.FileName",
                        "R.Replicate",
                        "PG.FastaHeaders",
                        "PG.ProteinAccessions",
                        "PG.NrOfPrecursorsUsedForQuantification",
                        "PG.Quantity")
  if (!all(required_columns %in% colnames(proteome_data))) {
    stop("The input file does not contain all required columns: ", paste(required_columns, collapse = ", "))
  }
  
  proteome_data$contaminant <- ifelse(grepl(CONTAMINANT_FLAG, unlist(proteome_data[,names(proteome_data) == CONTAMINANT_COLUMN, with = F]), fixed = T), TRUE, FALSE)
  
  proteome_data %>% 
    select(all_of(c(required_columns, "contaminant"))) %>% 
    rename(condition = R.Condition,
           filename = R.FileName,
           replicate = R.Replicate,
           pg_informations = PG.FastaHeaders,
           protein_group_accessions = PG.ProteinAccessions,
           nr_precursors_used_for_quantification = PG.NrOfPrecursorsUsedForQuantification,
           quantity = PG.Quantity) %>% 
    mutate(log2_quan = log2(quantity)) %>%
    distinct() 
  
}

#' Load a Spectronaut PTM Site BGS report
#'
#' Reads a tab-separated Spectronaut PTM Site BGS report, retains
#' phosphorylation events only (\code{Phospho (STY)}), collapses multiple
#' \code{PTM.CollapseKey} values into a single semicolon-separated
#' \code{site_ID} per PTM group, recodes zero quantities as \code{NA}, and
#' adds a log2-transformed quantity column.
#'
#' @param path_to_input_ptm Character. Path to the Spectronaut \code{.tsv}
#'   PTM Site BGS report file.
#'
#' @return A \code{tibble} in long format with one row per phosphosite
#'   group-run combination and the following columns:
#'   \describe{
#'     \item{condition}{Experimental condition (\code{R.Condition}).}
#'     \item{filename}{Raw file name (\code{R.FileName}).}
#'     \item{replicate}{Replicate label (\code{R.Replicate}).}
#'     \item{peptidoform}{PTM peptidoform group (\code{PTM.Group}).}
#'     \item{PTM_multiplicity}{PTM multiplicity (\code{PTM.Multiplicity}).}
#'     \item{quantity}{Raw PTM site quantity (\code{PTM.Quantity}); zeros
#'       are replaced with \code{NA}.}
#'     \item{site_ID}{Semicolon-separated collapse keys for all PTM sites in
#'       the group (\code{PTM.CollapseKey}).}
#'     \item{log2_quan}{Log2-transformed quantity.}
#'   }
#'
#' @importFrom data.table fread
#' @importFrom dplyr arrange select filter group_by mutate ungroup rename distinct all_of
#' @export
load_Spectronaut_ptm_input <- function(path_to_input_ptm) {
  
  # Load the PTM data from Spectronaut report
  ptm_data <- fread(path_to_input_ptm)
  
  # Check if the necessary columns are present
  required_columns <- c("R.Condition",
                        "R.FileName",
                        "R.Replicate",
                        "PTM.CollapseKey",
                        "PTM.Group",
                        "PTM.Multiplicity",
                        "PG.FastaHeaders",
                        "PTM.Quantity",
                        "PTM.ModificationTitle"
                        )
  if (!all(required_columns %in% colnames(ptm_data))) {
    stop("The input file does not contain all required columns: ", paste(required_columns, collapse = ", "))
  }
  
  cat("Keep only phorphorylations and combine the PTM.CollapseKey per PTM.Group...\n")
  
  ptm_data$contaminant <- ifelse(grepl(CONTAMINANT_FLAG, unlist(ptm_data[,names(ptm_data) == CONTAMINANT_COLUMN, with = F]), fixed = T), TRUE, FALSE)
  
  cat("Remove \\_number at the end of the site ID...\n")
  ptm_data$PTM.CollapseKey <- gsub("\\_\\d$", "", ptm_data$PTM.CollapseKey)
  
  ptm_data %>% 
    arrange(PTM.ProteinId, PTM.SiteLocation) %>%
    select(all_of(c(required_columns, "contaminant"))) %>% 
    filter(PTM.ModificationTitle == "Phospho (STY)") %>% 
    group_by(PTM.Group, R.FileName, R.Condition, R.Replicate, PTM.Multiplicity, PTM.Quantity, contaminant) %>%
    mutate(site_ID = toString(PTM.CollapseKey) %>% 
             gsub(", ", ";", .)) %>%
    ungroup()  %>% 
    select(!c(PTM.ModificationTitle, PTM.CollapseKey)) %>%
    rename(condition = R.Condition,
           filename = R.FileName,
           replicate = R.Replicate,
           pg_informations = PG.FastaHeaders,
           peptidoform = PTM.Group,
           PTM_multiplicity = PTM.Multiplicity,
           quantity = PTM.Quantity
           )  %>%
    distinct() -> ptm_data_filtered
  
  cat("Replace 0s with NAs...\n")
  
  ptm_data_filtered$quantity[ptm_data_filtered$quantity == 0] <- NA
  
  ptm_data_filtered %>% 
    mutate(log2_quan = log2(quantity)) 
    
}

#' Load and combine proteome and phosphoproteome Spectronaut reports
#'
#' Convenience wrapper that calls \code{\link{load_Spectronaut_protein_input}}
#' and \code{\link{load_Spectronaut_ptm_input}} and returns both tables in a
#' named list.
#'
#' @param path_to_input_protein Character. Path to the Spectronaut Protein
#'   Group (MLP) \code{.tsv} report.
#' @param path_to_input_ptm Character. Path to the Spectronaut PTM Site BGS
#'   \code{.tsv} report.
#'
#' @return A named list with two elements:
#'   \describe{
#'     \item{proteome}{Table returned by
#'       \code{\link{load_Spectronaut_protein_input}}.}
#'     \item{ptm}{Table returned by
#'       \code{\link{load_Spectronaut_ptm_input}}.}
#'   }
#'
#' @seealso \code{\link{load_Spectronaut_protein_input}},
#'   \code{\link{load_Spectronaut_ptm_input}}
#' @export
open_phos_prot_files <- function(path_to_input_protein, path_to_input_ptm) {
  list(
    proteome = load_Spectronaut_protein_input(path_to_input_protein),
    ptm = load_Spectronaut_ptm_input(path_to_input_ptm)
  )
}
