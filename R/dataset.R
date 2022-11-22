############################
# DataSet class (WIP)

#' An S4 class to represent the Data Set
#'
#' @slot summary_sets A list of SummarySet objects (default NA).
#' @slot overlap_SNPs among all SummarySets
#' @slot is_resized logical (default FALSE).
#' @slot is_harmonised logical (default FALSE).
#' @slot overall_dropped_SNPs
#' @slot dropped_SNPs
#' @slot palindromic_SNPs
#' @slot ambiguous_SNPs
#' @slot incompatible_alleles_SNPs
#' @slot ld_matrices A list of LD correlation matrices (default NA).
#' @slot is_harmonisedLD logical (default FALSE).
#' @slot is_converted logical (default FALSE).
setClass("DataSet",
  slots = c(
    summary_sets = "list",
    overlap_SNPs = "character",
    is_resized = "logical",
    is_harmonised = "logical",
    overall_dropped_SNPs = "character",
    dropped_SNPs = "list",
    palindromic_SNPs = "list",
    ambiguous_SNPs = "list",
    incompatible_alleles_SNPs = "list",
    ld_matrices = "list",
    is_harmonisedLD = "logical",
    is_converted = "logical"
  ),
  prototype = prototype(
    sumset = list(NA_character_),
    overlap_SNPs = NA_character_,
    is_resized = FALSE,
    is_harmonised = FALSE,
    overall_dropped_SNPs = NA_character_,
    dropped_SNPs = list(NA_character_),
    palindromic_SNPs = list(NA_character_),
    ambiguous_SNPs = list(NA_character_),
    incompatible_alleles_SNPs = list(NA_character_),
    ld_matrices = list(NA_character_),
    is_harmonisedLD = FALSE,
    is_converted = FALSE
  ),
  contains = c(class(tibble()))
)


#' DataSet function
#'
#' @param ...
#'
#' @return
DataSet <- function(...) {
  new("DataSet", summary_sets = list(...))
}

