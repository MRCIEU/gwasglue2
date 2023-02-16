############################
# DataSet class (WIP)

#' An S4 class to represent the Data Set
#'
#' @slot summary_sets A list of SummarySet objects (default NA).
#' @slot overlap_SNPs among all SummarySets
#' @slot is_resized logical (default FALSE).
#' @slot is_harmonised logical (default FALSE).
#' @slot overall_dropped_SNPs A vector of RSIDs that were removed from the summary_sets.
#' @slot dropped_SNPs A list of pairwise harmonising ouptput (SNPs removed from the summary_sets )
#' @slot palindromic_SNPs A list of pairwise harmonising ouptput.
#' @slot ambiguous_SNPs A list of pairwise harmonising ouptput.
#' @slot incompatible_alleles_SNPs A list of pairwise harmonising ouptput.
#' @slot ld_matrices A list of LD correlation matrices (default NA).
#' @slot is_harmonisedLD logical (default FALSE).
#' @slot is_converted logical (default FALSE).
#' @export 
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
    zscores = "list",
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
    zscores = list(NA_real_),
    is_converted = FALSE
  ),
  contains = c(class(dplyr::tibble()))
)


#' DataSet function
#'
#' @param ... Array of SummarySet object names.
#' @importFrom methods new
#' @return  A DataSet S4 object
DataSet <- function(...) {
  new("DataSet", summary_sets = list(...))
}

# Get Methods for summary set (similar in Summaryset class)
setGeneric("getData", function(object,...) standardGeneric("getData"))
setMethod("getData", "DataSet",
          function(object,index) {
            return(object@summary_sets[[index]]@ss)
          })


# Get Methods for length of DAtaset
setGeneric("getLength", function(object) standardGeneric("getLength"))
setMethod("getLength", "DataSet",
          function(object) {
            return(object@summary_sets[[index]]@ss)
          })

# Set and get methods for zscores
setGeneric("setZscores", function(object) standardGeneric("setZscores"))
setMethod( "setZscores", "DataSet",function(object) {
  message("Calculating zscores")
  for (i in seq_along(object@summary_sets)){
  object@zscores[[i]] <- object@summary_sets[[i]]@ss$beta/object@summary_sets[[i]]@ss$se
  }
 return(object)
    }
)

setGeneric("getZscores",function(object,...) standardGeneric("getZscores"))
setMethod("getZscores", "DataSet",
          function(object,index) {
            return(object@zscores[[index]])
          })

