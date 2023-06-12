# DataSet class (WIP)

#' An S4 class to represent the Data Set
############################
#'
#' @slot summary_sets A list of SummarySet objects (default NA).
#' @slot overlap_variants among all SummarySets
#' @slot is_resized logical (default FALSE).
#' @slot is_harmonised logical (default FALSE).
#' @slot overall_dropped_SNPs A vector of RSIDs that were removed from the summary_sets.
#' @slot dropped_SNPs A list of pairwise harmonising ouptput (SNPs removed from the summary_sets )
#' @slot palindromic_SNPs A list of pairwise harmonising ouptput.
#' @slot ambiguous_SNPs A list of pairwise harmonising ouptput.
#' @slot incompatible_alleles_SNPs A list of pairwise harmonising ouptput.
#' @slot ld_matrix LD matrix buoild from reference population
#' @slot is_harmonisedLD logical (default FALSE).
#' @slot zscores vector of calculated z-scores
#' @slot susie_marginalised logical (default FALSE).
#' @slot susieR susieR::susie_rss() output
#' @slot is_converted logical (default FALSE).
#' @export 
setClass("DataSet",
  slots = c(
    summary_sets = "list",
    overlap_variants = "character",
    is_resized = "logical",
    is_harmonised = "logical",
    overall_dropped_SNPs = "character",
    dropped_SNPs = "list",
    palindromic_SNPs = "list",
    ambiguous_SNPs = "list",
    incompatible_alleles_SNPs = "list",
    ld_matrix = "matrix",
    is_harmonisedLD = "logical",
    zscores = "list",
    susie_marginalised = "logical",
    susieR = "list",
    is_converted = "logical"
  ),
  prototype = prototype(
    sumset = list(NA_character_),
    overlap_variants = NA_character_,
    is_resized = FALSE,
    is_harmonised = FALSE,
    overall_dropped_SNPs = NA_character_,
    dropped_SNPs = list(NA_character_),
    palindromic_SNPs = list(NA_character_),
    ambiguous_SNPs = list(NA_character_),
    incompatible_alleles_SNPs = list(NA_character_),
    ld_matrix = matrix(NA_real_),
    is_harmonisedLD = FALSE,
    zscores = list(NA_real_),
    susie_marginalised = FALSE,
    susieR = list(NA_character_),
    is_converted = FALSE
  ),
  contains = c(class(dplyr::tibble()))
)


#' DataSet function
#'
#' @param ... Array of gwasglue2 SummarySet object names.
#' @importFrom methods new
#' @return  A gwasglue2 DataSet object
#' @export 
DataSet <- function(...) {
  new("DataSet", summary_sets = list(...))
}


#' Get Method to retrieve the GWAS Summary Statistics
#'
#' @param object A gwasglue2 DataSet object
#' @param index Index of gwasglue2 SummarySet objects within DataSet
#' @return A tibble with GWAS summary statistics
#' @seealso Similar to [getSummaryData()]
#' @export
#' @docType methods
#' @rdname getData-methods
setGeneric("getData", function(object,index) standardGeneric("getData"))
#' @rdname getData-methods
setMethod("getData", "DataSet",
          function(object,index) {
            return(object@summary_sets[[index]]@ss)
          })



#' Get Method to retrieve the gwasglue2 SummarySet object
#'
#' @param object A gwasglue2 DataSet objec
#' @param index Index of gwasglue2 SummarySet objects within DataSet
#' @return summarySet gwasglue2 SummarySet object
#' @export
#' @docType methods
#' @rdname getSummarySet-methods
setGeneric("getSummarySet", function(object,index) standardGeneric("getSummarySet"))

#' @rdname getSummarySet-methods
setMethod("getSummarySet", "DataSet",
          function(object,index) {
            return(object@summary_sets[[index]])
          })


#' Size of the DataSet 
#'
#' @param object A gwasglue2 DataSet objec
#' @return Number of gwasglue2 SummarySet objects within the DataSet
#' @export
#' @docType methods
#' @rdname getLength-methods
setGeneric("getLength", function(object) standardGeneric("getLength"))
#' @rdname getLength-methods
setMethod("getLength", "DataSet",
          function(object) {
            return(length(object@summary_sets))
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

