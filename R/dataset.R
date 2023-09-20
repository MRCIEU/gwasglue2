# DataSet class (WIP)

#' An S4 class to represent the Data Set
############################
#'
#' @slot summary_sets A list of SummarySet objects (default NA).
#' @slot overlap_variants among all SummarySets
#' @slot is_resized logical (default FALSE).
#' @slot is_harmonised logical (default FALSE).
#' @slot overall_dropped_SNPs A vector of RSIDs that were removed from the summary_sets.
#' @slot dropped_SNPs A list of pairwise harmonising output (SNPs removed from the summary_sets )
#' @slot palindromic_SNPs A list of pairwise harmonising output.
#' @slot ambiguous_SNPs A list of pairwise harmonising output.
#' @slot incompatible_alleles_SNPs A list of pairwise harmonising output.
#' @slot ld_matrix LD matrix from reference population
#' @slot is_harmonisedLD logical (default FALSE).
#' @slot zscores vector of calculated z-scores
#' @slot susie_marginalised logical (default FALSE).
#' @slot susieR susieR::susie_rss() output
#' @slot is_converted logical (default FALSE).
#' @slot describe A description of the DataSet (default NA).
#' @slot trait_organisation A list with the trait organization within the DataSet (default NA).
#' @export 
#' @rdname DataSet
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
    is_converted = "logical",
    describe = "list",
    trait_organisation = "list"
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
    is_converted = FALSE,
    describe = list(NA_character_),
    trait_organisation = list(NA_character_)
  ),
  contains = c(class(dplyr::tibble()))
)


#' DataSet function
#'
#' @param ... Array of gwasglue2 SummarySet object names.
#' @importFrom methods new
#' @return  A gwasglue2 DataSet object
#' @export 
#' @rdname DataSet
DataSet <- function(...) {
  new("DataSet", summary_sets = as.list(...))
}


#' Get Method to retrieve the GWAS Summary Statistics
#'
#' @param dataset A gwasglue2 DataSet object
#' @param index Index of gwasglue2 SummarySet objects within DataSet
#' @return A tibble with GWAS summary statistics
#' @seealso Similar to [getSummaryData()]
#' @export
#' @docType methods
#' @rdname getData-methods
setGeneric("getData", function(dataset,index) standardGeneric("getData"))
#' @rdname getData-methods
setMethod("getData", "DataSet",
          function(dataset,index) {
            return(dataset@summary_sets[[index]]@ss)
          })



#' Get Method to retrieve the gwasglue2 SummarySet object
#'
#' @param dataset A gwasglue2 DataSet object
#' @param index Index of gwasglue2 SummarySet objects within DataSet
#' @return summarySet gwasglue2 SummarySet object
#' @export
#' @docType methods
#' @rdname getSummarySet-methods
setGeneric("getSummarySet", function(dataset,index) standardGeneric("getSummarySet"))

#' @rdname getSummarySet-methods
setMethod("getSummarySet", "DataSet",
          function(dataset,index) {
            return(dataset@summary_sets[[index]])
          })


#' Size of the DataSet 
#'
#' @param dataset A gwasglue2 DataSet objec
#' @return Number of gwasglue2 SummarySet objects within the DataSet
#' @export
#' @docType methods
#' @rdname getLength-methods
setGeneric("getLength", function(dataset) standardGeneric("getLength"))
#' @rdname getLength-methods
setMethod("getLength", "DataSet",
          function(dataset) {
            return(length(dataset@summary_sets))
          })

#' Calculating Z-scores
#'
#' @param dataset A gwasglue2 DataSet object
#' @return An extra '"zscores"' column in the  GWAS summary statistics tibble. 
#' @export
#' @docType methods
#' @rdname setZscores-methods
setGeneric("setZscores", function(dataset) standardGeneric("setZscores"))
#' @rdname setZscores-methods
setMethod( "setZscores", "DataSet",function(dataset) {
  message("Calculating zscores")
  for (i in seq_along(dataset@summary_sets)){
  dataset@zscores[[i]] <- dataset@summary_sets[[i]]@ss$beta/dataset@summary_sets[[i]]@ss$se
  }
 return(dataset)
    }
)

# #' Get Method to retrieve  Z-scores stored in the SummarySet
# #'
# #' @param summary_set A gwasglue2 SummarySet object
# #' @return The z-scores
# #' @export
# #' @docType methods
# #' @rdname  getZ-scores-methods
# setGeneric("getZscores",function(dataset,index) standardGeneric("getZscores"))
# #' @rdname  getZ-scores-methods
# setMethod("getZscores", "DataSet",
#           function(dataset,index) {
#             return(dataset@zscores[[index]])
#           })


#'  Set the trait organisation of the gwasglue2 DataSet
#' @param dataset A gwasglue2 DataSet object
#' @param ... The organisation of the DataSet
#' @return The gwasglue2 object with the trait organisation stored
#' @export
#' @docType methods
#' @rdname setTraitOrg-methods
setGeneric("setTraitOrg",function(dataset, ...) standardGeneric("setTraitOrg"))
#' @rdname setTraitOrg-methods
setMethod("setTraitOrg", "DataSet", function(dataset,...) {
  # TODO: check if the trait organisation is valid
  dataset@trait_organisation <- list(...)
  return(dataset)
})

#' Get the trait organisation of the gwasglue2 DataSet
#' @param dataset A gwasglue2 DataSet object.
#' @return The trait organisation of the gwasglue2 object
#' @export
#' @docType methods
#' @rdname getTraitOrg-methods
setGeneric("getTraitOrg",function(dataset) standardGeneric("getTraitOrg"))
#' @rdname getTraitOrg-methods
setMethod("getTraitOrg", "DataSet", function(dataset) {
  return(dataset@trait_organisation)
})

#' Assert if the shapes of SummarySets in the gwasglue2 DataSet are the same
#' @param dataset A gwasglue2 DataSet object.
#' @return logical TRUE/FALSE
#' @export 
#' @docType methods
#' @rdname assertSameShape-methods
setGeneric("assertSameShape",function(dataset) standardGeneric("assertSameShape"))
#' @rdname assertSameShape-methods
setMethod("assertSameShape","DataSet", function(dataset){

  # getshapes
  n_sumsets <- getLength(dataset)
  shapes <- lapply(1:n_sumsets, function(i){
    # read summaryset
    shape <- getSummarySet(dataset, i) %>%
            getShape()
  }) %>% unlist()

  # no shape defined within the DataSet
  if(all(is.na(shapes))){
   stop(" No SummarySets inside the DataSet have shapes defined.")
  }

  # no NAs
  if(!any(is.na(shapes))){
    # one shape and no NAS, return TRUE
    if(length(unique(shapes)) == 1){
      return(TRUE)
    }
    # more than one shape and no NAS, return FALSE
    if(length(unique(shapes)) > 1){
      return(FALSE)
    }
  }

  # one shape and NAS, fill NAS with present shape
  if(unique(shapes) %>% is.na() %>% length() == 2){
    stop(" Only one SummarySet has the shape defined in the DataSet.")
    
    # TODO: to be implemented later, maybe elsewhere
    # # the only shape defined
    # shape <- unique(shapes)
    # shape <- shape[!is.na(shape)]

    # for(i in 1:n_sumsets){
    #   # read summaryset
    #   summaryset <- getSummarySet(dataset, i)
    #   # set shape
    #   summaryset <- setShape(summaryset, shape)
    #   # write back
    #   dataset@summary_sets[[i]] <- summaryset
    # }    
    # message(" Only one SummarySet has the shape defined in the DataSet. gwasglue2 automatically assigned all the SummarySets with the '",shape,"'' shape.")
    # return(dataset)    
  }

  # more than one shape and NAS
  if(unique(shapes) %>% is.na() %>% length() > 2){
    return(FALSE)
  }


})


# show method 
setMethod(f = "show", signature="DataSet", definition = function(object) {
  # set description of DataSet
  length <- getLength(object)
  overlap <- object@describe$overlap_variants
  overlap_2_3 <- object@describe$variants_after_harmonization_action2_3
  action <- object@describe$action
  refpop_variants <- object@describe$refpop_variants_avail
  variants_afterLD <- object@describe$variants_after_LDharmonization
  is_harm <- isHarmonised(object)
  is_LDharm <- isHarmonisedLD(object)
  trait_organisation <- getTraitOrg(object)
  
  # write
  cat("A DataSet with", length, "SummarySets.\n")
  
  # trait organisation
  if (is.na(trait_organisation)){
  cat("\nTrait organisation: No trait organisation defined. Use the setTraitOrg() function to add it to the DataSet.\n")
 cat("NOTE: This feature is not fully implemented yet. Analyses can continue without defining it.\n")
  } else{
      cat("\nTrait organisation:", trait_organisation, "\n")
  }
  # harmonisation
  if(isTRUE(is_harm)){
    cat("\nHarmonisation:\n")
    if(action == 1){
    cat("All SummarySets are assumed to be on the forward strand.\n")
    cat(overlap, "variants remaining.\n")
    }
    if(action == 2 || action == 3){
    cat("The SummarySets are not assumed to be all on the forward strand and corrections were made to try to harmonise the data\n.")
    cat(overlap_2_3, "variants remaining.\n")
    }
  } else{
    cat("\nThe DataSet is not harmonised\n")
  }
  # LD harmonisation
  if(isTRUE(is_LDharm)){
    cat("\nHarmonisation against a reference population:\n")
    cat(refpop_variants, "available variants in the reference population data.\n")
    cat(variants_afterLD, "variants remaining.\n")
  } else{
    cat("\nThe DataSet is not harmonised against a reference population.\n")
  }

  cat("\nTo access the GWAS summary data in each of the SummarySets use getData(dataset, index).\n")

})