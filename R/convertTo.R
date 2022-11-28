#' @include dataset.R
NULL

# Set and get methods for convert methods

#IsConverted
setGeneric("isConverted",function(object) standardGeneric("isConverted"))
setMethod("isConverted","DataSet",
          function(object) {
            return(object@is_converted)
          }
)

# Convert for TwoSampleMR
setGeneric("convertForTwoSampleMR",function(object) standardGeneric("convertForTwoSampleMR"))
setMethod("convertForTwoSampleMR", "DataSet", function(object) {
  # TODO add warnings for when the data  is already converted
  message("Gwasglue is now converting the data to TwoSampleMR format!")
  for (i in seq_along(object@summary_sets)){
    # change column names
    object@summary_sets[[i]]@ss <- dplyr::rename(object@summary_sets[[i]]@ss, c(SNP = rsid,effect_allele=ea, other_allele = nea, pval=p, samplesize=n)) # nolint
    #  add extra column named "mr_keep"
    object@summary_sets[[i]]@ss$mr_keep <- !is.na(object@summary_sets[[i]]@ss$beta) &
      !is.na(object@summary_sets[[i]]@ss$se) &
      !is.na(object@summary_sets[[i]]@ss$effect_allele) &
      !is.na(object@summary_sets[[i]]@ss$SNP)
    # Exposure
    if (object@summary_sets[[i]]@mr_label == "exposure"){
      c <-which(colnames(object@summary_sets[[i]]@ss) != c("SNP","mr_keep"))
      colnames(object@summary_sets[[i]]@ss)[c] <- paste0(colnames(object@summary_sets[[i]]@ss)[c], ".exposure")
      #  add extra column named "exposure"
      object@summary_sets[[i]]@ss$exposure <- "exposure" #TODO change how to fill this column?
    }
    # Outcome
    if (object@summary_sets[[i]]@mr_label == "outcome"){
      c <-which(colnames(object@summary_sets[[i]]@ss) != c("SNP","mr_keep"))
      colnames(object@summary_sets[[i]]@ss)[c] <- paste0(colnames(object@summary_sets[[i]]@ss)[c], ".outcome")
      #  add extra column named "outcome"
      object@summary_sets[[i]]@ss$outcome <- "outcome"
    }

  }
  object@is_converted <- TRUE
  return(object)
}
)



