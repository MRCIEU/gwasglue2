

# Set and get methods for convert methods

#' Title
#'
#' @param object
#'
#' @return
#' @export
#'
#' @examples
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
  for (i in seq_along(object@sumset)){
    # change column names
    object@sumset[[i]]@ss <- object@sumset[[i]]@ss %>%
      rename(c(SNP = rsid,effect_allele=ea, other_allele = nea, pval=p, samplesize=n))
    #  add extra column named "mr_keep"
    object@sumset[[i]]@ss$mr_keep <- !is.na(object@sumset[[i]]@ss$beta) &
      !is.na(object@sumset[[i]]@ss$se) &
      !is.na(object@sumset[[i]]@ss$effect_allele) &
      !is.na(object@sumset[[i]]@ss$SNP)
    # Exposure
    if (object@sumset[[i]]@mr_label == "exposure"){
      c <-which(colnames(object@sumset[[i]]@ss)!="SNP")
      colnames(object@sumset[[i]]@ss)[c] <- paste0(colnames(object@sumset[[i]]@ss)[c], ".exposure")
      #  add extra column named "exposure"
      object@sumset[[i]]@ss$exposure <- "exposure" #TODO change how to fill this column?
    }
    # Outcome
    if (object@sumset[[i]]@mr_label == "outcome"){
      c <-which(colnames(object@sumset[[i]]@ss)!="SNP")
      colnames(object@sumset[[i]]@ss)[c] <- paste0(colnames(object@sumset[[i]]@ss)[c], ".outcome")
      #  add extra column named "outcome"
      object@sumset[[i]]@ss$outcome <- "outcome"
    }

  }
  object@is_converted <- TRUE
  return(object)
}
)

