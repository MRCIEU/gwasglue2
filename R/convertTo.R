#' @include dataset.R
NULL

# Set and get methods for convert methods

#IsConverted
setGeneric("isConverted",function(dataset) standardGeneric("isConverted"))
setMethod("isConverted","DataSet",
          function(dataset) {
            return(dataset@is_converted)
          }
)


# Convert for TwoSampleMR
#' Look for overlapped variants between SummarySets in the DataSet and Resize
#' 
#' @param dataset The gwasglue2 DataSet object
#' @return The gwasglue2 SummarySet object converted to TwoSampleMR format
#' @export
#' @docType methods
#' @rdname convertForTwoSampleMR
setGeneric("convertForTwoSampleMR", function(dataset) standardGeneric("convertForTwoSampleMR"))
#' @rdname convertForTwoSampleMR
setMethod("convertForTwoSampleMR", "DataSet", function(dataset) {
  # TODO add warnings for when the data  is already converted
  message("Gwasglue is now converting the data to TwoSampleMR format!")
  for (i in seq_along(dataset@summary_sets)){
    # change column names
    dataset@summary_sets[[i]]@ss <- dplyr::rename(dataset@summary_sets[[i]]@ss, c(SNP = rsid,effect_allele=ea, other_allele = nea, pval=p, samplesize=n)) # nolint
    #  add extra column named "mr_keep"
    dataset@summary_sets[[i]]@ss$mr_keep <- !is.na(dataset@summary_sets[[i]]@ss$beta) &
      !is.na(dataset@summary_sets[[i]]@ss$se) &
      !is.na(dataset@summary_sets[[i]]@ss$effect_allele) &
      !is.na(dataset@summary_sets[[i]]@ss$SNP)
    # Exposure
    if (dataset@summary_sets[[i]]@attributes$mr_label == "exposure"){
      c <- colnames(dataset@summary_sets[[i]]@ss) %ni% c("SNP","mr_keep")
      colnames(dataset@summary_sets[[i]]@ss)[c] <- paste0(colnames(dataset@summary_sets[[i]]@ss)[c], ".exposure")
      #  add extra column named "exposure"
      dataset@summary_sets[[i]]@ss$exposure <- "exposure" #TODO change how to fill this column?
    }
    # Outcome
    if (dataset@summary_sets[[i]]@attributes$mr_label == "outcome"){
      c <- colnames(dataset@summary_sets[[i]]@ss) %ni% c("SNP","mr_keep")
      colnames(dataset@summary_sets[[i]]@ss)[c] <- paste0(colnames(dataset@summary_sets[[i]]@ss)[c], ".outcome")
      #  add extra column named "outcome"
      dataset@summary_sets[[i]]@ss$outcome <- "outcome"
    }

  }
  dataset@is_converted <- TRUE
  return(dataset)
}
)