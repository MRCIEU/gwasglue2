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


#' Convert tool to TwoSampleMR format
#' @description Converts SummarySets within a Dataset to a format that can be read by TwoSampleMR
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


#' Convert tool to TwoSampleMR format
#' @description Converts SummarySets within a Dataset to a format that can be read by TwoSampleMR
#' @param dataset The gwasglue2 DataSet object
#' @return A TwoSampleMR data.frame input format
#' @export
#' @docType methods
#' @rdname convertToTwoSampleMR
setGeneric("convertToTwoSampleMR", function(dataset) standardGeneric("convertToTwoSampleMR"))
#' @rdname convertToTwoSampleMR
setMethod("convertToTwoSampleMR", "DataSet", function(dataset) {
  # TODO add warnings for when the data  is already converted


  message("Gwasglue is now converting the data to TwoSampleMR format!")
  for (i in seq_along(dataset@summary_sets)){
    # Exposure
    if (dataset@summary_sets[[i]]@attributes$mr_label == "exposure"){
      
      if(i==1){
          # change column names
      mr_exposure <- dplyr::rename(dataset@summary_sets[[i]]@ss, c(SNP = rsid,effect_allele=ea, other_allele = nea, pval=p, samplesize=n)) # nolint
      #  add extra column named "mr_keep"
      mr_exposure$mr_keep <- !is.na( mr_exposure$beta) &
        !is.na( mr_exposure$se) &
        !is.na( mr_exposure$effect_allele) &
        !is.na( mr_exposure$SNP)

        c <- colnames(mr_exposure) %ni% c("SNP","mr_keep")
        names(mr_exposure)[c] <- paste0(colnames(mr_exposure[c]), ".exposure")
        #  add extra column named "exposure"
        mr_exposure$exposure <- "exposure" #TODO change how to fill this column?
      } else{
         ss <- dplyr::rename(dataset@summary_sets[[i]]@ss, c(SNP = rsid,effect_allele=ea, other_allele = nea, pval=p, samplesize=n)) # nolint
      #  add extra column named "mr_keep"
        ss$mr_keep <- !is.na(ss$beta) &
        !is.na( ss$se) &
        !is.na(ss$effect_allele) &
        !is.na(ss$SNP)
        c <- colnames(ss) %ni% c("SNP","mr_keep")
        names(ss)[c] <- paste0(colnames(ss[c]), ".exposure")
      
      mr_exposure <- dplyr::bind_rows(mr_exposure,ss)
      }
    }

    # Outcome
    if (dataset@summary_sets[[i]]@attributes$mr_label == "outcome"){
        
      if(i==1){
          # change column names
      mr_outcome <- dplyr::rename(dataset@summary_sets[[i]]@ss, c(SNP = rsid,effect_allele=ea, other_allele = nea, pval=p, samplesize=n)) # nolint
      #  add extra column named "mr_keep"
      mr_outcome$mr_keep <- !is.na( mr_outcome$beta) &
        !is.na( mr_outcome$se) &
        !is.na( mr_outcome$effect_allele) &
        !is.na( mr_outcome$SNP)

        c <- colnames(mr_outcome) %ni% c("SNP","mr_keep")
        names(mr_outcome)[c] <- paste0(colnames(mr_outcome[c]), ".outcome")
        #  add extra column named "exposure"
        mr_outcome$exposure <- "outcome" #TODO change how to fill this column?
      } else{
         ss <- dplyr::rename(dataset@summary_sets[[i]]@ss, c(SNP = rsid,effect_allele=ea, other_allele = nea, pval=p, samplesize=n)) # nolint
      #  add extra column named "mr_keep"
        ss$mr_keep <- !is.na(ss$beta) &
        !is.na( ss$se) &
        !is.na(ss$effect_allele) &
        !is.na(ss$SNP)
        c <- colnames(ss) %ni% c("SNP","mr_keep")
        names(ss)[c] <- paste0(colnames(ss[c]), ".outcome")
      
      mr_outcome <- dplyr::bind_rows(mr_outcome,ss)
      }
    }

  }

  mr_input <- merge(mr_exposure,mr_outcome, by = c("SNP", "mr_keep")) 
  return(mr_input)
}
)
