# WIP


  # create a new slot for dataset-specific SNPs
  # get intersect of SNPs across all datasets
  # only keep intersect





# Notes
# - If you have multiple datasets, it could be that 1v2 has a different set in common than 1v3 and 2v3. So we're potentially dropping more SNPs than we need to
# - A single dataset containing 1-3 the idea is analysis requires overlap
# - If you wanna do 2v3 with maximum SNPs then create a new dataset of 2v3 only, excluding 1, to avoid unnecessarily dropping SNPs for 2v3
# - At some point in the future we could implement summary imputation on the fly to maximise SNP inclusion



# pipeline
# - 1. Extract data for requested IDs at requested SNPs
# - 2. Identify overlapping SNPs
# - 3. Harmonise
# - 4. Identify any new SNPs that need to be removed due to harmonisation issues (strand issues, allele issues), and remove them from all summary sets



setGeneric("overlapSNP",function(object) standardGeneric("overlapSNP"))
setMethod("overlapSNP", "DataSet", function(object) {
            # Find overlapped variants among all summarySets
            variants_list <- sapply(1:length(object@summary_sets), function(i) list(object@summary_sets[[i]]@variants))
            overlap <- Reduce(intersect, variants_list)

            object@overlap_SNPs <- overlap
            message("\nThere are ", length(object@overlap_SNPs), " SNPs in common among all SummarySets")

            #TODO if length(object@overlap_SNPs) == dim(object@summary_sets[[i]]@ss)[1] message(no resizing)

              message("We are  now resizing the SumarySets")
              for (i in seq_along(object@summary_sets)) {
                object@summary_sets[[i]]@ss <- object@summary_sets[[i]]@ss[which(object@summary_sets[[i]]@ss$rsid %in% object@overlap_SNPs), ]
              }
              message("Done!\n")
              object@is_resized <- TRUE

            return(object)
          })


# Set and get methods for harmonise data

#' Harmonise the alleles and effects between two summary sets
#' @param object The DataSet object
#' @param action Level of strictness in dealing with SNPs.
#' * `action = 1`: Assume all alleles are coded on the forward strand, i.e. do not attempt to flip alleles
#' * `action = 2`: Try to infer positive strand alleles, using allele frequencies for palindromes (default, conservative);
#' * `action = 3`: Correct strand for non-palindromic SNPs, and drop all palindromic SNPs from the analysis (more conservative).
#' @param tolerance Tolerance value (default 0.08).
setGeneric("harmoniseData", function(object, tolerance, action) standardGeneric("harmoniseData"))
setMethod( "harmoniseData", "DataSet", function(object,tolerance = 0.08,action = 2){
  dat1 <- object@summary_sets[[1]]@ss
  message("Gwasglue is now harmonising! ")
  for (i in seq_along(object@summary_sets)[-1]){
    count <- i - 1
    dat2 <-object@summary_sets[[i]]@ss
    message("\nHarmonising ", dat1$id[1], " and ", dat2$id[1],"\n")

    rsid <- dat1$rsid
    A1 <- dat1$ea
    A2 <- dat1$nea
    B1 <- dat2$ea
    B2 <- dat2$nea
    betaA <- dat1$beta
    betaB <- dat2$beta
    fA <- dat1$eaf
    fB <- dat2$eaf

    h <- harmonise(rsid, A1, A2, B1, B2, betaA, betaB, fA, fB, tolerance=tolerance, action=action)
    # remove extra columns
    h1 <-h[,c("rsid","A1","A2","betaA","fA")]
    h1 <- dplyr::rename(h1, c(ea = A1, nea = A2, beta = betaA, eaf = fA))
    h2 <-h[,c("rsid","B1","B2","betaB","fB")]
    h2 <- dplyr::rename(h2, c(ea = B1, nea = B2, beta = betaB, eaf = fB))
    # replace c(ea, nea, beta, eaf) columns
    dat1 <- merge(subset(dat1, select=-c(ea, nea, beta, eaf)), h1, by="rsid")
    dat2 <- merge(subset(dat2, select=-c(ea, nea, beta, eaf)), h2, by="rsid")

    object@summary_sets[[1]]@ss <- dplyr::as_tibble(dat1) #TODO check if there is any condition where dat1 changes?
    object@summary_sets[[i]]@ss <- dplyr::as_tibble(dat2)

    object@dropped_SNPs[[count]] <-  h$rsid[h$keep == FALSE]
    names(object@dropped_SNPs)[[count]] <- paste0(dat1$id[1],"_vs_",dat2$id[1])
    object@palindromic_SNPs[[count]] <-  h$rsid[h$palindromic == TRUE]
    names(object@palindromic_SNPs)[[count]] <- paste0(dat1$id[1],"_vs_",dat2$id[1])
    object@ambiguous_SNPs[[count]] <-  h$rsid[h$ambiguous == TRUE]
    names(object@ambiguous_SNPs)[[count]] <- paste0(dat1$id[1],"_vs_",dat2$id[1])
    object@incompatible_alleles_SNPs[[count]] <-  h$rsid[h$remove == TRUE]
    names(object@incompatible_alleles_SNPs)[[count]] <- paste0(dat1$id[1],"_vs_",dat2$id[1])
     }
  message("\nDone harmonising!")

  # TODO harmonise_make_snp_effects_positive?

  # at the end of harmonisation remove all union(object@dropped_SNPs)
  dropped_SNPs <- sapply(1:length(object@dropped_SNPs), function(i) c(object@dropped_SNPs[[i]]))
  dropped_SNPs <- Reduce(union, dropped_SNPs)

  object@overall_dropped_SNPs <- dropped_SNPs
  message("\nThere are ", length(dropped_SNPs), " SNPs to remove among all SumarySets after harmonising")

  #TODO if length(object@overall_dropped_SNPs) == dim(object@summary_sets[[i]]@ss)[1] message(no resizing)
  message("We are  now resizing the SumarySets")
  "%ni%" <- Negate("%in%")
  for (i in seq_along(object@summary_sets)) {
    object@summary_sets[[i]]@ss <- object@summary_sets[[i]]@ss[which(object@summary_sets[[i]]@ss$rsid %ni% dropped_SNPs), ]
  }
  message("Done!\n")

    object@is_harmonised <- TRUE
    return(object)
}
)


setGeneric("isHarmonised",function(object) standardGeneric("isHarmonised"))
setMethod("isHarmonised","DataSet",function(object) {
  return(object@is_harmonised)
}
)
