#  The S4 methods in this file use the functions in `harmonise.R`


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



# Rules

# if not assuming forward strand, don't handle multi-allelics, arrange the dataset alphabetically (rsid, a1, a2), then drop duplicates, then only merge on rsid
# if assuming forward strand, make internal identifier variantid that looks like chr_pos_a1_a2
# a1 and a2 are alphabetical
# if a1 or a2 is longer than 10 characters, create hash
# if chr_pos not available then rsid_a1_a2 ? TODO
# switch effect allele to be the alphabetical allele
# arrange dataset by variant_id

#' Look for overlapped variants between SummarySets in the DataSet and Resize
#' 
#' @param dataset The gwasglue2 DataSet object
#' @param action Level of strictness in dealing with SNPs during harmonisation.T#' * `action = 1`: Assume all alleles are coded on the forward strand, i.e. do not attempt to flip alleles
#' * `action = 2`: Try to infer positive strand alleles, using allele frequencies for palindromes (default, conservative);
#' * `action = 3`: Correct strand for non-palindromic SNPs, and drop all palindromic SNPs from the analysis (more conservative).
#' @return The gwasglue2 DataSet object resized
#' @export 
#' @docType methods
#' @rdname overlapVariants-methods
setGeneric("overlapVariants",function(dataset, action) standardGeneric("overlapVariants"))
#' @rdname overlapVariants-methods
setMethod("overlapVariants", "DataSet", function(dataset, action) {
  # SummarySet is already  standardised in create_summaryset()
   
  # clean describe slot and set action
  dataset@describe <- dataset@describe[-1]
  dataset@describe$action <- action

  if(action == 1){
    # Find overlapped variants among all summarySets using variantid
    variants_list <- sapply(1:length(dataset@summary_sets), function(i) list(dataset@summary_sets[[i]]@ss$variantid))
    overlap <- Reduce(intersect, variants_list)


    #TODO if length(dataset@overlap_variants) == dim(dataset@summary_sets[[i]]@ss)[1] message(no resizing)
    #TODO if chr_pos not available then rsid_a1_a2 ?

    message("We are  now resizing the SummarySets")
    for (i in seq_along(dataset@summary_sets)) {
      dataset@summary_sets[[i]]@ss <- dataset@summary_sets[[i]]@ss[which(dataset@summary_sets[[i]]@ss$variantid %in% overlap), ] %>% dplyr::arrange(., variantid, chr, position, ea, nea)
    }
    
     dataset@overlap_variants <- unique(dataset@summary_sets[[1]]@ss$rsid)
     dataset@is_harmonised <- TRUE
    message("\nThere are ", length(overlap), " variants in common among all SummarySets. Data is harmonised.")
    dataset@describe$overlap_variants <- length(overlap)
   
  }
  
  if(action == 2 || action == 3){
  
    # Find overlapped variants among all summarySets
     #  check for multiallelic variants and drop them
    #  TODO probably will vave to change this to use chr:position instead of rsid
    variants_list <- sapply(1:length(dataset@summary_sets), function(i){
      variants <- list(dataset@summary_sets[[i]]@ss$rsid[which(dataset@summary_sets[[i]]@ss$rsid %ni% names(which(table(dataset@summary_sets[[i]]@ss$rsid) > 1)))])
    })
  
    overlap <- Reduce(intersect, variants_list)
    message("We are  now resizing the SummarySets")
              
    for (i in seq_along(dataset@summary_sets)) {
      dataset@summary_sets[[i]]@ss <- dataset@summary_sets[[i]]@ss[which(dataset@summary_sets[[i]]@ss$rsid %in% overlap), ] %>% dplyr::arrange(., variantid, chr, position, ea, nea)
    }
   

  dataset@overlap_variants <- overlap
    message("\nThere are ", length(dataset@overlap_variants), " variants in common among all SummarySets")
    dataset@describe$overlap_variants <- length(overlap)

  }

  message("Done!\n")
  dataset@is_resized <- TRUE

  return(dataset)
})



# Set and get methods for harmonise data

#' Harmonise the alleles and effects between two summary sets
#' 
#' @param dataset The gwasglue2 DataSet object
#' @param action Level of strictness in dealing with SNPs.
#' * `action = 1`: Assume all alleles are coded on the forward strand, i.e. do not attempt to flip alleles
#' * `action = 2`: Try to infer positive strand alleles, using allele frequencies for palindromes (default, conservative);
#' * `action = 3`: Correct strand for non-palindromic SNPs, and drop all palindromic SNPs from the analysis (more conservative).
#' @param tolerance Tolerance value.
#' @return The gwasglue2 DataSet object harmonised
#' @export 
#' @docType methods
#' @rdname harmoniseData-methods
setGeneric("harmoniseData", function(dataset, tolerance, action) standardGeneric("harmoniseData"))
#' @rdname harmoniseData-methods
setMethod( "harmoniseData", "DataSet", function(dataset, tolerance, action){
  
  

  message("Gwasglue is now harmonising! ")
  for (i in seq_along(dataset@summary_sets)[-1]){
      count <- i - 1
      
      # We are always going to harmonise agaist the firts SummarySet and use variantid instead of rsid
      dat1 <- dataset@summary_sets[[1]]@ss
      rsid <- dat1$variantid
      A1 <- dat1$ea
      A2 <- dat1$nea
      fA <- dat1$eaf
      betaA <- dat1$beta

      dat2 <-dataset@summary_sets[[i]]@ss
      message("\nHarmonising ", dat1$id[1], " and ", dat2$id[1],"\n")

      B1 <- dat2$ea
      B2 <- dat2$nea
      betaB <- dat2$beta
      fB <- dat2$eaf

      h <- harmonise(rsid, A1, A2, B1, B2, betaA, betaB, fA, fB, tolerance=tolerance, action=action)

      # remove extra columns
      h1 <-h[,c("rsid","A1","A2","betaA","fA")]
      h1 <- dplyr::rename(h1, c(variantid = rsid, ea = A1, nea = A2, beta = betaA, eaf = fA))
      h2 <-h[,c("rsid","B1","B2","betaB","fB")]
      h2 <- dplyr::rename(h2, c(variantid = rsid,ea = B1, nea = B2, beta = betaB, eaf = fB))
      # replace c(ea, nea, beta, eaf) columns
      dat1 <- merge(subset(dat1, select=-c(ea, nea, beta, eaf)), h1, by="variantid")
      dat2 <- merge(subset(dat2, select=-c(ea, nea, beta, eaf)), h2, by="variantid")

      dataset@summary_sets[[1]]@ss <- dplyr::as_tibble(dat1) #TODO check if there is any condition where dat1 changes?
      dataset@summary_sets[[i]]@ss <- dplyr::as_tibble(dat2)

      # TODO fill slots
      dataset@dropped_SNPs[[count]] <-  h$variantid[h$keep == FALSE]
      # # names(dataset@dropped_SNPs)[[count]] <- paste0(dat1$id[1],"_vs_",dat2$id[1])
      # dataset@palindromic_SNPs[[count]] <-  h$variantid[h$palindromic == TRUE]
      # # names(dataset@palindromic_SNPs)[[count]] <- paste0(dat1$id[1],"_vs_",dat2$id[1])
      # dataset@ambiguous_SNPs[[count]] <-  h$variantid[h$ambiguous == TRUE]
      # # names(dataset@ambiguous_SNPs)[[count]] <- paste0(dat1$id[1],"_vs_",dat2$id[1])
      # dataset@incompatible_alleles_SNPs[[count]] <-  h$variantid[h$remove == TRUE]
      # # names(dataset@incompatible_alleles_SNPs)[[count]] <- paste0(dat1$id[1],"_vs_",dat2$id[1]) 
  }
  message("\nDone harmonising!")
  
  # TODO harmonise_make_snp_effects_positive?

  # at the end of harmonisation remove all union(dataset@dropped_SNPs)
  if(length(dataset@dropped_SNPs) != 0){

    dropped_SNPs <- sapply(1:length(dataset@dropped_SNPs), function(i) c(dataset@dropped_SNPs[[i]]))
    dropped_SNPs <- Reduce(union, dropped_SNPs)

    dataset@overall_dropped_SNPs <- dropped_SNPs
    message("\nThere are ", length(dropped_SNPs), " variants to remove among all SummarySets after harmonising")

    #TODO if length(dataset@overall_dropped_SNPs) == dim(dataset@summary_sets[[i]]@ss)[1] message(no resizing)
    message("We are  now resizing the SummarySets")
  
    for (i in seq_along(dataset@summary_sets)) {
        dataset@summary_sets[[i]]@ss <- dataset@summary_sets[[i]]@ss[which(dataset@summary_sets[[i]]@ss$variantid %ni% dropped_SNPs), ]
    }

     dataset@describe$variants_after_harmonization_action2_3 <- nrow(dataset@summary_sets[[1]]@ss)
    message("Done!\n")
  }

  dataset@is_harmonised <- TRUE
  return(dataset)
}
)


#' Check if the DataSet is harmonised
#' 
#' @param dataset A gwasglue2 DataSet object
#' @return TRUE/FALSE
#' @export 
#' @docType methods
#' @rdname isHarmonised-methods
setGeneric("isHarmonised",function(dataset) standardGeneric("isHarmonised"))
#' @rdname isHarmonised-methods
setMethod("isHarmonised","DataSet",function(dataset) {
  return(dataset@is_harmonised)
}
)

