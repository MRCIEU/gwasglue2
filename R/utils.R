#' Merge Datasets 
#'
#' @param ... A list of DataSet S4 objects
#' @return  A DataSet S4 object  with input DataSets merged
merge_datasets <- function(...) {
   	datasets <- list(...)
    n_datasets <- length(datasets)
	ds <- DataSet()
    count <- 1
    for(i in 1:n_datasets){
        n_ss <- length(datasets[[i]]@summary_sets)
        for(j in 1:n_ss){
          ds@summary_sets[[count]] <- datasets[[i]]@summary_sets[[j]]
          count <- count + 1
        }
    }
	# For now I am assuming that in all datasets the following slots are the same
    # ds@overlap_SNPs <- datasets[[1]]@overlap_SNPs
    # ds@is_resized <- datasets[[1]]@is_resized
    # ds@is_harmonised <- datasets[[1]]@is_harmonised
    # ds@overall_dropped_SNPs <- datasets[[1]]@overall_dropped_SNPs
    # ds@dropped_SNPs <- datasets[[1]]@dropped_SNPs
    # ds@palindromic_SNPs <- datasets[[1]]@palindromic_SNPs
    # ds@ambiguous_SNPs <- datasets[[1]]@ambiguous_SNPs
    # ds@incompatible_alleles_SNPs <- datasets[[1]]@incompatible_alleles_SNPs
    ds@ld_matrices <- datasets[[1]]@ld_matrices
    # ds@is_harmonisedLD <- datasets[[1]]@is_harmonisedLD
    # ds@zscores <- datasets[[1]]@zscores
    # ds@is_converted <- datasets[[1]]@is_converted

    return(ds)
}



