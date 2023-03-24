#' Merge Datasets 
#'
#' @param ... A list of DataSet S4 objects
#' @return  A DataSet S4 object  with input DataSets merged
merge_datasets <- function(datasets) {
    n_datasets <- length(datasets)
	  ds <- DataSet()
    count <- 1
    for(i in 1:n_datasets){
      if (isS4(datasets[[i]]) == FALSE){
         next
         }
       n_ss <- length(datasets[[i]]@summary_sets)
       ds@ld_matrix <- datasets[[i]]@ld_matrix
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
    
    # ds@is_harmonisedLD <- datasets[[1]]@is_harmonisedLD
    # ds@zscores <- datasets[[1]]@zscores
    # ds@is_converted <- datasets[[1]]@is_converted

    return(ds)
}






# @title Utility function to display warning messages as they occur (from susieR)
# @param ... warning message
# @param style either "warning" or "hint"
#'@importFrom crayon combine_styles
warning_message = function(..., style=c("warning", "hint")) {
  style = match.arg(style)
  if (style=="warning" && getOption("warn")>=0) {
    alert <- combine_styles("bold", "underline", "red")	
    message(alert("WARNING:"), " ", ...)
  } else {
    alert <- combine_styles("bold", "underline", "magenta")	
    message(alert("HINT:"), " ", ...)
  }
}
