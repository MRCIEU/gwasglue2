# Wrapper functions for other tools


#' DataSet to hyprcoloc
#' 
#' dataset_to_hyprcoloc is a wrapper function used inside ritarasteiro/hyprcoloc::hyprcoloc() to read DataSet objects
#' @param dataset gwasglue2 DataSet object
#' @return parameters needed to run hyprcoloc
#' @export 
dataset_to_hyprcoloc <- function(dataset){

	ntraits <- length(dataset@summary_sets)
 	trait.names <- 	unlist(lapply(1:ntraits, function(i){
		t <- dataset@summary_sets[[i]]@metadata$id
	}))

    snp.id <- dataset@summary_sets[[1]]@ss$rsid 
    ld.matrix <- dataset@ld_matrices[[1]]
    effect.est <- matrix(ncol = length(dataset@summary_sets),nrow=length(snp.id))
    effect.se <- matrix(ncol = length(dataset@summary_sets),nrow=length(snp.id))
    
    for (i in seq_along(trait.names)){
      effect.est[,i] <- dataset@summary_sets[[i]]@ss$beta
      effect.se[,i] <- dataset@summary_sets[[i]]@ss$se
    }
return(list(trait.names, snp.id, ld.matrix, effect.est, effect.se))
}


summaryset_to_susieR <- function(summaryset){

}

susieR_to_dataset <- function(){



}
    