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
    





# To make the summary set

# - chr, pos, a1, a2, n, eaf, etc are all the same as the
# - beta = from lbf_variable
# - se = from lbf_variable


#' Convert log Bayes Factor to summary stats
#'
#' @param lbf p-vector of log Bayes Factors for each SNP
#' @param n Overall sample size
#' @param af p-vector of allele frequencies for each SNP
#' @param prior_v Variance of prior distribution. SuSiE uses 50
#'
#' @return tibble with lbf, af, beta, se, z
#' @export 
lbf_to_z_cont <- function(lbf, n, af, prior_v=50){
  se = sqrt(1 / (2 * n * af * (1-af)))
  r = prior_v / (prior_v + se^2)
  z = sqrt((2 * lbf - log(sqrt(1-r)))/r)
  beta <- z * se
  return(cbind(lbf, af, z, beta, se))
}



create_summary_set_from_lbf <- function(summaryset, lbf, L){
  n <- summaryset@ss$n
  af <- summaryset@ss$eaf
  
  lbf_conv <- lbf_to_z_cont(lbf, n, af)
  
   # replace the beta and se columns in summaryset
  summaryset@ss$beta <- bf_conv$beta
  summaryset@ss$se <-bf_conv$se

 
  # update metadata to explain which credible set this is
  summaryset@metadata$id <- paste0(summaryset@metadata$id, "_L",i)
  # - trait name?
  # - id?
  # - notes?
}



