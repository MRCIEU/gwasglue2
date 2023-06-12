# Wrapper functions for other tools


#' DataSet to hyprcoloc
#' 
#' dataset_to_hyprcoloc is a wrapper function used inside ritarasteiro/hyprcoloc::hyprcoloc() to read DataSet objects
#' @param dataset gwasglue2 DataSet object
#' @return parameters needed to run hyprcoloc
#' @export 
dataset_to_hyprcoloc <- function(dataset){
  message("hyprcoloc is using gwasglue2 DataSet class object as input")
	ntraits <- length(dataset@summary_sets)
 	trait.names <- 	unlist(lapply(1:ntraits, function(i){
		t <- dataset@summary_sets[[i]]@metadata$id
	}))

    snp.id <- dataset@summary_sets[[1]]@ss$variantid
    ld.matrix <- dataset@ld_matrix
    effect.est <- matrix(ncol = length(dataset@summary_sets),nrow=length(snp.id))
    effect.se <- matrix(ncol = length(dataset@summary_sets),nrow=length(snp.id))
    
    for (i in seq_along(trait.names)){
      effect.est[,i] <- dataset@summary_sets[[i]]@ss$beta
      effect.se[,i] <- dataset@summary_sets[[i]]@ss$se
    }
return(list(trait.names, snp.id, ld.matrix, effect.est, effect.se))
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
lbf_to_z_cont <- function(lbf, n, af, prior_v = 50){
  se = sqrt(1 / (2 * n * af * (1-af)))
  r = prior_v / (prior_v + se^2)
  z = sqrt((2 * lbf - log(sqrt(1-r)))/r)
  beta <- z * se
  return(data.frame(lbf, af, z, beta, se))
}


#' Create SummarySet from log Bayes Factor
#'
#' @param summaryset gwasglue2 SummarySet object
#' @param lbf p-vector of log Bayes Factors for each SNP
#' @param L credible set index number
#' @return marginalised summaryset (beta, se and trait id)
#' @export
#' 
create_summary_set_from_lbf <- function(summaryset, lbf, L){
  if (is.null(getMetadata(summaryset)$sample_size) || is.na(getMetadata(summaryset)$sample_size)){
      stop("No sample size information in metadata for this SummarySet. More details on how to add to metadata in 'help(gwasglue2::create_metadata)' and 'help(gwasglue2::getMetadata)'." )
  }
  af <- summaryset@ss$eaf
  n <- summaryset@metadata$sample_size
  
  lbf_conv <- lbf_to_z_cont(lbf, n, af)
   # replace the beta and se columns in summaryset
  summaryset@ss$beta <- lbf_conv$beta
  summaryset@ss$se <- lbf_conv$se

 
  # update metadata to explain which credible set this is
  summaryset@metadata$id <- paste0(summaryset@metadata$id, "_L",L)
  # - trait name?
  # - id?
  # - notes?
  return(summaryset)
}

#' SusieR to DataSet
#'
#' #' susie_to_dataset is a wrapper function used inside ritarasteiro/susieR::susie_rss() to create a marginalised DataSet object
#' Converts ABFs to summary statistics and creates a new SummarySet for each credible set. Returns object is a gwasglue2 DataSet class object. 
#' @param summaryset gwasglue2 SummarySet object
#' @param s susieR object
#' @param R lD matrix
#' @return DataSet object
#' @export 
#' 
susie_to_dataset <- function(summaryset, s, R){
 ncredible_sets <- length(s$sets$cs)
    if(ncredible_sets == 0){
      warning_message(paste0("There is no credible sets for this trait (",summaryset@metadata$id,"), with the parameter values used. The summary statistics beta and se will not be marginalised."))

    } else{
       ds <- DataSet(list())
      
      for(i in 1:ncredible_sets){
        ds@summary_sets[[i]] <- create_summary_set_from_lbf(summaryset, s$lbf_variable[i,], L = i)
        ds@ld_matrix <- R
      }

      ds@susie_marginalised <- TRUE

    # ds@susieR <- s
    s <- ds
    
    }
   return(s)
}


