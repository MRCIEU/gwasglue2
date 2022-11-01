# ieugwasr Functions

#' A function to clump the top hit variants
#'
#' @param traits
#' @param clump
#' @param source
#'
#' @return
clumpTophits <- function(traits,  clump = TRUE, source = ieugwasr::check_access_token()) {
  ch  <- ieugwasr::tophits(traits[1], clump = clump)
  snps <- ch$rsid
  return(snps)
}

#' PHEWAS_ids
#' @param variants Array of variants
#' @param pval p-value threshold. Default = `0.00001`. Iherited from ieugwasr::phewas
#' @param batch Vector of batch IDs to search across. If `c()` (default) then returns all batches. Iherited fromieugwasr::phewas
#' @return id_list Array with traits ids
phewas_ids <- function(variants, batch, pval) {
  id_list <- ieugwasr::phewas(pval = pval, variants = variants, batch = batch)$id %>%
    unique()
  return(id_list)
}

#' A function to call the summary statistics associated with traits and variants chosen
#'
#' @param traits
#' @param variants
#'
#' @return
createSumset <- function(traits, variants) {
  x <-  ieugwasr::associations(variants = variants, id = traits)
  return(x)
}
