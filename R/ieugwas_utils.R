# ieugwasr Functions

#' A function to clump the top hit variants
#'
#' @param traits ID of GWAS studies to query
#' @param clump Logical (default TRUE)
#' @param source  IEU (default)
#'
#' @return Array of SNPs rsids
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
phewasIDs <- function(variants, batch, pval) {
  id_list <- unique(ieugwasr::phewas(pval = pval, variants = variants, batch = batch)$id) 
  return(id_list)
}

#' A function to call the summary statistics associated with specific traits and variants chosen
#'
#' @param traits ID of GWAS studies to query
#' @param variants Array of SNPs rsids
#'
#' @return A tibble with the summary statistics for variant association.
createSumset <- function(traits, variants) {
  x <-  ieugwasr::associations(variants = variants, id = traits)
  return(x)
}
