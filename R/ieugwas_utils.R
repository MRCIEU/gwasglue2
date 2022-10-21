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
