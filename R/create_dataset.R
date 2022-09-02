#' CLUMPED HITS
#' @param clump Inherited from ieugwasr
#'     If clump=FALSE, extracts the data without clumping
#' @param source By default, OpenGWAS, Inherited from ieugwasr. 
#'     access_token Google OAuth2 access token. Used to authenticate level of access to data.
#' @return snps Array of variants
clumped_hits <- function(traits,  clump = TRUE, source = ieugwasr::check_access_token()) {
    # ch  <- ieugwasr::tophits(traits[1], clump=clump, source=source) 
        
    #use this line code instead if there is an "unused argument (source = source)" error. It happens if the access was already opened
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

#' CREATE DATASET
#' @param traits An array with the study ids (eg. "ieu-a-2")
#' @param variants An array of variants to be analysed
#' @return data_set class(dataset) List of dataframes of summary statistics.
create_dataset <- function(traits, variants) {
       data_set <- list()
   for (i in seq_along(traits)) {

        data_set[[i]] <-  ieugwasr::associations(variants = variants, id = traits[i])
        names(data_set)[i] <-  unique(data_set[[i]]$id)
    }
    # source ("harmonise.R") #wip: this is modified from the harmonise function of TwoSampleMR
        # harmonise_data(data1,data2)
    # S3 class
     class(data_set) <- "dataset"
    return(data_set)
}
