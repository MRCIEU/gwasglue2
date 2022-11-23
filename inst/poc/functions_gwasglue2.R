# 28/07/2022  created
# 29/07/2022
    # annotate function updated to allow more than one type of analysis_type,
    # 1:length() replaced by seq_along()
    # some spaces and lint solved
    #  using Roxygen comments for documenting functions
# 9/08/2022
    # changes to annotate function (create LD matrices)


library(magrittr)
library(dplyr)
library(ieugwasr)

######################################################################################
# GWASGLUE2 FUNCTIONS

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
#' @return dataset List of dataframes of summary statistics.
create_dataset <- function(traits, variants) {
       data_set <- list()
   for (i in seq_along(traits)) {

        data_set[[i]] <-  ieugwasr::associations(variants = variants, id = traits[i])
        names(data_set)[i] <-  unique(data_set[[i]]$id)
    }
    # source ("harmonise.R") #wip: this is modified from the harmonise function of TwoSampleMR::mv_harmonise_data()
    # pairwise? comparing to one  or all combinations?
        # harmonise_data(data1,data2)
    return(data_set)
}



#' ANNOTATE function (add extra columns to data_set)
#' @param data_set Output of create_dataset function()
#' @params exposure and outcome  For now it must be the dataset$id
#' @param ld_ref Name of the reference dataset to be used to build the LD matrix
    #' @param pop TRUE/FALSE
    #' Default FALSE. If TRUE the pop parameter will be used in ieugwasr::ld_matrix
    #' @param bfile TRUE/FALSE
    #' Default TRUE. If TRUE the bfile parameter will be used in ieugwasr::ld_matrix
    #' @param plink_bin Path to where the plink executable is in the OS
#' (eg "EUR" on OpenGWAS or if using locally corresponds to the prefix of the plink files)
#' @param analysis_type  Array of labels used for the convert tools functions (eg mr, finemap, coloc, etc)
#'  It creates the same number of columns as the length of the array (Eg, analysis_type1, analysis_type2, etc)
#' @return data_set A list of annotated dataframes
annotate <- function(data_set, exposure = NULL, outcome = NULL, analysis_type = NULL, ld_ref = NULL, pop = FALSE, with_alleles = FALSE, bfile = TRUE, plink_bin = NULL) {
    # sweep through the data_sets
    for (i in seq_along(data_set)) {

        # specific for MR analysis
        if (!is.null(exposure) && !is.null(outcome)) {
            # we named previously each dataset with $id (see create_dataset ())
            if (names(data_set)[i] == exposure) {
                data_set[[i]]$exposure <- exposure
            } else {
                data_set[[i]]$outcome <- outcome
            }
        }

        # LD reference data_set name
        if (!is.null(ld_matrix)) {
            # TODO
            #   LD matrices (possible for each trait to have a different LD matrix, but not multiple LD matrices per trait) : ld_ref as to be an array eg rep("EUR",length(data_set))
            # read the ld_ref column for the ref name (for now it is just an array of length 1)
            if(length(ld_ref) == 1) {
                ldr <- ld_ref 
            }
            else{
                ldr <- ld_ref[[i]]
            }

            message("Building LD matrix")
            if (pop == TRUE) {
                ld <- suppressWarnings(ieugwasr::ld_matrix(data_set[[i]]$rsid, pop = ldr, with_alleles = with_alleles, bfile = FALSE,plink_bin = plink_bin))
            }
            if (bfile == TRUE){
                ld <- suppressWarnings(ieugwasr::ld_matrix(data_set[[i]]$rsid, pop=FALSE, with_alleles=with_alleles, bfile=ldr,plink_bin=plink_bin))   
            }
            rsid_avail <- rownames(ld)
            message("Data available for ", length(rsid_avail), " variants")
            # subseting
            data_set[[i]] <- subset(data_set[[i]], data_set[[i]]$rsid %in% rsid_avail)
            index <- match(data_set[[i]][["rsid"]], rsid_avail)

            ld <- ld[index, index]
            stopifnot(all(data_set[[i]][["rsid"]] == rownames(ld)))

            # n <- data_set[[i]][["n"]]
            # if(all(is.na(n))){
            # 	g <- ieugwasr::gwasinfo(id[i])
            # 	n <- g[["sample_size"]]

            # harmonise data
            harmonise_ld_dat(data_set[[i]],ld)

            # add objects to the data_set list
            data_set[[i]][["ld"]]  <- ld

            
        }

        # Type of analysis
        if (!is.null(analysis_type)) {
            for (j in seq_along(analysis_type)) {
            data_set[[i]][, paste0("analysis_type", j)] <- analysis_type[j]
            }
        }
     }
    return(data_set)
}
