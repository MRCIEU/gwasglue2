# 28/07/2022

library(magrittr)
library(dplyr)
library(ieugwasr)

######################################################################################
# GWASGLUE2 FUNCTIONS

# CLUMPED HITS
#  @param clump Inherited from ieugwasr
#     If clump=FALSE, extracts the data without clumping
# @param source By default, OpenGWAS, Inherited from ieugwasr. 
#     access_token Google OAuth2 access token. Used to authenticate level of access to data.
# @return snps Array of variants
clumped_hits <- function(traits,  clump = TRUE,source = ieugwasr::check_access_token()){
    # ch  <- ieugwasr::tophits(traits[1], clump=clump, source=source) 
        
    #use this line code instead if there is an "unused argument (source = source)" error. It happens if the access was already opened
    ch  <- ieugwasr::tophits(traits[1], clump=clump)
    snps <- ch$rsid
    return(snps)
}

# PHEWAS_ids
# @param variants 
# @param variants Array of variants 
# @param pval p-value threshold. Default = `0.00001`. Iherited from ieugwasr::phewas
# @param batch Vector of batch IDs to search across. If `c()` (default) then returns all batches. Iherited fromieugwasr::phewas
# @return id_list Array with traits ids
phewas_ids <- function(variants, batch, pval){
    id_list <- ieugwasr::phewas(pval = pval, variants=variants, batch=batch)$id %>%
        unique()
    return(id_list)
}

# CREATE DATASET
# @param traits An array with the study ids (eg. "ieu-a-2")
# @param variants An array of variants to be analysed
# @return dataset List of dataframes of summary statistics. 
create_dataset <- function (traits, variants){
       data_set <- list() 
   for (i in 1:length(traits)) {

        data_set[[i]] <-  ieugwasr::associations(variants=variants, id=traits[i])
        names(data_set)[i] <-  unique(data_set[[i]]$id)
    }
    # source ("harmonise.R") #wip: this is modified from the harmonise function of TwoSampleMR
        # harmonise_data(data1,data2)
    return(data_set)
}



# ANNOTATE function (add extra columns to data_set)
# @param data_set Output of create_dataset function()
# @params exposure and outcome  For now it must be the dataset$id
# @param ld_ref Name of the reference dataset to be used to build the LD matrix
    # (eg "EUR" on OpenGWAS or if using locally corresponds to the prefix of the plink files)
# @param Type of analysis  Label to be used the the convert tools functions (eg mr, finemap, coloc, etc) 
# @return data_set A list of annotated dataframes
annotate <- function (data_set, exposure=NULL, outcome=NULL, analysis_type=NULL, ld_ref=NULL){
    # sweep through the data_sets
    for (i in 1:length(data_set)) {

        # specific for MR analysis
        if (!is.null(exposure) & !is.null(outcome)){
            # we named previously each dataset with $id (see create_dataset ())
            if (names(data_set)[i] == exposure){
                data_set[[i]]$exposure <- exposure
            } else {
                data_set[[i]]$outcome <- outcome
            }
        }

        # LD reference data_set name 
        if (!is.null(ld_matrix)){
        data_set[[i]]$ld_ref <- ld_ref
        }

        # Type of analysis
        if (!is.null(analysis_type)){
        data_set[[i]]$analysis_type <- analysis_type
        }
     }
    return(data_set)
}
