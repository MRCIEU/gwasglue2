# this is a proof-of-concept based on a Tutorial created by Gib on the 22nd of June 2022

## 1. TwoSampleMR analysis

# We will perform an MR analysis of body mass index (exposure - ieu-a-2) against coronary heart disease (outcome - ieu-a-7).

# Obtain the data

# ```{r}
# dataset <- gwasglue2::clumped_hits(traits = "ieu-a-2", source = "OpenGWAS") %>%
#   gwasglue2::create_dataset(traits = c("ieu-a-2", "ieu-a-7"), variants = .) %>%
#   gwasglue2::annotate(exposure = "ieu-a-2", outcome = "ieu-a-7")
# ```

# Perform the analysis

# ```{r}
# result <- dataset %>%
#   TwoSampleMR(dataset=., method_list="mr_ivw")
# ```

# INSTALL ieugwasr and TwoSampleMR packages
install.packages("remotes")
remotes::install_github("mrcieu/ieugwasr")
remotes::install_github("MRCIEU/TwoSampleMR@0.4.26")


library(magrittr)
library(dplyr)
library(ieugwasr)
library(TwoSampleMR)


ieugwasr::api_status()

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

# CREATE DATASET
# @param traits An array with the study ids (eg. "ieu-a-2")
# @param variants An array of variants to be analysed
# @return dataset A list of dataframes of summary statistics. ÷
create_dataset <- function (traits, variants){  
       data_set <- list() 
   for (i in 1:length(traits)) {

        data_set[[i]] <-  ieugwasr::associations(variants=variants, id=traits[i])
        data_set[[i]] <- change_colnames(data_set[[i]])
    }
    # source ("harmonise.R") #wip: this is modified from the harmonise function of TwoSampleMR
    # For now I am doing the harmonisation during the mr analyses(see below)
    return(data_set)
}




# CHANGE COLNAMES 
# (makes it easier to work in MR: it is inside create_dataset, but it can be moved to convert_mr())
change_colnames <- function (data_set){
d <- data_set %>% dplyr::rename(c(SNP = rsid,effect_allele=ea, other_allele = nea, pval=p, samplesize=n))
return (d)
}


# ANNOTATE function (For now, I am using the format for MR analysis)
# @param data_set Output of create_dataset function()
# @param labels Array with labels for each dataframe
#   add extra column named by labels
# @return data_set A list of annotated dataframes
annotate <- function (data_set, labels) {
    for (i in 1:length(data_set)) {
        #  add extra column named by labels (I am using a format I saw in one of the MR examples)
        data_set[[i]][,labels[i]] <- paste(data_set[[i]]$trait,"|| id:", data_set[[i]]$id)
        # but it wil also work if we use:
        # data_set[[i]][,labels[i]] <- labels[i]
     }
    return(data_set)
}


# CONVERT MR function Converts data_set in a format to be used in MR
# @param labels Array with labels for each dataframe 
#   add sufixes ".labels" to col names
#   add  mr_keep col
# @return dataset Modified list of  dataframes to be used in TwoSampleMR()
convert_mr <- function(data_set, labels){
    for (i in 1:length(data_set)) {
    #  add extra column named "mr_keep"
    data_set[[i]][,"mr_keep"] <- !is.na(data_set[[i]]$beta) & !is.na(data_set[[i]]$se) & !is.na(data_set[[i]]$effect_allele) & !is.na(data_set[[i]]$SNP)

    # add sufix ".label" to all column names except to  columns SNP and label
        c <-which(colnames(data_set[[i]])!="SNP" & colnames(data_set[[i]])!=labels[i])
        colnames(data_set[[i]])[c] <- paste0(colnames(data_set[[i]])[c], ".", labels[i])
    }
    return(data_set)
}

######################################################################################
# ANALYSIS

# Obtain the data
x<-clumped_hits(traits="ieu-a-2")  %>%
    create_dataset(traits=c("ieu-a-2","ieu-a-7"),variants=.) %>%
    annotate (.,labels=c("exposure","outcome"))%>% 
    convert_mr(.,labels=c("exposure","outcome"))


# Perform the MR analysis (with harmonisation)
mr_result<- TwoSampleMR::harmonise_data(x[[1]],x[[2]])  %>%
    mr(., method_list="mr_ivw")

mr_result




######################################################################################
# OTHER TESTS


#  OLD CREATE DATASET (For two traits, with clumped_hits inside)
# @param traits An array with the study ids (eg. "ieu-a-2")
# @param nb_traits Number of traits analised ( 1 or 2, eg. for MR analysis =2)
# @param clump Inherited from ieugwasr
#     When clump= TRUE corresponds to clummped_hits in vignette
#     traits[1] is always the one to be clumped (exposure in MR)
#     traits[2] is treated as outcome for MR analysis
#     If clump=FALSE, extracts the data without clumping
# @param source By default, OpenGWAS, Inherited from ieugwasr. 
#     access_token Google OAuth2 access token. Used to authenticate level of access to data.
# @return Dataframe of summary statistics. 
#     If nb_traits=2, list of two dataframes
create_dataset <- function (traits, nb_traits, clump = TRUE,source = ieugwasr::check_access_token()){  # nolint
        if (nb_traits==1) {
        # data1  <- ieugwasr::tophits(traits[1], clump=clump, source=source) 
        
        # used this line code instead if there is an "unused argument (source = source)"" error. It happens if the access was already opened
        data_set  <- ieugwasr::tophits(traits[1], clump=clump)
        
        data_set <- change_colnames(ds)
        return (data_set)
    }
    if (nb_traits==2) {       
        # data1  <- ieugwasr::tophits(traits[1], clump=clump, source=source) 
        
        #used this line code instead if there is an "unused argument (source = source)"" error. It happens if the access was already opened
        data1  <- ieugwasr::tophits(traits[1], clump=clump)
        
        snps <- data1$rsid
        data1 <- change_colnames(data1)
        data2 <- ieugwasr::associations(variants=snps, id=traits[2], proxies=0)
        data2 <- change_colnames(data2)

        # source ("harmonise.R") #this is modified from the harmonise function of TwoSampleMR
        # harmonise_data(data1,data2)

        data_set <- list(data1,data2)
        return(data_set)
    }
}



# OLD ANNOTATE function (For now, I am using the format for MR analysis)
# @param data_set Output of create_dataset function()
# @param nb_traits Number of traits analised ( 1 or 2, eg. for MR analysis =2)
# @param labels Array with labels for each dataframe (add ".labels" to col names)
# @return Annotated dataframe

OLD_annotate <- function (data_set, nb_traits, labels) {
    if (nb_traits == 1) {
        colnames(data_set) <- paste0(colnames(data_set), ".", labels)
        return(data_set)
    }
    if (nb_traits == 2) {
        colnames(data_set[[1]]) <- paste0(colnames(data_set[[1]]), ".", labels[1])
        colnames(data_set[[2]]) <- paste0(colnames(data_set[[2]]), ".", labels[2])
        return(data_set)
    }
}

