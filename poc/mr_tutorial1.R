# this proof-of-concept is based on the Tutorial from the 22nd of June 2022

# 28/07/2022

# call gwasglue2 functions
source("functions_gwasglue2.R")

# INSTALL ieugwasr and TwoSampleMR packages
install.packages("remotes")
remotes::install_github("mrcieu/ieugwasr")
remotes::install_github("MRCIEU/TwoSampleMR")


library(magrittr)
library(dplyr)
library(ieugwasr)
library(TwoSampleMR)


ieugwasr::api_status()

######################################################################################
# GWASGLUE2 MR FUNCTIONs 


# CHANGE COLNAMES 
# (makes it easier to work in MR: it is inside create_dataset, but it can be moved to convert_mr)
change_colnames <- function (data_set){
d <- data_set %>% rename(c(SNP = rsid,effect_allele=ea, other_allele = nea, pval=p, samplesize=n))
return (d)
}


# CONVERT MR function Converts data_set in a format to be used in MR
# @param data_set list of  dataframes
# @return dataset Modified list of  dataframes to be used in TwoSampleMR
convert_mr <- function(data_set){
    for (i in 1:length(data_set)) {
        # change column names
        data_set[[i]] <- change_colnames(data_set[[i]])
        #  add extra column named "mr_keep"
        data_set[[i]][,"mr_keep"] <- !is.na(data_set[[i]]$beta) & !is.na(data_set[[i]]$se) & !is.na(data_set[[i]]$effect_allele) & !is.na(data_set[[i]]$SNP)

        # add sufix ".exposure" or ".outcome" to all column names except to  columns SNP and label
        if("exposure" %in% colnames(data_set[[i]])){
            c <-which(colnames(data_set[[i]])!="SNP" & colnames(data_set[[i]])!="exposure")
            colnames(data_set[[i]])[c] <- paste0(colnames(data_set[[i]])[c], ".exposure")
        } else{
            c <-which(colnames(data_set[[i]])!="SNP" & colnames(data_set[[i]])!="outcome")
            colnames(data_set[[i]])[c] <- paste0(colnames(data_set[[i]])[c], ".outcome")
        }
    }
    return(data_set)
}

# CONVERT TOOL function Converts data_set in different formats depending of the method
# @param data_set list of  dataframes
# @return dataset Modified list of  dataframes
convert_tool <- function(data_set){
    # need to work more on this function to read all datasets and not just the first
    if (unique(data_set[[1]]$analysis_type) == "mr"){
        convert_mr(data_set)
    }
}


#####################################################################################
# TwoSampleMR analysis

# We will perform an MR analysis of body mass index (exposure - ieu-a-2) against coronary heart disease (outcome - ieu-a-7).

# call gwasglue2 functions
source("functions_gwasglue2.R")

# Obtain the data
x<-clumped_hits(traits="ieu-a-2")  %>%
    create_dataset(traits=c("ieu-a-2","ieu-a-7"),variants=.) %>%
    annotate (.,exposure = "ieu-a-2",outcome ="ieu-a-7",analysis_type ="mr") %>% 
    convert_tool(.)

# Perform the MR analysis (with harmonisation)
mr_result<- TwoSampleMR::harmonise_data(x[[1]],x[[2]])  %>%
    TwoSampleMR::mr(., method_list="mr_ivw")

mr_result
