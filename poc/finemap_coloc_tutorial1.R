# this is a proof-of-concept based on a Tutorial created by Gib on the 22nd of June 2022 

# call gwasglue2 functions
source("functions_gwasglue2.R")

## 2. Regional genotype-phenotype map

install.packages("Rfast")  # optional, for better performance of susieR
devtools::install_github("cnfoley/hyprcoloc", build_opts = c("--resave-data", "--no-manual"), build_vignettes = TRUE)
library(susieR)

ieugwasr::api_status()



#' CONVERT FINE MAPPING
#' @param data_set Output of create_dataset function

#' @return data_set A list of annotated dataframes with extra elements for each trait id (z$snp,z$zscore ld, n )
convert_finemap <- function(data_set) {
   for (i in seq_along(data_set)) {
        # id <- unique(x$id)
        dat <- list()
        # calculate some sumstats
		dat[["z"]] <- dplyr::tibble(snp = data_set[[i]][["rsid"]], zscore = data_set[[i]][["beta"]] / data_set[[i]][["se"]])
		
		# add objects to the data_set list
         dat_l <- append(data_set[[i]],dat)
         data_set[[i]] <- dat_l
        }
return(data_set)
}

#' run susieR
#' @param data_set 
#' @return credible sets and fmset added to data_set  
run_susieR<- function (data_set){
    for (i in seq_along(data_set)) {
            res <- susieR::susie_rss(
                z= data_set[[i]]$z$zscore, 
                R=data_set[[i]]$ld, 
                n= unique(data_set[[i]]$n),
                estimate_residual_variance = TRUE, 
                coverage=0.95 # default 0.95
            )
            res$fmset <- sapply(res$sets$cs, function(x) {
                data_set[[i]]$z$snp[x[which.max(res$pip[x])]]
                })
            res_l <- list()
            if(is.null(res$sets$cs)) {
                res_l[["credible_sets"]] <- NA
                res_l[["fmset"]] <- NA
            }else {
                res_l[["credible_sets"]] <-res$sets$cs
                res_l[["fmset"]]  <- res$fmset
                }
            dat_l <- append(data_set[[i]],res_l)
            data_set[[i]] <- dat_l
    }
    return(data_set)
}


#' CONVERT COLOC
#' @param data_set Output of annotate function
#' @return data_set A list of annotated dataframes (Beta and se matrices for colocalisation)
convert_coloc <- function(data_set){
    for (i in seq_along(data_set)) {

        if (is.na(data_set[[i]]$credible_sets)) {
            data_set[[i]]$beta_coloc <- NA
            data_set[[i]]$se_coloc <- NA
            } 

    else{
        for(j in seq_along(data_set[[i]]$credible_sets)){
            pos <- data_set[[i]]$credible_sets[paste0("L",j)]
            pos <-  pos[[paste0("L",j)]]

            beta_mat <- matrix(,ncol=length(data_set),nrow=length(pos))
            rownames(beta_mat)<-data_set[[i]]$rsid[pos]
            colnames(beta_mat)<-  names(data_set)
            se_mat <- matrix(,ncol=length(data_set),nrow=length(pos))
            rownames(se_mat)<-data_set[[i]]$rsid[pos]
            colnames(se_mat)<-  names(data_set)
            
            c=1
            for (z in seq_along(data_set)){
                beta_mat[,c] <- data_set[[z]]$beta[pos]
                se_mat[,c] <- data_set[[z]]$se[pos]
                c=c+1
            } 

            beta_coloc <- list()
            beta_coloc[[j]] <- beta_mat

            se_coloc <- list()
            se_coloc[[j]] <- se_mat
        }
        coloc <- list()
        coloc[["beta_coloc"]] <- beta_coloc
        coloc[["se_coloc"]] <- beta_coloc
        dat_c <- append(data_set[[i]],coloc)
        data_set[[i]] <- dat_c  
    }
    }
    return(data_set)
}



#' CONVERT TOOL function Converts data_set in different formats depending of the method
#' TODO NEEDS more work to include params from the several convert functions
#' @param data_set list of  dataframes
#' @param data_set analysis_order (?) (corresponds to analysis_type columns, eg nb_analyses=1 will read analysis_type1)
#' Usefull when we need to convert data in different steps
#' @return dataset Modified list of  dataframes
convert_tool <- function(data_set, analysis_order = 1, pop=FALSE, with_alleles=FALSE, bfile=TRUE,plink_bin=NULL) {
    # need to work more on this function to read all datasets and not just the first
    if (unique(data_set[[1]][paste0("analysis_type",analysis_order)]) == "mr") {
        message("Converting to mr format")
        convert_mr(data_set)
    }
    if (unique(data_set[[1]][paste0("analysis_type",analysis_order)]) == "finemap"){
        message("Converting to finemap format")
        convert_finemap(data_set,pop = pop, with_alleles = with_alleles, bfile = bfile, plink_bin =plink_bin)
    }
    if (unique(data_set[[1]][paste0("analysis_type",analysis_order)]) == "coloc"){
        message("Converting to coloc format")
        convert_coloc(data_set)
    }
    return(data_set)
}

######################################################################################
# ANALYSIS (I am working inside the poc folder)
# call gwasglue2 functions
source("functions_gwasglue2.R")

# The EUR plink files can be downloaded here (multiple populations):
# http://fileserve.mrcieu.ac.uk/ld/1kg.v3.tgz

# plink binaries can be found here: https://www.cog-genomics.org/plink/

library(magrittr)
library(dplyr)
source("functions_gwasglue2.R")

path_to_plink <- "plink"
# Obtain the data
dt1 <- phewas_ids(variants = "22:43714200-44995307", batch = "ukb-a", pval = 5e-6) %>%
    create_dataset(traits = ., variants = "22:43714200-44995307") %>%
    annotate(data_set = ., ld_ref="EUR", analysis_type = c("finemap","coloc"), bfile = TRUE, plink_bin = path_to_plink)

# convert to finemap format and run susieR
dt2<- dt1 %>%
    convert_finemap(data_set =.) %>% 
    run_susieR(data_set = .)

# convert to coloc format
dt3 <- dt2 %>%
     convert_coloc(data_set=.)

# coloc analysis (TODO create a function to run hypercoloc)
# assuming independent traits
for (i in seq_along(dt3)){
    if (is.na(dt3[[i]]$beta_coloc)) next
    else {
        # sweep L regions
        for(j in seq_along(dt3[[i]]$beta_coloc)){
        traits <- colnames(dt3[[i]]$beta_coloc[[j]])
        rsid <- rownames(dt3[[i]]$beta_coloc[[j]])
        betas <- dt3[[i]]$beta_coloc[[j]]
        ses <- dt3[[i]]$se_coloc[[j]]
        res <- hyprcoloc::hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid,bb.selection = "regional")
        print(res)
        }
    } 
}
