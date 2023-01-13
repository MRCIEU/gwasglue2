# This tutorial is based on vignettes/Tutorial1.Rmd

## 2. Regional genotype-phenotype map

# We need to first define a genomic range for the analysis to be performed. This is a known LD block on chromosome 22 (build hg19)
#
# ```{r}
# genomicrange <- "22:43714200-44995307"
#
# ```
#
# We need to scan OpenGWAS for traits that have a significant associations in this region
#
# ```{r}
# id_list <- ieugwasr::phewas(pval = 5e-6, variants=genomicrange, batch="ukb-a")$id %>%
#   unique()
# ```
#
# Obtain data
#
# ```{r}
# dataset <- gwasglue2::create_dataset(traits=id_list, variants=genomicrange) %>%
#   gwasglue2::annotate(ld_matrix="1000genomes_EUR")
# ```
#
# Next, for each trait we need to perform finemapping. This will involve using SusieR to identify credible sets in the region for the trait, and then generate the conditional summary datasets for each credible set. Explanation:
#
#   Credible set: If a trait is inferred to have e.g. 3 causal variants in the genomic region then there will be 3 credible sets. Each credible set is essentially a cluster of SNPs close by, it will include 1 or more SNPs that are candidates to be the causal variant.
#
# Conditional summary datasets: Once credible sets are identified, we can try to estimate what the summary data results would have looked like for each credible set after accounting for all the other credible sets. So in this example there would be 3 conditional summary datasets, each with a single peak in a different location.
#
# Use susieR to generate the credible sets for each dataset (e.g. some of them will have only 1, some will have more than 1. But all of them will have at least 1 because they have been selected to have an association in the region)
#
# ```{r}
# dataset <- gwasglue2::annotate(dataset = ., credible_sets="susieR")
# ```
#
# Finally, we can perform the colocalisation. Here we want to jointly colocalise across all the conditional summary datasets
#
#
# ```{r}
# result <- dataset %>%
#   hyprcoloc::hyprcoloc(dataset=.)
# ```


library(dplyr)


#install.packages("Rfast")  # optional, for better performance of susieR
#devtools::install_github("cnfoley/hyprcoloc", build_opts = c("--resave-data", "--no-manual"), build_vignettes = TRUE)
library(susieR)
library(hyprcoloc)

ieugwasr::api_status()

## 2. Regional genotype-phenotype map

# Obtain the data

ids <- ieugwasr::phewas(pval = 5e-6, variants = "22:43714200-44995307", batch = "ukb-a")$id %>%
  unique()

# create s4 SummarySet objects (named summary_set1...n) and fill metadata slot using the createSummarySets() in ieugwas_utils

for (i in seq_along(ids)){
  s <- createSummarySets(traits=ids[i],
                          variants = "22:43714200-44995307",
                          tools = c("finemap","coloc"),
                          source = "IEUopenGWAS",
                          ld_ref = "../data/reference/ld/EUR")
  assign(paste0("summary_set", i) , s)

}




# create S4 DataSet object
dataset <- DataSet(summary_set1,summary_set2,summary_set3, summary_set4) %>%
  overlapSNP(.) %>%
  harmoniseData(.,tolerance = 0.08,action = 1)


d <- buildLDMatrix(dataset, bfile = TRUE, plink_bin = "plink")
d1 <- harmoniseLDMatrix(d)  
%>%
d2 <- setZscores(d1)

res_susie <- run_susieR(d2, ids)


#' run susieR
#' @param data object class DataSet
#' @param ids array with trait ids
#' @return credible sets and fmset
run_susieR <- function (data,ids){
  res_f <- list()
  for (i in seq_along(ids)) {
    res <- susieR::susie_rss(
      z = getZscores(data, index = i),
      R = getLDMatrix(data,i),
      n = unique(getData(data, i)$n),
      estimate_residual_variance = FALSE,
      coverage=0.95 # default 0.95
    )
    res$fmset <- sapply(res$sets$cs, function(x) {
      getData(data,i)$rsid[x[which.max(res$pip[x])]]
    })
    res_l <- list()
    if(is.null(res$sets$cs)) {
      res_l[["credible_sets"]] <- NA
      res_l[["fmset"]] <- NA
    }else {
      res_l[["credible_sets"]] <-res$sets$cs
      res_l[["fmset"]]  <- res$fmset
    }
    res_f <- append(res_f,res_l)
    names(res_f)[i] <-ids[i]
    message(paste("susieR analyses for trait", ids[i], "finished"))
  }
  
  return(res_f)

}


#' CONVERT COLOC
#' @param data_set Output of annotate function
#' @return data_set A list of annotated dataframes (Beta and se matrices for colocalisation)
convert_coloc <- function(data_set){
  for (i in seq_along(data_set)) {

    if ( is.na(data_set[[i]]$credible_sets[1])) {
      data_set[[i]]$beta_coloc <- NA
      data_set[[i]]$se_coloc <- NA
    }

    else{
      beta_coloc <- list()
      se_coloc <- list()

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
        beta_coloc[[j]] <- beta_mat
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
  message(paste("Colocalisation convertion for trait", i, "finished"))
}



######################################################################################
# plink binaries can be found here: https://www.cog-genomics.org/plink/

library(magrittr)
library(dplyr)


path_to_plink <- "plink"
# Obtain the data
dt1 <- phewas_ids(variants = "22:43714200-44995307", batch = "ukb-a", pval = 5e-6) %>%

  dt1 <- phewas_ids(variants = "22:43714200-44995307", batch = "ukb-b", pval = 5e-6) %>%
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
  if (is.na(dt3[[i]]$beta_coloc[1])) next
  else {
    # sweep L regions
    for(j in seq_along(dt3[[i]]$beta_coloc)){
      print(paste0("L",j))
      traits <- colnames(dt3[[i]]$beta_coloc[[j]])
      rsid <- rownames(dt3[[i]]$beta_coloc[[j]])
      betas <- dt3[[i]]$beta_coloc[[j]]
      ses <- dt3[[i]]$se_coloc[[j]]
      res <- hyprcoloc::hyprcoloc(betas, ses, trait.names=traits, snp.id=rsid,bb.selection = "regional")
      print(res)
    }
  }
}
