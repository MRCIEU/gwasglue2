# YAML constructors


#' YAML constructor for the coloc job
#' @param traits A character vector of traits to be analysed. Each trait should be prefixed with the source of the data ("opengwas" if querring the IEU OpenGWAS database, "vcf", "text", "robject" and "summaryset") followed by their location. For "opengwas", the location will be the IEU OpenGWAS ID, for "vcf" and "text" will be the name of the files, for the "robject" will be the name of a R object containing a dataset dataframe and "summaryset"  be the name of a `SummarySet()` object. Eg. "opengwas:ieu-a-2", "vcf:ieu-a-2.vcf.gz", "text:ieu-a-2.txt", "robject:dt_ieu-a-2", "summaryset:sumset_ieu-a-2". Note that files should be in the work session directory.
#' @param region A character vector of regions to be analysed.
#' @param build  A character vector  of the reference genome assemblies to generate the vcf file. Default is "GRCh37". 
#' * Options are `"NCBI34"`, `"NCBI35"`, `"NCBI36"`, `"GRCh37"` or "GRCh38".
#' @param bfile It corresponds to the path and prefix of the plink files used to build the LD correlation matrix. 
#' @param plink_bin Path to the plink executable
#' @param write A boolean indicating whether the yaml object should be written to a file.
#' @param outfile A character vector indicating the path to write the yaml file.
#' @return A yaml object.
#' @importFrom yaml write_yaml as.yaml
#' @export 
yaml_coloc <- function(traits, region, build = "GRCh37", bfile = NULL, plink_bin = NULL, write = FALSE, outfile = "config.yaml"){
  
  # read number of DataSets/jobs to  build and analyse
  njobs <- length(region)
  
  out <- lapply(region, \(r){
    l <- list()
    l$variants <- list(shape = "single_region", variant_list = r)
    l$summarydata <- strsplit(traits, ":") %>% lapply(., \(x) list(source=x[1], location=x[2]))
    
    if(!is.null(build)){
      l$summarydata <- l$summarydata %>% lapply(., \(x) {x$build <- build; x})
    }
    if(!is.null(bfile) & !is.null(plink_bin)){
      l$lddata <- list(location = bfile, plink_bin=plink_bin)
    }
    
    return(l)
  })
    
  
  names(out) <- c(paste0("job", 1:njobs))
  out <- append(out, list(analyses=list(type = "coloc")))

  
  # write yaml to file
  if (isTRUE(write)){yaml::write_yaml(out, outfile)}

  # return a yaml object
  out <- yaml::as.yaml(out)
  return(out)
}   




#' YAML constructor for the LD scores job
#' @param traits A character vector of traits to be analysed. Each trait should be prefixed with the source of the data ("opengwas" if querring the IEU OpenGWAS database, "vcf", "text", "robject" and "summaryset") followed by their location. For "opengwas", the location will be the IEU OpenGWAS ID, for "vcf" and "text" will be the name of the files, for the "robject" will be the name of a R object containing a dataset dataframe and "summaryset"  be the name of a `SummarySet()` object. Eg. "opengwas:ieu-a-2", "vcf:ieu-a-2.vcf.gz", "text:ieu-a-2.txt", "robject:dt_ieu-a-2", "summaryset:sumset_ieu-a-2". Note that files should be in the work session directory.
#' @param variants A character vector of variants to be analysed.
#' @param build  A character vector  of the reference genome assemblies to generate the vcf file. Default is "GRCh37". 
#' * Options are `"NCBI34"`, `"NCBI35"`, `"NCBI36"`, `"GRCh37"` or "GRCh38".
#' @param write A boolean indicating whether the yaml object should be written to a file.
#' @param outfile A character vector indicating the path to write the yaml file.
#' @return A yaml object.
#' @importFrom yaml write_yaml as.yaml
#' @export 
yaml_ldscores <- function(traits, variants, build = "GRCh37", write = FALSE, outfile = "config.yaml"){
  # read number of DataSets/jobs to  build and analyse
  njobs <- 1
  
  out <- lapply(njobs, \(r){
    l <- list()
    l$variants <- list(shape = "scatterd", variant_list = variants)
    l$summarydata <- strsplit(traits, ":") %>% lapply(., \(x) list(source=x[1], location=x[2]))

    if(!is.null(build)){
      l$summarydata <- l$summarydata %>% lapply(., \(x) {x$build <- build; x})
    }

    return(l)
    })
  
  names(out) <- c(paste0("job", 1))
  out <- append(out, list(analyses=list(type = "LDScores")))

  
  # write yaml to file
  if (isTRUE(write)){yaml::write_yaml(out, outfile)}

  # return a yaml object
  out <- yaml::as.yaml(out)
  return(out)
}   





#' YAML constructor for the MR job
#' @param exposure A character vector of exposure traits to be analysed. Each trait should be prefixed with the source of the data ("opengwas" if querring the IEU OpenGWAS database, "vcf", "text", "robject" and "summaryset") followed by their location. For "opengwas", the location will be the IEU OpenGWAS ID, for "vcf" and "text" will be the name of the files, for the "robject" will be the name of a R object containing a dataset dataframe and "summaryset"  be the name of a `SummarySet()` object. Eg. "opengwas:ieu-a-2", "vcf:ieu-a-2.vcf.gz", "text:ieu-a-2.txt", "robject:dt_ieu-a-2", "summaryset:sumset_ieu-a-2". Note that files should be in the work session directory.
#' @param outcome A character vector of outcome traits to be analysed. The same rules as for the exposure apply.
#' @param build  A character vector  of the reference genome assemblies to generate the vcf file. Default is "GRCh37". 
#' * Options are `"NCBI34"`, `"NCBI35"`, `"NCBI36"`, `"GRCh37"` or "GRCh38".
#' @param pval A numeric vector of p-value thresholds to be used to select the top hits. Default is 5e-08.
#' @param clump A boolean indicating whether the top hits should be clumped. Default is TRUE.
#' @param r2 A numeric vector of r2 thresholds to be used to clump the top hits. Default is 0.001.
#' @param kb A numeric vector of kb thresholds to be used to clump the top hits. Default is 10000.
#' @param variant_col A character vector indicating the name of the column containing the variant id in the dataset. Default is "rsid".
#' @param pval_col A character vector indicating the name of the column containing the p-value in the dataset. Default is "p".
#' @param chr_col Column name for chromosome . Default is `"chr"`.
#' @param position_col Column name for the position. Together, with @param chr gives the physical coordinates of the variant. Default is `"position"`.
#' @param bfile It corresponds to the path and prefix of the plink files used when `clump = TRUE`. 
#' @param plink_bin Path to the plink executable
#' @param write A boolean indicating whether the yaml object should be written to a file.
#' @param outfile A character vector indicating the path to write the yaml file.
#' @return A yaml object.
#' @importFrom yaml write_yaml as.yaml
#' @export 

yaml_mr <- function(exposure, outcome, pop = "EUR", build = "GRCh37", pval = 5e-08, clump = TRUE, r2 = 0.001, kb = 10000, variant_col="rsid", pval_col = "p", chr_col = "chr", position_col = "position", bfile=NULL, plink_bin=NULL,write = FALSE, outfile = "config.yaml"){
  # read number of DataSets/jobs to  build and analyse
  njobs <- length(exposure)
  noutcomes <-length(outcome)

  out <- lapply(1:njobs, \(i){
    l <- list()
    # get instruments for each exposure
    if(grepl("^opengwas:", exposure[i])) {
       if (!requireNamespace("ieugwasr", quietly =TRUE)){
            stop("The MRC IEU R package `ieugwasr` needs to be installed.")}
      tophits <- ieugwasr::tophits(gsub("^opengwas:", "", exposure[i]), pop = pop)
      tophits <- paste0(tophits$chr, ":", tophits$pos)
    } 

    if(grepl("^vcf:", exposure[i])) {
       if (!requireNamespace("gwasvcf", quietly =TRUE)){
            stop("The MRC IEU R package `gwasvcf` needs to be installed.")}
      dat <- gwasvcf::readVcf(gsub("^vcf:", "", exposure[i])) %>% gwasvcf::vcf_to_tibble()
      tophits <- get_tophits_from_data(dat, pval = pval, clump = clump, r2 = r2, kb = kb, variant_col="ID", pval_col = "p", bfile=NULL, plink_bin=NULL)
      tophits <- paste0(tophits$seqnames, ":", tophits$start)
     }

    if(grepl("^text:", exposure[i])) {
      dat <- if (!requireNamespace("fread", quietly =TRUE)){
            stop("The CRAN package `fread` needs to be installed.")}
      data <- data.table::fread(gsub("^text:", "", exposure[i]), header = TRUE)
      
      tophits <- get_tophits_from_data(dat, pval = pval, clump = clump, r2 = r2, kb = kb, variant_col=variant_col, pval_col = pval_col, bfile=plink_bin, plink_bin=plink_bin)
      tophits <- paste0(tophits[chr_col], ":", tophits[position_col])
      }
    
    if(grepl("^robject:", exposure[i])) {
      dat <- get(gsub("^robject:", "", exposure[i]))
      if(!inherits(t,c("tbl", "tbl_df", "data.frame"))){
          stop("The object ", exposure[i], " is not a `data.frame`.`")
          }
      tophits <- get_tophits_from_data(dat, pval = pval, clump = clump, r2 = r2, kb = kb, variant_col=variant_col, pval_col = pval_col, bfile=bfile, plink_bin=plink_bin)
      tophits <- paste0(tophits[chr_col], ":", tophits[position_col])
      }
    if(grepl("^summaryset:", exposure[i])) {
      summaryset <- get(gsub("^summaryset:", "", exposure[i]))
      if(!isS4(summary_sets)){
         stop("The object ", exposure[i], " is not a SummarySet.")
        }
      tophits <- summaryset %>% getSummaryData() %>% get_tophits_from_data(pval = pval, clump = clump, r2 = r2, kb = kb, bfile=bfile, plink_bin=plink_bin)
      tophits <- paste0(tophits$chr, ":", tophits$position)
      }

    l$variants <- list(shape = "scattered", variant_list = tophits)
     l$summarydata <- strsplit(outcome, ":") %>% lapply(., \(x) list(source=x[1], location=x[2], organization = "outcome"))
     l$summarydata[1+noutcomes] <- strsplit(exposure[i], ":") %>% lapply(., \(x) list(source=x[1], location=x[2], organization = "exposure"))

   
    if(!is.null(build)){
      l$summarydata <- l$summarydata %>% lapply(., \(x) {x$build <- build; x})
    }
    
    return(l)
  })
    
  
  names(out) <- c(paste0("job", 1:njobs))
  out <- append(out, list(analyses=list(type = "MR")))

  
  # write yaml to file
  if (isTRUE(write)){yaml::write_yaml(out, outfile)}

  # return a yaml object
  out <- yaml::as.yaml(out)
  return(out)
}   

