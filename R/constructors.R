# NOTES: sample_size Sample size.  (first look at sumstats and look for max. If Na then some analyses wiil not run. Have a message in analyses for user add meta info)
# nsnp Number of variants in the study. Optional (length of the summary stats)_
# build   genome build version. Optional ? message

#' @title Metadata object
#' @description Reads metadata and converts it to gwasglue2 format.
#' 
#' @param metadata A dataframe with metadata information. Not required.
#' @param id GWAS study ID.
#' @param sample_size Sample size.
#' @param nsnp Number of variants in the study.
#' @param trait  Phenotype name corresponding the the variant.
#' @param sd Trait standard deviation.
#' @param unit TODO 
#' @param ncontrol TODO
#' @param build   genome build version.
#' @param population  Study sample population.
#' @param ncase Number of cases in study.
#' @param ... Other metadata information
#' @return A metadata list.
create_metadata <- function(metadata = NULL,
                           id = NA,
                           sample_size = NA,
                           nsnp = NA,
                           trait = NA,
                           sd = NA,
                           unit = NA,
                           ncontrol = NA, 
                           build = NA,
                           population = NA,
                           ncase = NA,
                           ...) {

  if (!is.null(metadata)){
    metadata <- as.list(c(metadata, ...))
    } else{
    metadata <- list (id = id , sample_size = sample_size, nsnp = nsnp, trait = trait, sd = sd, unit = unit, ncontrol = ncontrol, build = build, population = population, ncase = ncase, ...)
    }
  return(metadata)
}

#' A function to create a gwasglue2 SummarySet object from different sources and formats
#' 
#' @param data GWAS summary statistics.
#' @param metadata A list with metadata information. If NULL, it creates metadata with information retrieved from the dataset
#' @seealso [create_metadata()] to create a metadata object
#' @param type Input @param data type. Default is '"tibble"'. 
#' @param beta_col Name of column with effect sizes. The default is `"beta"`.
#' @param se_col Name of column with standard errors. The default is `"se"`.
#' @param eaf_col Name of column with effect allele frequency. The default is `"eaf"`.
#' @param effect_allele_col Name of column with effect allele. Must contain only the characters "A", "C", "T" or "G". The default is `"ea"`.
#' @param other_allele_col Name of column with non effect allele. Must contain only the characters "A", "C", "T" or "G". The default is `"nea`.
#' @param pvalue_col Name of column with p-value. The default is `"p"`.
#' @param samplesize_col Column name for sample size. The default is `"n"`.
#' @param chr_col Column name for chromosome . The default is `"chr"`.
#' @param position_col Column name for the position. Together, with @param chr gives the physical coordinates of the variant. The default is `"position"`.
#' @param rsid_col Required name of column with variants rs IDs. The default is `"rsid"`.
#' @param id_col GWAS study ID column. The default is `"id"`.
#' @param trait_col Column name for the column with phenotype name corresponding the the variant. The default is `"trait"` 
#' @param ... Other columns
#' @return A gwasglue2 SummarySet object
#' @export 
create_summaryset <- function (data,
                              metadata = NULL,
                              type = "tibble",
                              beta_col = "beta",
                              se_col = "se",
                              samplesize_col = "n",
                              pvalue_col = "p",
                              chr_col = "chr",
                              position_col = "position",
                              rsid_col = "rsid",
                              effect_allele_col = "ea",
                              other_allele_col = "nea",
                              eaf_col = "eaf",
                              id_col = "id",
                              trait_col = "trait",
                              ...) {
  if (type == "tibble") {
    s <- create_summaryset_from_tibble(data = data,
                              metadata = metadata,
                              beta_col =  beta_col,
                              se_col = se_col,
                              samplesize_col = samplesize_col,
                              pvalue_col = pvalue_col,
                              chr_col = chr_col,
                              position_col = position_col,
                              rsid_col = rsid_col,
                              effect_allele_col = effect_allele_col,
                              other_allele_col = other_allele_col,
                              eaf_col = eaf_col,
                              id_col = id_col,
                              trait_col = trait_col,
                              ...)
  }
  return(s)
}           


#' A function to create a gwasglue2 SummarySet object from a tibble
#' 
#' @param data GWAS summary statistics. A tibble
#' @param metadata A list with metadata information. If NULL, it creates metadata with information retrieved from the dataset
#' @seealso [create_metadata()] to create a metadata object
#' @param beta_col Name of column with effect sizes. The default is `"beta"`.
#' @param se_col Name of column with standard errors. The default is `"se"`.
#' @param eaf_col Name of column with effect allele frequency. The default is `"eaf"`.
#' @param effect_allele_col Name of column with effect allele. Must contain only the characters "A", "C", "T" or "G". The default is `"ea"`.
#' @param other_allele_col Name of column with non effect allele. Must contain only the characters "A", "C", "T" or "G". The default is `"nea`.
#' @param pvalue_col Name of column with p-value. The default is `"p"`.
#' @param samplesize_col Column name for sample size. The default is `"n"`.
#' @param chr_col Column name for chromosome . The default is `"chr"`.
#' @param position_col Column name for the position. Together, with @param chr gives the physical coordinates of the variant. The default is `"position"`.
#' @param rsid_col Required name of column with variants rs IDs. The default is `"rsid"`.
#' @param id_col GWAS study ID column. The default is `"id"`.
#' @param trait_col Column name for the column with phenotype name corresponding the the variant. The default is `"trait"` 
#' @param ... Other columns
#' @return A gwasglue2 SummarySet object
#' @export 
create_summaryset_from_tibble <- function (data = tibble(),
                              metadata = NULL,
                              beta_col = "beta",
                              se_col = "se",
                              samplesize_col = "n",
                              pvalue_col = "p",
                              chr_col = "chr",
                              position_col = "position",
                              rsid_col = "rsid",
                              effect_allele_col = "ea",
                              other_allele_col = "nea",
                              eaf_col = "eaf",
                              id_col = "id",
                              trait_col = "trait",
                              ...) {
  

  # Check column names and change them to ieugwasr nomenclature
  data_cols <- names(data)  
  ieu_req <- c("beta", "se", "p", "chr",  "position", "rsid", "ea",  "nea")
  ieu_optional <- c("n", "rsid", "eaf", "id", "trait")
  req_cols <- c(beta_col, se_col, pvalue_col, chr_col, position_col, rsid_col, effect_allele_col, other_allele_col) # nolint: line_length_linter

    if (all(req_cols%in%data_cols) == FALSE){
      dif <- setdiff(req_cols, data_cols)
      
      stop("Column(s) ", list(dif), " is/are missing. Please specify.",call. = FALSE)
    }

    if (all(req_cols%in%data_cols) == TRUE && all(ieu_req%in%req_cols) == FALSE){
 
      # replace column names by ieu_cols
      data1 <- data[req_cols]
      names(data1) <- ieu_req
      data2 <- data[,which(names(data) %ni% req_cols)]
      data <- dplyr::bind_cols(data1, data2)
          }
  
  #  standardise and sort summarySet by
     data <- data %>% 
            standardise(.) %>% 
            dplyr::arrange(chr, position, ea, nea)

  # create metadata if not in input from user
  if (is.null(metadata)){
    metadata = create_metadata()
  }
  # checks on metadata using data info
  if ("id" %in% colnames(metadata) && "id" %in% colnames(data) && is.na(metadata$id)){ 
    metadata$id <- unique(data$id)
    }
  if ("sample_size" %in% colnames(metadata) && "n" %in% colnames(data) && is.na(metadata$sample_size) && !all(is.na(data$n))){
    metadata$sample_size <- max(data$n, na.rm = TRUE)
    }
  if ("nsnp" %in% colnames(metadata) && is.na(metadata$nsnp)){ 
    metadata$nsnp <- dim(data)[1]
    }
  if ("trait" %in% colnames(metadata) && "trait" %in% colnames(data) && is.na(metadata$trait)){ 
    metadata$trait <- unique(data$trait)
    }

  s <- SummarySet(sumstats = data) %>%
    setMetadata(., metadata) %>%
    setVariantid(.) %>%
    setRSID(.,.@ss$rsid)
    return(s)
}


standardise <- function(d)
{
    toflip <- d$ea > d$nea
    d$eaf[toflip] <- 1 - d$eaf[toflip]
    d$beta[toflip] <- d$beta[toflip] * -1
    temp <- d$nea[toflip]
    d$nea[toflip] <- d$ea[toflip]
    d$ea[toflip] <- temp
    # d$snpid <- paste0(d$chr, ":", d$pos, "_", d$ea, "_", d$nea)
    d
}

#' Creates a DataSet object using gwasglue2 SummarySet objects, and harmonise data against data
#'
#' @param summary_sets A list of gwasglue2 SummarySet objects
#' @param harmonise logical (default TRUE). It harmonises the summary sets in the DataSet against each other. 
#' @param tolerance Inherited from harmoniseData() (default 0.08)
#' @param action Inherited from harmoniseData() (Default 1)
#' * `action = 1`: Assume all alleles are coded on the forward strand, i.e. do not attempt to flip alleles
#' * `action = 2`: Try to infer positive strand alleles, using allele frequencies for palindromes (default, conservative);
#' * `action = 3`: Correct strand for non-palindromic SNPs, and drop all palindromic SNPs from the analysis (more conservative).
#' @param strand dna strand orientation  (Default "forward", other option "reverse")
#' @return A harmonised gwasglue2 DataSet object
#' @export 
create_dataset <- function(summary_sets=list(),
                          harmonise = TRUE,
                          tolerance = 0.08,
                          action = 1,
                          strand = "forward") {

  ds <- DataSet(summary_sets) %>% overlapVariants(., strand = strand)
    
  if (harmonise == TRUE) {
     ds <-  harmoniseData(ds, tolerance = tolerance, action = action)
  }
  return(ds)
}


#' Creates a DataSet object using GWAS summary statistics, and harmonise data against data
#'
#' @param data A list of GWAS summary data (tibles)
#' @param metadata A list with metadata information. If NULL, it creates metadata with information retrieved from the dataset. 
#' @param harmonise logical (default TRUE). It harmonises the summary sets in the DataSet against each other. 
#' @param tolerance Inherited from harmoniseData() (default 0.08)
#' @param action Inherited from harmoniseData() (Default 1)
#' * `action = 1`: Assume all alleles are coded on the forward strand, i.e. do not attempt to flip alleles
#' * `action = 2`: Try to infer positive strand alleles, using allele frequencies for palindromes (default, conservative);
#' * `action = 3`: Correct strand for non-palindromic SNPs, and drop all palindromic SNPs from the analysis (more conservative).
#' @param strand dna strand orientation  (Default "forward", other option "reverse")
#' @param beta_col Name of column with effect sizes. The default is `"beta"`.
#' @param se_col Name of column with standard errors. The default is `"se"`.
#' @param eaf_col Name of column with effect allele frequency. The default is `"eaf"`.
#' @param effect_allele_col Name of column with effect allele. Must contain only the characters "A", "C", "T" or "G". The default is `"ea"`.
#' @param other_allele_col Name of column with non effect allele. Must contain only the characters "A", "C", "T" or "G". The default is `"nea`.
#' @param pvalue_col Name of column with p-value. The default is `"p"`.
#' @param samplesize_col Column name for sample size. The default is `"n"`.
#' @param chr_col Column name for chromosome . The default is `"chr"`.
#' @param position_col Column name for the position. Together, with @param chr gives the physical coordinates of the variant. The default is `"position"`.
#' @param rsid_col Required name of column with variants rs IDs. The default is `"rsid"`.
#' @param id_col The default is `"id"`.
#' @param trait_col Column name for the column with phenotype name corresponding the the variant. The default is `"trait"`
#' @param ... Other columns
#' @seealso [create_metadata()] to create a metadata object
#' @return A harmonised gwasglue2 DataSet object
#' @export
create_dataset_from_tibble <- function(data = list(),
                          metadata = NULL,
                          harmonise = TRUE,
                          tolerance = 0.08,
                          action = 1,
                          strand = "forward",
                          beta_col = "beta",
                          se_col = "se",
                          samplesize_col = "n",
                          pvalue_col = "p",
                          chr_col = "chr",
                          position_col = "position",
                          rsid_col = "rsid",
                          effect_allele_col = "ea",
                          other_allele_col = "nea",
                          eaf_col = "eaf",
                          id_col = "id",
                          trait_col = "trait",
                          ...) {

  ds <- DataSet(list())
  if (is.null(metadata)){
    metadata = vector("list", length(data))
  }
  for (i in seq_along(data)){
    s <- create_summaryset(data= data[[i]],
        metadata = metadata[[i]],
        source = source,
        beta_col = beta_col,
        se_col = se_col,
        samplesize_col = samplesize_col,
        pvalue_col = pvalue_col,
        chr_col = chr_col,
        position_col = position_col,
        rsid_col = rsid_col,
        effect_allele_col = effect_allele_col,
        other_allele_col = other_allele_col,
        eaf_col = eaf_col,
        id_col = id_col,
        trait_col = trait_col)
    
    ds@summary_sets[[i]] <- s
  }

  ds <- ds %>% overlapVariants(., strand = strand)
    
  if (harmonise == TRUE) {
     ds <-  harmoniseData(ds, tolerance = tolerance, action = action)
  }

  return(ds)
}


#' Add a SummarySet to a DataSet 
#'
#' @param summary_sets one or more gwasglue2 Summarysets objects to add to an existent DataSet object. If more than one it should be a list
#' @param dataset The gwasglue2 DataSet object to add to
##' @param harmonise logical (default TRUE). It harmonises the summary sets in the DataSet against each other. 
#' @param tolerance Inherited from harmoniseData() (default 0.08)
#' @param action Inherited from harmoniseData() (Default 1)
#' * `action = 1`: Assume all alleles are coded on the forward strand, i.e. do not attempt to flip alleles
#' * `action = 2`: Try to infer positive strand alleles, using allele frequencies for palindromes (default, conservative);
#' * `action = 3`: Correct strand for non-palindromic SNPs, and drop all palindromic SNPs from the analysis (more conservative).
#' @param strand dna strand orientation  (Default "forward", other option "reverse")
#' @return A harmonised gwasglue2 DataSet object with input SummarySets added
#' @export 
add_summaryset <- function(summary_sets,
                          dataset,
                          harmonise = TRUE,
                          tolerance = 0.08,
                          action = 1,
                          strand = "forward") {
  
  l <- length(dataset@summary_sets)
  
  if(length(summary_sets) == 0){
     dataset@summary_sets[[l + 1]] <- summary_sets
  }else {
    for(i in seq_along(summary_sets)){
    dataset@summary_sets[[l + i]] <- summary_sets[[i]]
    }
  }

  message("\nSummarySet added to DataSet")
  
  ds <- dataset %>% overlapVariants(., strand = strand)
    
  if (harmonise == TRUE) {
     ds <-  harmoniseData(ds, tolerance = tolerance, action = action)
  }
  return(ds)

} 


#' Merge Datasets 
#'
#' @param datasets A list of gwasglue2 DataSet objects
#' @return  A gwasglue2 DataSet object  with input DataSets merged
merge_datasets <- function(datasets) {
    n_datasets <- length(datasets)
	  ds <- DataSet(list())
    count <- 1
    for(i in 1:n_datasets){
      if (isS4(datasets[[i]]) == FALSE){
         next
         }
       n_ss <- length(datasets[[i]]@summary_sets)
       ds@ld_matrix <- datasets[[i]]@ld_matrix
       for(j in 1:n_ss){
            ds@summary_sets[[count]] <- datasets[[i]]@summary_sets[[j]]
            count <- count + 1
            }     
    }
	# For now I am assuming that in all datasets the following slots are the same
    # ds@overlap_variants <- datasets[[1]]@overlap_variants
    # ds@is_resized <- datasets[[1]]@is_resized
    # ds@is_harmonised <- datasets[[1]]@is_harmonised
    # ds@overall_dropped_SNPs <- datasets[[1]]@overall_dropped_SNPs
    # ds@dropped_SNPs <- datasets[[1]]@dropped_SNPs
    # ds@palindromic_SNPs <- datasets[[1]]@palindromic_SNPs
    # ds@ambiguous_SNPs <- datasets[[1]]@ambiguous_SNPs
    # ds@incompatible_alleles_SNPs <- datasets[[1]]@incompatible_alleles_SNPs
    
    # ds@is_harmonisedLD <- datasets[[1]]@is_harmonisedLD
    # ds@zscores <- datasets[[1]]@zscores
    # ds@is_converted <- datasets[[1]]@is_converted

    return(ds)
}

