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
#' @param unit Unit 
#' @param ncontrol Nb of controls in study
#' @param build   genome build version.
#' @param population  Study sample population.
#' @param ncase Number of cases in study.
#' @param ... Other metadata information
#' @return A metadata list.
#' @export
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
#' @param type Input @param data type. Default is `"tibble"`. Other options: `"vcf"`
#' @param build Reference genome assembly to generate the genomic data. Default is NULL. 
#' * Options are `"NCBI34"`, `"NCBI35"`, `"NCBI36"`, `"GRCh37"` or "GRCh38".
#' @param qc Quality control. It checks the @param data and look for problems that can stop gwasglue2 from running. If TRUE gwasglue2 will try to solve the problems.  Default is FALSE
#' @param beta_col Name of column with effect sizes. The default is `"beta"` for @param type `"tibble"` and `"ES"`for @param type `"vcf"`..
#' @param se_col Name of column with standard errors. The default is `"se"` for @param type `"tibble"` and `"SE"`for @param type `"vcf"`.
#' @param eaf_col Name of column with effect allele frequency. The default is `"eaf"` for @param type `"tibble"` and `"AF"`for @param type `"vcf"`.
#' @param effect_allele_col Name of column with effect allele. Must contain only the characters "A", "C", "T" or "G". The default is `"ea"` for @param type `"tibble"` and `"ALT"`for @param type `"vcf"`.
#' @param other_allele_col Name of column with non effect allele. Must contain only the characters "A", "C", "T" or "G". The default is `"nea` for @param type `"tibble"` and `"REF"`for @param type `"vcf"`.
#' @param pvalue_col Name of column with p-value. The default is `"p"`.
#' @param logpvalue_col Name of column with log(p-value). The default is `"LP"` for @param type `"vcf"`.
#' @param samplesize_col Column name for sample size. The default is `"n"` for @param type `"tibble"` and `"SS"`for @param type `"vcf"`.
#' @param chr_col Column name for chromosome . The default is `"chr"` for @param type `"tibble"` and `"seqnames"`for @param type `"vcf"`.
#' @param position_col Column name for the position. Together, with @param chr gives the physical coordinates of the variant. The default is `"position"` for @param type `"tibble"` and `"start"`for @param type `"vcf"`.
#' @param rsid_col Required name of column with variants rs IDs. The default is `"rsid"` for @param type `"tibble"` and `"ID"`for @param type `"vcf"`.
#' @param id_col GWAS study ID column. The default is `"id"`.
#' @param trait_col Column name for the column with phenotype name corresponding the the variant. The default is `"trait"` 
#' @return A gwasglue2 SummarySet object
#' @export 
create_summaryset <- function (data,
                              metadata = NULL,
                              type = "tibble",
                              qc = FALSE,
                              beta_col = NULL,
                              se_col =NULL,
                              samplesize_col = NULL,
                              pvalue_col = NULL,
                              logpvalue_col = NULL,
                              chr_col = NULL,
                              position_col = NULL,
                              rsid_col = NULL,
                              effect_allele_col = NULL,
                              other_allele_col = NULL,
                              eaf_col = NULL,
                              id_col = NULL,
                              trait_col = NULL, 
                              build = NULL) {
  # create SummarySet from tibbles
  if (type == "tibble") {
    if(is.null(beta_col)){
      beta_col = "beta"
    }
    if(is.null(se_col)){
      se_col = "se"
    }
    if(is.null(samplesize_col)){
      samplesize_col = "n"
    }
    if(is.null(pvalue_col)){
      pvalue_col = "p"
    }
    if(is.null(chr_col)){
       chr_col = "chr"
    }
    if(is.null(position_col)){
      position_col = "position"
    }
    if(is.null(rsid_col)){
      rsid_col = "rsid"
    }
    if(is.null(effect_allele_col)){
      effect_allele_col = "ea"
    }
    if(is.null(other_allele_col)){
      other_allele_col = "nea"
    }
    if(is.null(eaf_col)){
      eaf_col = "eaf"
    }
    if(is.null(id_col)){
      id_col = "id"
    }
    if(is.null(trait_col)){
      trait_col = "trait"
    }

    s <- create_summaryset_from_tibble(data = data,
                              metadata = metadata,
                              qc = qc,
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
                              build = build)
  }
  # create SummarySet from vcf
  if (type == "vcf") {
    if(is.null(beta_col)){
      beta_col = "ES"
    }
    if(is.null(se_col)){
      se_col = "SE"
    }
    if(is.null(samplesize_col)){
     samplesize_col = "SS"
    }
    if(is.null(logpvalue_col)){
     logpvalue_col = "LP"
    }
    if(is.null(pvalue_col)){
      pvalue_col = "p"
    }
    if(is.null(chr_col)){
      chr_col = "seqnames"
    }
    if(is.null(position_col)){
       position_col = "start"
    }
    if(is.null(rsid_col)){
      rsid_col = "ID"
    }
    if(is.null(effect_allele_col)){
      effect_allele_col = "ALT"
    }
    if(is.null(other_allele_col)){
      other_allele_col = "REF"
    }
    if(is.null(eaf_col)){
      eaf_col = "AF"
    }
    if(is.null(id_col)){
      id_col = "id"
    }
    
    s <- create_summaryset_from_gwasvcf(data = data,
                              metadata = metadata,
                              qc = qc,
                              beta_col =  beta_col,
                              se_col = se_col,
                              samplesize_col = samplesize_col,
                              pvalue_col = pvalue_col,
                              logpvalue_col = logpvalue_col,
                              chr_col = chr_col,
                              position_col = position_col,
                              rsid_col = rsid_col,
                              effect_allele_col = effect_allele_col,
                              other_allele_col = other_allele_col,
                              eaf_col = eaf_col,
                              id_col = id_col,
                              build = build)
  }

  return(s)
}           


#' A function to create a gwasglue2 SummarySet object from a tibble or dataframe
#' 
#' @param data GWAS summary statistics. A tibble or a dataframe
#' @param metadata A list with metadata information. If NULL, it creates metadata with information retrieved from the dataset
#' @seealso [create_metadata()] to create a metadata object
#' @param build Reference genome assembly to generate the genomic data. Default is NULL. 
#' * Options are `"NCBI34"`, `"NCBI35"`, `"NCBI36"`, `"GRCh37"` or "GRCh38".
#' @param qc Quality control. It checks the @param data and look for problems that can stop gwasglue2 from runing. If TRUE gwasglue will try to solve the problems.  Default is FALSE
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
#' @return A gwasglue2 SummarySet object
#' @export 
create_summaryset_from_tibble <- function (data,
                              metadata = NULL,
                              qc = FALSE,
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
                              build = NULL) {
  
  data <- data %>% dplyr::as_tibble()                              
  # Check column names and change them to ieugwasr nomenclature
  data_cols <- names(data)  

  # IEU required columns
  ieu_req <- c("beta", "se", "p", "chr",  "position", "rsid", "ea",  "nea")
  req_cols <- c(beta_col, se_col, pvalue_col, chr_col, position_col, rsid_col, effect_allele_col, other_allele_col) # nolint: line_length_linter
  # IEU optional columns
  ieu_optional <- c("n", "eaf", "id","trait")
  optional_cols <- c(samplesize_col, eaf_col, id_col,trait_col)

# stop if any required column is missing
  if (isFALSE(all(req_cols%in%data_cols))){
    dif <- setdiff(req_cols, data_cols)
    stop("Column(s) ", list(dif), " is/are missing. Please specify.",call. = FALSE)
  }

  # replace column names by optional ieu_cols
  if (isTRUE(any(optional_cols %in% data_cols))){
    data1 <- data[optional_cols[optional_cols %in% data_cols]]
    names(data1) <- ieu_optional[optional_cols %in% data_cols]
    data2 <- data[,which(names(data) %ni% ieu_optional & names(data) %ni% optional_cols)]
    data <- dplyr::bind_cols(data1, data2)
  }

     # replace column names by required ieu_cols
  if (isTRUE(all(req_cols%in%data_cols))){
    data1 <- data[req_cols]
    names(data1) <- ieu_req
      
    data2 <- data[,which(names(data) %ni% ieu_req & names(data) %ni% req_cols)]
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
  if ("id" %in% names(metadata) && "id" %in% colnames(data) && is.na(metadata$id)){ 
    metadata$id <- unique(data$id)
    }
  if ("sample_size" %in% names(metadata) && "n" %in% colnames(data) && is.na(metadata$sample_size) && !all(is.na(data$n))){
    metadata$sample_size <- max(data$n, na.rm = TRUE)
    }
  if ("trait" %in% names(metadata) && "trait" %in% colnames(data) && is.na(metadata$trait)){ 
    metadata$trait <- unique(data$trait)
    }
  if ("build" %in% names(metadata) && is.na(metadata$build)){ 
    metadata$build <- build
    }

  s <- SummarySet(sumstats = data) %>%
    setMetadata(., metadata) %>%
    setVariantid(.) %>%
    setRSID(.,.@ss$rsid) 

  # sanity checks
  # check if there are repeated variantids
  variant_id <- s@ss$variantid
  variant_id_rep <- variant_id[(duplicated(variant_id) | duplicated(variant_id, fromLast = TRUE)) ]
  
  # warning message if there is problems run with option
  if (length(variant_id_rep) > 0){
    if (isFALSE(qc)){
    warning( "There are repeated variants in the data, that could stop gwasglue2 from running in some analyses. If you want gwasglue to deal with it, run again 'create_summaryset' with option 'qc = TRUE'\n" )
    print(s@ss[which(s@ss$variantid %in% variant_id_rep), ] %>% dplyr::select(., variantid, chr, position, ea, nea, eaf,beta, se))
    }
    # remove repeated variantids
    if (isTRUE(qc)){
      message("gwasglue2 is removing problematic variants\n")
      print(s@ss[which(s@ss$variantid %in% variant_id_rep), ] %>% dplyr::select(., variantid, chr, position, ea, nea, eaf,beta, se))
      s@ss <- s@ss[which(s@ss$variantid %ni% variant_id_rep), ]
    }
  }


  # set attributes
  s@attributes <- list("type" = "tibble", "creation" = Sys.time())
  
  return(s)
}

#' A function to create a gwasglue2 SummarySet object from a vcf file
#' 
#' @param data GWAS summary statistics. In the GWAS vcf dataframe format
#' @param metadata A list with metadata information. If NULL, it creates metadata with information retrieved from the dataset
#' @seealso [create_metadata()] to create a metadata object
#' @param build Reference genome assembly to generate the vcf file. Default is NULL. 
#' * Options are `"NCBI34"`, `"NCBI35"`, `"NCBI36"`, `"GRCh37"` or "GRCh38".
#' @param qc Quality control. It checks the @param data and look for problems that can stop gwasglue2 from runing. If TRUE gwasglue will try to solve the problems.  Default is FALSE
#' @param beta_col Name of column with effect sizes. The default is `"ES"`.
#' @param se_col Name of column with standard errors. The default is `"SE"`.
#' @param eaf_col Name of column with effect allele frequency. The default is `"AF"`.
#' @param effect_allele_col Name of column with effect allele. Must contain only the characters "A", "C", "T" or "G". The default is `"ALT"`.
#' @param other_allele_col Name of column with non effect allele. Must contain only the characters "A", "C", "T" or "G". The default is `"REF"`.
#' @param pvalue_col Name of column with p-value. The default is `"p"`.
#' @param logpvalue_col Name of column with log(p-value). The default is `"LP"`.
#' @param samplesize_col Column name for sample size. The default is `"SS"`.
#' @param chr_col Column name for chromosome . The default is `"seqnames"`.
#' @param position_col Column name for the position. Together, with @param chr gives the physical coordinates of the variant. The default is `"start"`.
#' @param rsid_col Required name of column with variants rs IDs. The default is `"ID"`.
#' @param id_col GWAS study ID column. The default is `"id"`.
#' @return A gwasglue2 SummarySet object
#' @export 
create_summaryset_from_gwasvcf <- function (data,
                              metadata = NULL,
                              qc = FALSE,
                              beta_col = "ES",
                              se_col = "SE",
                              samplesize_col = "SS",
                              logpvalue_col = "LP",
                              pvalue_col = "p",
                              chr_col = "seqnames",
                              position_col = "start",
                              rsid_col = "ID",
                              effect_allele_col = "ALT",
                              other_allele_col = "REF",
                              eaf_col = "AF",
                              id_col = "id", 
                              build = NULL) {

  data <- data %>% dplyr::as_tibble()

  # Check column names and change them to ieugwasr nomenclature
  data_cols <- names(data)  

  # transform the pvalue if it is in log10
  if(isTRUE(any(data_cols==logpvalue_col)) & isFALSE(any(data_cols==pvalue_col))){
    LP <- data[which(colnames(data)==logpvalue_col)]
    data$p <- 10^-data$LP
    data_cols <- names(data)
  }

  # IEU required columns
  ieu_req <- c("beta", "se", "p", "chr",  "position", "rsid", "ea",  "nea")
  req_cols <- c(beta_col, se_col, pvalue_col, chr_col, position_col, rsid_col, effect_allele_col, other_allele_col) # nolint: line_length_linter
  # IEU optional columns
  # TODO add trait?
  ieu_optional <- c("n", "eaf", "id")
  optional_cols <- c(samplesize_col, eaf_col, id_col)

# stop if any required column is missing
  if (isFALSE(all(req_cols%in%data_cols))){
    dif <- setdiff(req_cols, data_cols)
    stop("Column(s) ", list(dif), " is/are missing. Please specify.",call. = FALSE)
  }

  # replace column names by optional ieu_cols
  if (isTRUE(any(optional_cols%in%data_cols))){
    data1 <- data[optional_cols[optional_cols%in%data_cols]]
    names(data1) <- ieu_optional[optional_cols%in%data_cols]
    data2 <- data[,which(names(data) %ni% ieu_optional & names(data) %ni% optional_cols)]
    data <- dplyr::bind_cols(data1, data2)
  }

     # replace column names by required ieu_cols
  if (isTRUE(all(req_cols%in%data_cols))){
    data1 <- data[req_cols]
    names(data1) <- ieu_req
      
    data2 <- data[,which(names(data) %ni% ieu_req & names(data) %ni% req_cols)]
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
  if ("id" %in% names(metadata) && "id" %in% colnames(data) && is.na(metadata$id)){ 
    metadata$id <- unique(data$id)
    }
  if ("sample_size" %in% names(metadata) && "n" %in% colnames(data) && is.na(metadata$sample_size) && !all(is.na(data$n))){
    metadata$sample_size <- max(data$n, na.rm = TRUE)
    }
  if ("trait" %in% names(metadata) && "trait" %in% colnames(data) && is.na(metadata$trait)){ 
    metadata$trait <- unique(data$trait)
    }
  if ("build" %in% names(metadata) && is.na(metadata$build)){ 
    metadata$build <- build
    }

  s <- SummarySet(sumstats = data) %>%
    setMetadata(., metadata) %>%
    setVariantid(.) %>%
    setRSID(.,.@ss$rsid) 

  # sanity checks
  # check if there are repeated variantids
  variant_id <- s@ss$variantid
  variant_id_rep <- variant_id[(duplicated(variant_id) | duplicated(variant_id, fromLast = TRUE)) ]
  
  # warning message if there is problems run with option
  if (length(variant_id_rep) > 0){
    if (isFALSE(qc)){
    warning( "There are repeated variants in the data, that could stop gwasglue2 from running in some analyses. If you want gwasglue to deal with it, run again 'create_summaryset' with option 'qc = TRUE'\n")
    print(s@ss[which(s@ss$variantid %in% variant_id_rep), ] %>% dplyr::select(., variantid, chr, position, ea, nea, eaf,beta, se))
    }
    # remove repeated variantids
    if (isTRUE(qc)){
      message("gwasglue2 is removing problematic variants\n")
      print(s@ss[which(s@ss$variantid %in% variant_id_rep), ] %>% dplyr::select(., variantid, chr, position, ea, nea, eaf,beta, se))
      s@ss <- s@ss[which(s@ss$variantid %ni% variant_id_rep), ]
    }
  }
    

  # set attributes
  s@attributes <- list("type" = "vcf", "creation" = Sys.time())
  
  return(s)
}





#' Creates a DataSet object using gwasglue2 SummarySet objects, and harmonise data against data
#'
#' @param summary_sets A list of gwasglue2 SummarySet objects
#' @param harmonise logical (default TRUE). It harmonises the summary sets in the DataSet against each other. 
#' @param tolerance Inherited from harmoniseData() (default 0.08)
#' @param action Inherited from harmoniseData() (Default 1)
#' * `action = 1`: Assume all alleles are coded on the forward strand, i.e. do not attempt to flip alleles (default);
#' * `action = 2`: Try to infer positive strand alleles, using allele frequencies for palindromes ;
#' * `action = 3`: Correct strand for non-palindromic SNPs, and drop all palindromic SNPs from the analysis (more conservative).
#' @return A harmonised gwasglue2 DataSet object
#' @export 
create_dataset <- function(summary_sets=list(),
                          harmonise = TRUE,
                          tolerance = 0.08,
                          action = 1) {
  ds <- DataSet(summary_sets)

  if (harmonise == TRUE) {
    ds <- ds %>% overlapVariants(., action = action)
  }
  if (action != 1 && harmonise == TRUE){
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
    s <- create_summaryset_from_tibble(data= data[[i]],
        metadata = metadata[[i]],
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

  if (harmonise == TRUE) {
    ds <- ds %>% overlapVariants(., action = action)
  }
  if (action != 1 && harmonise == TRUE){
        ds <-  harmoniseData(ds, tolerance = tolerance, action = action)
  } 
  return(ds)
}


#' Add a SummarySet to a DataSet 
#'
#' @param summary_sets one or more gwasglue2 Summarysets objects to add to an existent DataSet object. If more than one it should be a list
#' @param dataset The gwasglue2 DataSet object to add to
#' @param harmonise logical (default TRUE). It harmonises the summary sets in the DataSet against each other. 
#' @param tolerance Inherited from harmoniseData() (default 0.08)
#' @param action Inherited from harmoniseData() (Default 1)
#' * `action = 1`: Assume all alleles are coded on the forward strand, i.e. do not attempt to flip alleles
#' * `action = 2`: Try to infer positive strand alleles, using allele frequencies for palindromes (default, conservative);
#' * `action = 3`: Correct strand for non-palindromic SNPs, and drop all palindromic SNPs from the analysis (more conservative).
#' @return A harmonised gwasglue2 DataSet object with input SummarySets added
#' @export 
add_summaryset <- function(summary_sets,
                          dataset,
                          harmonise = TRUE,
                          tolerance = 0.08,
                          action = 1) {
  
  l <- length(dataset@summary_sets)
  
  if(length(summary_sets) == 0){
     dataset@summary_sets[[l + 1]] <- summary_sets
  }else {
    for(i in seq_along(summary_sets)){
    dataset@summary_sets[[l + i]] <- summary_sets[[i]]
    }
  }

  message("\nSummarySet added to DataSet")
  
  ds <- dataset

  if (harmonise == TRUE) {
    ds <- ds %>% overlapVariants(., action = action)
  }
  if (action != 1 && harmonise == TRUE){
        ds <-  harmoniseData(ds, tolerance = tolerance, action = action)
  } 
  return(ds)

} 


#' Merge Datasets 
#'
#' @param datasets A list of gwasglue2 DataSet objects
#' @param harmonise logical (default TRUE). It harmonises the summary sets in the DataSet against each other. 
#' @param tolerance Inherited from harmoniseData() (default 0.08)
#' @param action Inherited from harmoniseData() (Default 1)
#' * `action = 1`: Assume all alleles are coded on the forward strand, i.e. do not attempt to flip alleles
#' * `action = 2`: Try to infer positive strand alleles, using allele frequencies for palindromes (default, conservative);
#' * `action = 3`: Correct strand for non-palindromic SNPs, and drop all palindromic SNPs from the analysis (more conservative).
#' @return  A gwasglue2 DataSet object  with input DataSets merged
#' @export
merge_datasets <- function(datasets,
                          harmonise = TRUE,
                          tolerance = 0.08,
                          action = 1) {
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
    # harmonise
    if (harmonise == TRUE) {
    ds <- ds %>% overlapVariants(., action = action)
    }
    if (action != 1 && harmonise == TRUE){
          ds <-  harmoniseData(ds, tolerance = tolerance, action = action)
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






# standardise function: it standardises the data and flips the alleles to be in alphabetical order. Thus, the effect allele will always be the first alphabetically. 
standardise <- function(d){
    toflip <- d$ea > d$nea
    d$eaf[toflip] <- 1 - d$eaf[toflip]
    d$beta[toflip] <- d$beta[toflip] * -1
    temp <- d$nea[toflip]
    d$nea[toflip] <- d$ea[toflip]
    d$ea[toflip] <- temp
    return(d)
}


#'gwasglue2 variant IDs system
#'@export
create_variantid <-function(chr,pos,a1,a2) {
  if (!requireNamespace("digest", quietly =TRUE)){
    stop("The CRAN package `digest` needs to be installed.")}
  alleles_sorted <- t(apply(cbind(a1,a2),1,sort)) 
  #  create variantid
  variantid <- paste0(chr,":", pos,"_",alleles_sorted[,1],"_",alleles_sorted[,2])

  # create hashes when alleles nchar > 10
  # allele ea
   if (all(nchar(alleles_sorted[,1]) <= 10) == FALSE){
    index = which(nchar(alleles_sorted[,1]) > 10)
    variantid[index] <- lapply(index, function(i){
      v <- paste0(chr[i],":", pos[i],"_#",digest::digest(alleles_sorted[i,1],algo= "murmur32"),"_",alleles_sorted[i,2],) 
    }) %>% unlist()
  }

  # allele nea
  if (all(nchar(alleles_sorted[,2]) <= 10) == FALSE){
    index = which(nchar(alleles_sorted[,2]) > 10)
    variantid[index] <- lapply(index, function(i){
      v <- paste0(chr[i],":", pos[i],"_",alleles_sorted[i,1],"_#",digest::digest(alleles_sorted[i,2],algo= "murmur32")) 
    }) %>% unlist()
  }
  
  return(variantid)
  }

