# ieugwasr Functions

#' A function to clump the top hit variants
#'
#' @param traits ID of GWAS studies to query
#' @param clump Logical (default TRUE)
#' @param source  IEU (default)
#'
#' @return Array of SNPs rsids
clumpTophits <- function(traits,  clump = TRUE, source = ieugwasr::check_access_token()) {
  ch  <- ieugwasr::tophits(traits[1], clump = clump)
  snps <- ch$rsid
  return(snps)
}

#' PHEWAS_ids
#' @param variants Array of variants
#' @param pval p-value threshold. Default = `0.00001`. Inherited from ieugwasr::phewas # nolint: line_length_linter.
#' @param batch Vector of batch IDs to search across. If `c()` (default) then returns all batches. Iherited fromieugwasr::phewas  # nolint: line_length_linter.
#' @return id_list Array with traits ids
phewasIDs <- function(variants, batch, pval) {
  id_list <- unique(ieugwasr::phewas(pval = pval, variants = variants, batch = batch)$id) 
  return(id_list)
}

#' A function to call the summary statistics associated
#' with specific traits and variants chosen
#'
#' @param traits ID of GWAS studies to query
#' @param variants Array of SNPs rsids
#'
#' @return A tibble with the summary statistics for variant association.
createSumset <- function(traits,variants) {
  x <- ieugwasr::associations(variants = variants, id = traits)
  x <- dplyr::arrange(x,rsid)
  return(x)
}

# NOTES: sample_size Sample size.  (first look at sumstats and look for max. If Na then some analyses wiil not run. Have a message in analyses for user add meta info)
# nsnp Number of variants in the study. Optional (length of the summary stats)_
# build   genome build version. Optional ? message

#' Reads metadata and converts it to gwasglue2 format
#' 
#' @param metadata A dataframe with metadata information.
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

            


#' A function to create a gwasglue2 SummarySet object
#' 
#' @param data GWAS summary statistics. A tibble
#' @param metadata A list with metadata information. If NULL, it creates metadata with information retrieved from the dataset
#' @seealso [create_metadata()] to create a metadata object
#' @param tools Array of methods that gwasglue2 is going to convert the summarySet to (eg. "mr") 
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


create_summaryset <- function (data = tibble(),
                              metadata = NULL,
                              tools,
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
  if (is.na(metadata$id)){ 
    metadata$id <- unique(data$id)
    }
  if (is.na(metadata$sample_size) && !all(is.na(data$n))){
    metadata$sample_size <- max(data$n, na.rm = TRUE)
    }
  if (is.na(metadata$nsnp)){ 
    metadata$nsnp <- dim(data)[1]
    }
  if (is.na(metadata$trait)){ 
    metadata$trait <- unique(data$trait)
    }

  s <- SummarySet(sumstats = data) %>%
    setMetadata(., metadata) %>%
    setVariantid(.) %>%
    setRSID(.,.@ss$rsid) %>%
    setTool(., tools = tools)

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



#' Creates a DataSet object and harmonise data against data
#'
#' @param data A list of gwasglue2 SummarySet objects
#' @param metadata A list with metadata information. If NULL, it creates metadata with information retrieved from the dataset. 
#' @param tools Array of methods that gwasglue2 is going to convert the summarySet to (eg. "mr") 
#' @param harmonise logical (default TRUE). It harmonises the summary sets in the DataSet against each other. 
#' @param tolerance Inherited from harmoniseData() (default 0.08)
#' @param action Inherited from harmoniseData() (Default 1)
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
#' @seealso [create_metadata()] to create a metadata object
#' @return A harmonised gwasglue2 DataSet object
#' @export
create_dataset <- function(data = list(),
                          metadata = list(NULL),
                          tools = NULL,
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
                          trait_col = "trait") {

  ds <- DataSet()
  for (i in seq_along(data)){
    s <- create_summaryset(data= data[[i]],
        metadata = metadata[[i]],
        tools = tools,
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



#' Plot
#'
#' @param dataset gwasglue2 DataSet object
#' @param type Type of plot (Only available "manhattan" plots at the moment)
#' @param title Main title for the plot
#' @export
#'
#' @return A plot
plot_gwasglue <- function(dataset, type, title){
  
  if(type == "manhattan"){
    ntraits <- getLength(dataset)
    nb_rows <- ceiling(ntraits/2)
     # Add main title
    
    
    par(mfrow=c(nb_rows, 2))
    
    
      for (i in 1:ntraits){
        plot(dataset@summary_sets[[i]]@ss$position, -log10(dataset@summary_sets[[i]]@ss$p), main = "", xlab = "position", ylab = "-log10(p-value)", cex=0.8, pch=20)
        mtext(dataset@summary_sets[[i]]@metadata$trait, side = 3, line = 0.5)
      }
    mtext(as.expression(bquote(bold(.(title)))),
          side = 3,
          line = - 2.5,
          outer = TRUE,
          cex = 1.3)
  }
}


