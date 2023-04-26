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
#' @param pval p-value threshold. Default = `0.00001`. Iherited from ieugwasr::phewas
#' @param batch Vector of batch IDs to search across. If `c()` (default) then returns all batches. Iherited fromieugwasr::phewas
#' @return id_list Array with traits ids
phewasIDs <- function(variants, batch, pval) {
  id_list <- unique(ieugwasr::phewas(pval = pval, variants = variants, batch = batch)$id) 
  return(id_list)
}

#' A function to call the summary statistics associated with specific traits and variants chosen
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


#' A function to create a gwasglue2 SummarySet object
#' 
#' @param data GWAS summary statistics. A tibble
#' @param id  GWAS study ID
#' @param tools Array of methods that gwasglue2 is going to convert the summarySet to (eg. "mr") 
#' @param source ID of Data Source to set the DataSet @slot metadata (Default and only option: IEUopenGWAS). 
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
#' @return A gwasglue2 SummarySet object
#' @export 


constructSummarySet <- function (data = tibble(), id, tools, source = "IEUopenGWAS", beta_col = "beta", se_col = "se", samplesize_col = "n", pvalue_col = "p",
chr_col = "chr", position_col = "position", rsid_col = "rsid", effect_allele_col = "ea", other_allele_col = "nea", eaf_col = "eaf", id_col = "id", trait_col = "trait"){

  # Check column names and change them to ieugwasr nomenclature 
  data_cols <- names(data)  
  ieu_cols <- c("beta", "se", "n",  "p", "chr",  "position", "rsid", "ea",  "nea",  "eaf", "id", "trait")
  cols = c(beta_col, se_col, samplesize_col, pvalue_col, chr_col, position_col, rsid_col, effect_allele_col, other_allele_col, eaf_col, id_col, trait_col)
    
    if (all(cols%in%data_cols) == FALSE){
      dif <- setdiff(cols, data_cols)
      
      stop("The columns ", list(dif), " are not present. Please specify them.",call. = FALSE)
    }

    if (all(cols%in%data_cols) == TRUE && all(ieu_cols%in%cols) == FALSE){
      # order data by cols variable 
      # TODO: what if data has more columns?
      data <- data[cols]
      # replace column names by ieu_cols (same order as cols)
      names(data) <- ieu_cols
    }

  s <- SummarySet(sumstats = data) %>%
    setMetadata(., source = source, id = id) %>%
    setRSID(.,.@ss$rsid) %>%
    setTool(., tools = tools)

    return(s)
   }


#' Creates a DataSet object and harmonise data against data
#'
#' @param data A list of gwasglue2 SummarySet objects
#' @param ids ID of GWAS studies to query
#' @param tools Array of methods that gwasglue2 is going to convert the summarySet to (eg. "mr") 
#' @param source ID of Data Source to set the DataSet @slot metadata (Default and only option: IEUopenGWAS). 
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
#' @export
#'
#' @return A harmonised gwasglue2 DataSet object

createDataSet <- function(data=list(), ids, tools, source="IEUopenGWAS", harmonise = TRUE, tolerance = 0.08,action = 1, strand = "forward", beta_col = "beta", se_col = "se", samplesize_col = "n", pvalue_col = "p", chr_col = "chr", position_col = "position", rsid_col = "rsid", effect_allele_col = "ea", other_allele_col = "nea", eaf_col = "eaf", id_col = "id", trait_col = "trait"){

  ds <- DataSet()
  for (i in seq_along(data)){
    s <- constructSummarySet(data= data[[i]],
                            id = ids[i],
                            tools = tools,
                            source = source,
                            beta_col = beta_col, se_col = se_col, samplesize_col = samplesize_col, pvalue_col = pvalue_col, chr_col = chr_col, position_col = position_col, rsid_col = rsid_col, effect_allele_col = effect_allele_col, other_allele_col = other_allele_col, eaf_col = eaf_col, id_col = id_col, trait_col = trait_col)
  ds@summary_sets[[i]] <- s
  }
  ds <- overlapVariants(ds, strand = strand)
    
  if (harmonise == TRUE && strand == "forward"){
     ds <-  harmoniseData(ds,tolerance = tolerance, action = action, strand = "forward")
  }

  if (harmonise == TRUE && strand == "reverse"){
     ds <-  harmoniseData(ds,tolerance = tolerance, action = action, strand = "reverse")
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


