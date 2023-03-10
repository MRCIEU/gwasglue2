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
createSumset <- function(traits, variants) {
  x <- ieugwasr::associations(variants = variants, id = traits)
  x <- dplyr::arrange(x,rsid)
  return(x)
}


# create s4 SummarySet objects and fill metadata slotand fill metadata slot
createSummarySets <- function (traits,variants, tools, source, ld_ref){
s <- SummarySet(traits = traits, variants = variants, tools = tools) %>%
    setMetadata(., source = source, traits = traits) %>%
    setRSID(.,.@ss$rsid)

    return(s)
   }

#' A function to create a DataSet object and harmonise data against data
#'
#' @param traits ID of GWAS studies to query
#' @param variants Array of SNPs rsids
#' @param tools Array of methods that gwasglue2 is going to convert the summarySet to (eg. "mr") 
#' @param source ID of Data Source to set the DataSet @slot metadata (Default and only option: IEUopenGWAS). 
#' @param harmonise logical (default TRUE). It harmonises the summary sets in the DataSet against each other. 
#' @param tolerance Inherited from harmoniseData() (default 0.08)
#' @param action Inherited from harmoniseData() (Default 1)
#'
#' @return A harmonised gwasglue2 DataSet object

createDataSet <- function(traits, variants, tools, source="IEUopenGWAS", harmonise = TRUE, tolerance = 0.08,action = 1){

  ds <- DataSet()
  for (i in seq_along(traits)){
    s <- createSummarySets(traits=traits[i],
                            variants = variants,
                            tools = tools,
                            source = source)
  ds@summary_sets[[i]] <- s
  }

  if (harmonise == TRUE){
    ds <-  ds %>%
    overlapSNP(.) %>%
    harmoniseData(.,tolerance = tolerance,action = action)
  }

  return(ds)
}

#' Plot
#'
#' @param dataset gwasglue2 DataSet object
#' @param type Type of plot (Only available "manhattan" plots at the moment)
#' @param region Main title for the plot
#'
#' @return A plot
plot_gwasglue <- function(dataset, type, region){
  
  if(type == "manhattan"){
    ntraits <- getLength(dataset)
    nb_rows <- ceiling(ntraits/2)
     # Add main title
    
    
    par(mfrow=c(nb_rows, 2))
    
    
      for (i in 1:ntraits){
      plot(dataset@summary_sets[[i]]@ss$position, -log10(dataset@summary_sets[[i]]@ss$p), main = dataset@summary_sets[[i]]@metadata$trait, xlab = "position", ylab = "-log10(p-value)", pch=20 )
      }
    mtext(region,                  
          side = 3,
          line = - 2,
          outer = TRUE)
  }
}


