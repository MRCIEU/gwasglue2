

#' gwasglue 
#' @param yaml A yaml R object.
#' @param read A boolean indicating whether the yaml object is a path to a yaml file or a yaml R object. IF TRUE, yaml is a path.
#' @param path_to_yaml A character vector indicating the path to the yaml file.
#' @return A gwasglue2 DataSet object.
#' @importFrom yaml read_yaml yaml.load
#' @export 
gwasglue <- function(yaml, read=FALSE, path_to_yaml="config.yaml"){
  if (!requireNamespace("yaml", quietly =TRUE)){
    stop("The R package `yaml` needs to installed.")
    }
  # read or load yaml
  if (isTRUE(read)) { yaml <- yaml::read_yaml(path_to_yaml) }
  else { yaml <- yaml::yaml.load(yaml) }

  # get nb of datasets
  ndatasets <- length(yaml) -1 # -1 because of the analyses part

  dataset <- lapply(1:ndatasets, \(i){
    dataset <- yaml[[i]]
          
    # get variants
    variants <- dataset$variants$variant_list
    if(grepl(".RData", variants)){
     load(variants)
    }

    # get shape
    shape <- dataset$variants$shape
    
    # get the number of SummarySets and create a list
    nsummary_sets <- length(dataset$summarydata)

    summary_sets <- lapply (1:nsummary_sets, \(j){
      # read the parameters for each SummarySet
      source <- dataset$summarydata[[j]]$source
      location <- dataset$summarydata[[j]]$location
      build <- dataset$summarydata[[j]]$build
      organization <- dataset$summarydata[[j]]$organization
          
      if(source == "opengwas"){
          if (!requireNamespace("ieugwasr", quietly =TRUE)){
            stop("The MRC IEU R package `ieugwasr` needs to be installed.")}
        
        summary_sets <- ieugwasr::associations(variants = variants, id = location, gwasglue = TRUE) %>% setShape(.,shape = shape)
      }
      if(source == "gwascatalog"){
        if (!requireNamespace("readr", quietly =TRUE)){
              stop("The CRAN package `readr` needs to be installed.")}
        data <- readr::read_tsv(get(location))
        summary_sets <- create_summaryset(data, variants = variants, qc =TRUE, build = build) %>% setShape(.,shape = shape)
      }
      if(source == "vcf"){
        if (!requireNamespace("gwasvcf", quietly =TRUE)){
          stop("The MRC IEU R package `gwasvcf` needs to be installed.")}
        
      summary_sets <- gwasvcf::readVcf(location) %>% 
                gwasvcf::query_gwas(chrompos = variants) %>% gwasvcf::gwasvcf_to_summaryset()%>% setShape(.,shape = shape)
      }
      if(source == "summaryset"){
        summary_sets <- get(location) %>% setShape(.,shape = shape)
        if(!inherits(summary_sets,"SummarySet")){
          stop("The object ", location, " is not a SummarySet.")
        }
      }
      if(source == "robject"){
        data <- dplyr::as_tibble(get(location))
        if(!inherits(data, c("tbl", "tbl_df", "data.frame"))){
            stop("The object ", location, " is not a `data.frame`.`")
          }
        
      summary_sets <- create_summaryset(data, variants = variants, qc =TRUE, build = build) %>% setShape(.,shape = shape)
      }
      if(source == "text"){
        if (!requireNamespace("fread", quietly =TRUE)){
              stop("The CRAN package `fread` needs to be installed.")}
        data <- data.table::fread(location, header = TRUE)
        summary_sets <- create_summaryset(data, variants = variants, qc =TRUE, build = build) %>% setShape(.,shape = shape)
      }
      else{
        stop("The source ", source, " is not recognized.")}


      # set attributes
      if(yaml$analyses$type == "MR"){
        summary_sets <- setAttributes(summary_sets, mr_label = organization)
      }

      return(summary_sets)
      })

    # create dataset
    dataset <- create_dataset(summary_sets=summary_sets) 
      
    if(yaml$analyses$type == "LDScores"){
        dataset <- dataset %>% 
              # set zscores and chisq
              setZscores() %>%
              setChisq()
      }

    # get lddata info and harmonise againts reference panel
    if(!is.null(dataset$lddata)){
        bfile <-dataset$lddata$location
        plink_bin <- dataset$lddata$plink_bin
        dataset <- harmonise_ld(dataset, bfile = bfile, plink_bin = plink_bin)
    }
    return(dataset)
  })

  if(ndatasets > 1){
    return(ListDataSets(dataset))}
  else{
    return(dataset)
  }
    
}