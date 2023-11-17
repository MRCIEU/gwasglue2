#' YAML constructor (WIP)
#' @param traits A character vector of traits to be analysed.
#' @param region A character vector of regions to be analysed.
#' @param build  A character vector  of reference genome assembly to generate the vcf file. Default is NULL. 
#' * Options are `"NCBI34"`, `"NCBI35"`, `"NCBI36"`, `"GRCh37"` or "GRCh38".
#' @param bfile It corresponds to the path and prefix of the plink files used to build the LD correlation matrix. 
#' @param plink_bin Path to the plink executable
#' @param write A boolean indicating whether the yaml object should be written to a file.
#' @param outfile A character vector indicating the path to write the yaml file.
#' @return A yaml object.
#' @importFrom yaml write_yaml as.yaml

plan_coloc <- function(traits, region, build = NULL, bfile = NULL, plink_bin = NULL, write = FALSE, outfile = "config.yaml") {
  
  # read number of DataSets to  build and analyse TODO: should we allow for more than 1 dataset?
  ndatasets <- length(region)
  
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
    # TODO: 
    l$analysis <- "coloc"
    return(l)
  })
   
  names(out) <- paste0("plan", 1:length(region))

  # write yaml to file
  if (isTRUE(write)){ yaml::write_yaml(out, outfile) }
# return a yaml object
  return(yaml::as.yaml(out))
}


#' gwasglue  (WIP)
#' @param yaml A yaml R object.
#' @param read A boolean indicating whether the yaml object is a path to a yaml file or a yaml R object. IF TRUE, yaml is a path.
#' @param path_to_yaml A character vector indicating the path to the yaml file.
#' @return A gwasglue2 DataSet object.
#' @importFrom yaml read_yaml yaml_load
gwasglue <- function(yaml, read=FALSE, path_to_yaml="config.yaml"){

  # read or load yaml
  if (isTRUE(read)) { yaml <- yaml::read_yaml(path_to_yaml) }
  else { yaml <- yaml::yaml.load(yaml) }

  # get nb of plans
  nplans <- length(yaml)
 
  for(i in 1:nplans){
    plan <- yaml[[i]]
        
    # get variants
    chr <- (strsplit(plan$variants$variant_list, ":") [[1]][1] %>% strsplit(., "chr"))[[1]][2]
    start_end <- strsplit(plan$variants$variant_list, ":") [[1]][2]
    variants <- paste0(chr,":",start_end)
  
    # get the number of SummarySets and create a list
    nsummary_sets <- length(plan$summarydata)
    summary_sets <- list()

    for (i in 1:nsummary_sets){
      # read the parameters for each SummarySet
      source <- plan$summarydata[[i]]$source
      location <- plan$summarydata[[i]]$location
      build <- plan$summarydata[[i]]$build
        
      if(source == "ieugwasr"){
       data <- ieugwasr::associations(variants = variants, id = location)
      }else{stop("source not supported yet")}

    # create summarysets
    summary_sets[[i]] <- create_summaryset(data,
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
                    build = build)

        }
  
  # get lddata info
  bfile <-plan$lddata$location
  plink_bin <- plan$lddata$plink_bin
  
  # create dataset
  dataset <- create_dataset(summary_sets=summary_sets) %>% 
            harmonise_ld(bfile = bfile, plink_bin = plink_bin)

    return(dataset)
  # TODO: For now, we only return the one plan, one DataSet (to change this we should return a list of DataSets). But should we allow for more than 1 dataset be created by the gwasglue function?
    
}

