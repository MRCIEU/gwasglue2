

#' gwasglue 
#' @param yaml A yaml R object.
#' @param read A boolean indicating whether the yaml object is a path to a yaml file or a yaml R object. IF TRUE, yaml is a path.
#' @param path_to_yaml A character vector indicating the path to the yaml file.
#' @return A gwasglue2 DataSet object.
#' @importFrom yaml read_yaml yaml.load
#' @export 
gwasglue <- function(yaml, read=FALSE, path_to_yaml="config.yaml"){
    if (!requireNamespace("yaml", quietly =TRUE)){
    stop("The R package `yaml` needs to installed to perform the genome build liftover.")}
  # read or load yaml
  if (isTRUE(read)) { yaml <- yaml::read_yaml(path_to_yaml) }
  else { yaml <- yaml::yaml.load(yaml) }

  # get nb of jobs
  njobs <- length(yaml) -1 # -1 because of the analyses part


dataset <- lapply(1:njobs, \(i){
    job <- yaml[[i]]
        
    # get variants
    if(job$variants$tophits ==FALSE || is.null(job$variants$tophits)){
      variants <- job$variants$variant_list
    }else{
      exposures <- job$organization$lhs
      variants <- get_tophits_from_datasource(exposures, summarydata=job$summarydata)
    }

 
    shape <- job$variants$shape
  
    # get the number of SummarySets and create a list
    nsummary_sets <- length(job$summarydata)

    summary_sets <- lapply (1:nsummary_sets, \(j){
      # read the parameters for each SummarySet
      source <- job$summarydata[[j]]$source
      location <- job$summarydata[[j]]$location
      build <- job$summarydata[[j]]$build
      label <- job$summarydata[[j]]$label
        
      if(source == "opengwas"){
       data <- ieugwasr::associations(variants = variants, id = location)
      }else{stop("source not supported yet")}

      # create summarysets
      summary_sets <- create_summaryset(data,
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
  
    # set attributes
    if(yaml$analyses$type == "TwoSampleMR"){
      # get the label of the exposure and outcome
      lhs <- yaml$organization$lhs
      rhs <- yaml$organization$rhs

      if(label %in% lhs){
        summary_sets <- setAttributes(summary_sets, mr_label = "exposure")
      } 
      if(label %in% rhs){
      summary_sets <- setAttributes(summary_sets, mr_label = "outcome")
      }
    }
    return(summary_sets)

    })
    # create dataset
    dataset <- create_dataset(summary_sets=summary_sets)
    # get lddata info and harmonise againts reference panel
    if(!is.null(job$lddata)){
      bfile <-job$lddata$location
      plink_bin <- job$lddata$plink_bin
      dataset <- harmonise_ld(dataset, bfile = bfile, plink_bin = plink_bin)
    }
    return(dataset)
  })



  if(njobs > 1){
    return(listDatasets(dataset))}
  else{
    return(dataset)
  }
    
} 





get_tophits_from_datasource <- function(exposures, summarydata){

  exposure <- lapply(summarydata, \(s){
    if(exposures %in% s$label){
      return(s$location[which(exposures %in% s$label)])
    }
    # else{
    #   message("Exposure label ", s, " not found in the YAML summarydata.")
    # }
    
  }) %>% unlist()

  tophits <- ieugwasr::tophits(exposure[1])$rsid

  # TODO: deal with more than one exposure
  return(tophits)
 }


# -------------------------------------------------------------------------



# ------------------------------------------------------------------------