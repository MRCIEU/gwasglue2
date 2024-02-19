# YAML constructors


#' YAML constructor for the coloc job
#' @param traits A character vector of traits to be analysed.
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

job_coloc <- function(traits, region, build = "GRCh37", bfile = NULL, plink_bin = NULL, write = FALSE, outfile = "config.yaml"){
  
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


