
#' Chain files downloader 
#' @description Human genome chain files are download from UCSC
#' @param from genome assembly to which GWAS summary data is currently mapped. Default "GRCh37". 
#' * Other options: "NCBI34", "NCBI35", "NCBI36" and "GRCh38"
#' @param  to genome  assembly to which should be mapped. Default "GRCh38"
#' * Other options: "NCBI34", "NCBI35", "NCBI36" and "GRCh37"
#' @return A chain file
download_chainfile <- function(from = "GRCh37", to = "GRCh38") {
  file <- tempfile()
  utils::download.file(paste0("https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/",from,"_to_", to,".chain.gz"),file)
   
# TODO: Need a more efficient way to decompress
  gz <- gzfile(file)
  writeLines(readLines(gz), file)
  close.connection(gz)

  return(file)
}


 


