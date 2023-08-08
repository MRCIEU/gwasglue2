
#' Chain files downloader 
#' @description Human genome chain files are download from UCSC
#' @param from genome assembly to which GWAS summary data is currently mapped. Default "GRCh37". 
#' * Other options: "NCBI34", "NCBI35", "NCBI36" and "GRCh38"
#' @param  to genome  assembly to which should be mapped. Default "GRCh38"
#' * Other options: "NCBI34", "NCBI35", "NCBI36" and "GRCh37"
#' @return A chain file
#' @export
download_chainfile <- function(from = "GRCh37", to = "GRCh38") {
  file <- tempfile()
  utils::download.file(paste0("https://ftp.ensembl.org/pub/assembly_mapping/homo_sapiens/",from,"_to_", to,".chain.gz"),file)
   
# TODO: Need a more efficient way to decompress
  gz <- gzfile(file)
  writeLines(readLines(gz), file)
  close.connection(gz)

  return(file)
}


#' Remap genomic coordinates to a different genome assembly
#' @description Converts SummarySet GWAS summary data to a different genome assembly. Human genome chain files are download from ENSEMBL
#' @param summaryset A gwasglue2 SummarySet object
#' @param chainfile The chainfile used to remap the genomic coordinates. If `NULL` a chainfile is downloaded using the \code{to}  and the \code{summaryset}metadata to check the genome assembly to which GWAS summary data is currently mapped. 
#' @seealso [create_metadata()] and [addToMetadata()] on how to create or add to metadata.
#' Note that if using \code{chainfile} the analyses are not restricted to Human GWAS. 
#' @param  to genome  assembly to which should be mapped. Default "GRCh38"
#' * Other options: "NCBI34", "NCBI35", "NCBI36", GRCh37" and "GRCh38"
#' @return A gwasglue2 SummarySet object with GWAS summary data genomic coordinates remapped. 
#' @export 
liftover <- function(summaryset, chainfile = NULL, to = "GRCh38") {


  # check if bioconductor packages are installed
  if (!requireNamespace("rtracklayer", quietly =TRUE)  | !requireNamespace("IRanges", quietly =TRUE) | !requireNamespace("GenomicRanges", quietly =TRUE)){
    stop("The Bioconductor packages `rtracklayer`,`IRanges` and `GenomicRanges` need to installed to perform the genome build liftover.")}

  human_builds <- c("NCBI34", "NCBI35", "NCBI36","GRCh37","GRCh38")
  
  # check if the genome build assembly is present in the metadata
  from <- summaryset@metadata$build
  
  if (is.null(from) | is.na(from)){
      stop("There is no information about the SummarySet build in the metadata. Create again the SummarySet using the `build` argument or use the `addTOMetadata()` function to add the build to the metadata. Ex. `addTOMetadata(summaryset, build = 'GRCh37')` ")
  }
  
  if (!is.null(from) & !is.na(from) & from %ni% human_builds){
    stop( "Check if the  options for the `build` in the metadata are `NCBI34`, `NCBI35`, `NCBI36`, `GRCh37` or `GRCh38`. Create again the SummarySet using the `build` argument or use the `addTOMetadata()` function to add the build to the metadata. Ex. `addTOMetadata(summaryset, build = 'GRCh37')`")
  }

  # check if `chainfile` and `to` are not both NULL
  if(is.null(chainfile) && is.null(to)){
    stop("The `chainfile` and `to` arguments cannot be both NULL.")
  } 
  if(!is.null(chainfile) && is.null(to)){
    stop("The  `to` argument needs to be given to add to the metadata")
  }
   
   message("The  build to 'lift to' is ", to ,". Change the `to` argument if it is another build.")


  #  Download chainfile if not given
  if(is.null(chainfile) && !is.null(to)){
    if (is.null(to) || to %ni% human_builds){
      stop( "The options for argument `build_to` are `NCBI34`, `NCBI35`, `NCBI36`, `GRCh37` or `GRCh38`.")
    }
    chainfile <- download_chainfile(from = from, to = to)
  }


  message("Lifting build: ", from, " to ", to)
  
  # get GWAS summary data from summaryset
  dat <- as.data.frame(getSummaryData(summaryset))
  chr_col = "chr"
  pos_col = "position"
  ea_col = "ea"
  oa_col = "nea"

  message("Loading chainfile")
  ch <- rtracklayer::import.chain(chainfile)

  message("Converting chromosome codings")
  if(!grepl("chr", dat[[chr_col]][1])){
    dat[[chr_col]] <- paste0("chr", dat[[chr_col]])
  }

  dat[[chr_col]][dat[[chr_col]] == "chr23"] <- "chrX"
  dat[[chr_col]][dat[[chr_col]] == "chr24"] <- "chrY"
  dat[[chr_col]][dat[[chr_col]] == "chr25"] <- "chrXY"
  dat[[chr_col]][dat[[chr_col]] == "chr26"] <- "chrM"
  dat[[chr_col]][dat[[chr_col]] == "chrMT"] <- "chrM"

  message("Organising")
  datg <- GenomicRanges::GRanges(
    seqnames=dat[[chr_col]],
    ranges=IRanges::IRanges(start=dat[[pos_col]], end=dat[[pos_col]]),
    ind=1:nrow(dat)
  )

  message("Lifting")
  d19 <- rtracklayer::liftOver(datg, ch) %>% unlist()

  message("Organising again")

  dat <- dat[d19$ind,]
  dat[[chr_col]] <- d19@seqnames
  dat[[pos_col]] <- d19@ranges@start
  dat[[chr_col]] <- gsub("chr", "", dat[[chr_col]])

  message("Reordering")
  dat <- dat[order(dat[[chr_col]], dat[[pos_col]]), ]

  if(!is.null(ea_col) & !is.na(oa_col)){
    message("Removing duplicates")
    nom <- names(dat)

    if(is.numeric(chr_col)) chr_col <- nom[chr_col]
    if(is.numeric(pos_col)) pos_col <- nom[pos_col]
    if(is.numeric(ea_col)) ea_col <- nom[ea_col]
    if(is.numeric(oa_col)) oa_col <- nom[oa_col]

    dat <- dplyr::distinct(dat, .data[[chr_col]], .data[[pos_col]], .data[[ea_col]], .data[[oa_col]], .keep_all=TRUE)
  }

# replace columns 
  summaryset@ss$chr <- dat$chr
  summaryset@ss$position <- dat$position
  summaryset@ss$rsid <- dat$rsid
  summaryset@ss$ea <- dat$ea
  summaryset@ss$nea <- dat$nea

		
  # set variantid again
  summaryset <- summaryset %>% setVariantid(.)
  # add ifo to attributes
  summaryset@attributes$lift_from <- from
  summaryset@attributes$build <- to
  summaryset@attributes$time_liftover <-  Sys.time()

  # add info to metadata
  summaryset@metadata$lift_from <- from
  summaryset@metadata$build <- to
 
  message("Done")
  return(summaryset)

}

 


