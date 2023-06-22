
#' An S4 class to represent the Summary Set
#'
#' @slot source Eg.IEUopenGWAS (default NA).
#' @slot ss A tibble with the GWAS summary statistics (default NA).
#' @slot metadata  A list with the metadata associated to ss (default NA).
#' @slot variants The RSID/variants associated with ss (default NA).
#' @slot mr_label Exposure/Outcome (default NA).
setClass("SummarySet",
  slots = c(
    source = "list",
    ss = "tbl_df",
    metadata = "list", 
    variants = "character",
    mr_label = "character"
  ),
  prototype = prototype(
    source = list(NA),
    ss = NA_character_,
    metadata = list(NA),
    variants = NA_character_,
    mr_label = NA_character_
  ),
  contains = class(dplyr::tibble())
)


#' SummarySet function
#'
#' @param sumstats GWAS summary statistics
#' @importFrom methods new
#' @return A gwasglue2 SummarySet object.
SummarySet <- function(sumstats) {
  new("SummarySet",
    ss = sumstats
  )
}

#' Get Method to retrieve the GWAS Summary Statistics  from the SummarySet 
#'
#' @param summary_set A gwasglue2 SummarySet object.
#' @return A tibble with GWAS summary statistics
#' @export
#' @docType methods
#' @rdname getSummaryData-methods
setGeneric("getSummaryData", function(summary_set) standardGeneric("getSummaryData"))
#' @rdname getSummaryData-methods
setMethod("getSummaryData", "SummarySet",
          function(summary_set) {
            return(summary_set@ss)

          })


#' Set Method to add metadata to the SummarySet 
#'
#' @param summary_set A gwasglue2 SummarySet object.
#' @param metadata A list with metadata information. 
#' @docType methods
#' @rdname setMetadata-methods
setGeneric("setMetadata", function(summary_set, 
                           metadata) standardGeneric("setMetadata"))
#' @rdname setMetadata-methods
setMethod("setMetadata", "SummarySet",
          function(summary_set,
            metadata) {
            
 summary_set@metadata <- metadata

  return(summary_set)
})

#' Method to add metadata to the SummarySet 
#'
#' @param summary_set A gwasglue2 SummarySet object.
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
#' @return gwasglue2 SummarySet object with metadata stored
#' @export
#' @docType methods
#' @rdname setMetadata-methods
setGeneric("addToMetadata", function(summary_set, 
                           id = getMetadata(summary_set)$id,
                           sample_size = getMetadata(summary_set)$sample_size,
                           nsnp = getMetadata(summary_set)$nsnp,
                           trait = getMetadata(summary_set)$trait,
                           sd = getMetadata(summary_set)$sd,
                           unit = getMetadata(summary_set)$unit,
                           ncontrol = getMetadata(summary_set)$ncontrol, 
                           build = getMetadata(summary_set)$build,
                           population = getMetadata(summary_set)$population,
                           ncase = getMetadata(summary_set)$ncase) standardGeneric("addToMetadata"))
#' @rdname setMetadata-methods
setMethod("addToMetadata", "SummarySet",
          function(summary_set,
            id,
            sample_size,
            nsnp,
            trait,
            sd,
            unit,
            ncontrol, 
            build,
            population,
            ncase) {
            
 summary_set@metadata$id <- id
 summary_set@metadata$sample_size <- sample_size
 summary_set@metadata$nsnp <- nsnp
 summary_set@metadata$trait <- trait
 summary_set@metadata$sd <- sd
 summary_set@metadata$unit <- unit
 summary_set@metadata$ncontrol <- ncontrol
 summary_set@metadata$build <- build
 summary_set@metadata$population <- population
 summary_set@metadata$ncase <- ncase

  return(summary_set)
})


#' Get Method to retrieve the metadata stored in the SummarySet 
#'
#' @param summary_set A gwasglue2 SummarySet object.
#' @return The gwasglue2 SummarySet metadata
#' @export
#' @docType methods
#' @rdname getMetadata-methods
setGeneric("getMetadata", function(summary_set) standardGeneric("getMetadata"))
#' @rdname getMetadata-methods
setMethod("getMetadata", "SummarySet",
          function(summary_set) {
            return(summary_set@metadata)
          })




setGeneric("getSource", function(summary_set) standardGeneric("getSource"))
setMethod("getSource", "SummarySet",
          function(summary_set) {
            return(summary_set@source)
          })

#' Set Method to create an internal Variant ID for the SummarySet 
#'
#' @param summary_set A gwasglue2 SummarySet object
#' @return A extra '"variantid"' column in the  GWAS summary statistics tibble. The @rdname getSummaryData-methods can be used to retrieve it.
#' @export
#' @docType methods
#' @rdname setVariantid-methods
setGeneric("setVariantid", function(summary_set) standardGeneric("setVariantid"))
#' @rdname setVariantid-methods
setMethod("setVariantid", "SummarySet", function(summary_set) {
  sumstats <- getSummaryData(summary_set)
  nvariants <- dim(sumstats)[1] 
  variantid <- lapply(1:nvariants, function(i){
    x <- sort(c(sumstats[i,]$ea,sumstats[i,]$nea))
    if (nchar(x[1]) > 10 || nchar(x[2]) <= 10){
      id <- paste0(sumstats[i,]$chr,":", sumstats[i,]$position,"_#",digest::digest(x[1],algo= "murmur32"),"_",x[2]) 
    }
    if (nchar(x[1]) <= 10 || nchar(x[2]) > 10){
      id <- paste0(sumstats[i,]$chr,":", sumstats[i,]$position,"_",x[1],"_#",digest::digest(x[2],algo= "murmur32")) 
    }
    if (nchar(x[1]) > 10 || nchar(x[2]) > 10){
      id <- paste0(sumstats[i,]$chr,":", sumstats[i,]$position,"_#",digest::digest(x[1],algo= "murmur32"),"_#",digest::digest(x[2],algo= "murmur32")) 
    } else {
      id <- paste0(sumstats[i,]$chr,":", sumstats[i,]$position,"_",x[1],"_",x[2]) 
    }
  }) 

  summary_set@ss[,"variantid"] <- unlist(variantid)
  return(summary_set)
})
 

#' Set Method to store RSID/variants in the SummarySet
#'
#' @param summary_set A gwasglue2 SummarySet object
#' @param variants The RSID/variants associated with the GWAS summary statistics
#' @seealso Similar to [setVariants()]
#' @return The gwasglue2 SummarySet object with RSID/variants stored
#' @export
#' @docType methods
#' @rdname setRSID-methods
setGeneric("setRSID",function(summary_set,variants) standardGeneric("setRSID"))
#' @rdname setRSID-methods
setMethod( "setRSID", "SummarySet",
           function(summary_set,variants) {
             summary_set@variants <- variants
             return(summary_set)
           }
)

#' Get Method to retrieve RSID/variants stored in the SummarySet
#'
#' @param summary_set A gwasglue2 SummarySet object
#' @seealso Similar to [getVariants()]
#' @return The RSID/variants
#' @export
#' @docType methods
#' @rdname getRSID-methods
setGeneric("getRSID",function(summary_set) standardGeneric("getRSID"))
#' @rdname getRSID-methods
setMethod("getRSID","SummarySet",
          function(summary_set) {
            return(summary_set@variants)
          }
)

#' Set Method to store RSID/variants in the SummarySet
#'
#' @param summary_set A gwasglue2 SummarySet object
#' @param variants The RSID/variants associated with the GWAS summary statistics
#' @return The gwasglue2 SummarySet object with RSID/variants stored
#' @seealso Similar to [setRSID()]
#' @export
#' @docType methods
#' @rdname setVariants-methods
setGeneric("setVariants",function(summary_set,variants) standardGeneric("setVariants"))
#' @rdname setVariants-methods
setMethod( "setVariants", "SummarySet",
           function(summary_set,variants) {
             summary_set@variants <- variants
             return(summary_set)
           }
)

#' Get Method to retrieve RSID/variants stored in the SummarySet
#'
#' @param summary_set A gwasglue2 SummarySet object
#' @seealso Similar to [getRSID()]
#' @return The RSID/variants
#' @export
#' @docType methods
#' @rdname getVariants-methods
setGeneric("getVariants",function(summary_set) standardGeneric("getVariants"))
#' @rdname getVariants-methods
setMethod("getVariants","SummarySet",
          function(summary_set) {
            return(summary_set@variants)
          }
)


#' Set Method to store the Mendelian randomization (MR) label in the SummarySet
#'
#' @param summary_set A gwasglue2 SummarySet object
#' @param mr_label MR label for SummarySet. It can be either `"exposure"` or `"outcome"`
#' @return The gwasglue2 SummarySet object with the MR labels stored
#' @export
#' @docType methods
#' @rdname setMRlabel-methods
setGeneric("setMRlabel",function(summary_set,mr_label) standardGeneric("setMRlabel"))
#' @rdname setMRlabel-methods
setMethod( "setMRlabel", "SummarySet",
           function(summary_set,mr_label) {
             summary_set@mr_label <- mr_label
             return(summary_set)
           }
)

#' Get Method to retrieve the Mendelian randomization (MR) label linked to the SummarySet
#'
#' @param summary_set A gwasglue2 SummarySet object
#' @return The MR label associated with the SummarySet
#' @export
#' @docType methods
#' @rdname getMRlabel-methods
setGeneric("getMRlabel",function(summary_set) standardGeneric("getMRlabel"))
#' @rdname getMRlabel-methods
setMethod("getMRlabel","SummarySet",
          function(summary_set) {
            return(summary_set@mr_label)
          }
)




# setGeneric("getZscores",function(summary_set) standardGeneric("getZscores"))
# setMethod("getZscores","SummarySet",
#           function(summary_set) {
#             return(summary_set@zscores)
#           }
# )


#' Dimensions of the GWAS Summary Statistics data
#'
#' @param summary_set A gwasglue2 SummarySet object
#' @return The dimensions of the GWAS Summary Statistics data
#' @export
#' @docType methods
#' @rdname dimData-methods
setGeneric("dimData",function(summary_set) standardGeneric("dimData"))
#' @rdname dimData-methods
setMethod("dimData", "SummarySet", function(summary_set) {
    return(dim(summary_set@ss))
})





setMethod(f = "show", signature="SummarySet", definition = function(object) {
  cat("A SummarySet with ", nrow(object@ss), " variants\n")
  print(object@ss)
})
