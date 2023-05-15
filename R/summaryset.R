
#' An S4 class to represent the Summary Set
#'
#' @slot source Eg.IEUopenGWAS (default NA).
#' @slot ss A tibble with the GWAS summary statistics (default NA).
#' @slot metadata  A list with the metadata associated to ss (default NA).
#' @slot variants The RSID/variants associated with ss (default NA).
#' @slot tools The tools that gwasglue2 is going to convert to (default NA).
#' @slot mr_label Exposure/Outcome (default NA).
setClass("SummarySet",
  slots = c(
    source = "character",
    ss = "tbl_df",
    metadata = "list", 
    variants = "character",
    tools = "character",
    mr_label = "character"
  ),
  prototype = prototype(
    source = NA_character_,
    ss = NA_character_,
    metadata = list(NA),
    variants = NA_character_,
    tool = NA_character_,
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
#' @param object A gwasglue2 SummarySet object.
#' @return A tibble with GWAS summary statistics
#' @export
#' @docType methods
#' @rdname getSummaryData-methods
setGeneric("getSummaryData", function(object) standardGeneric("getSummaryData"))
#' @rdname getSummaryData-methods
setMethod("getSummaryData", "SummarySet",
          function(object) {
            return(object@ss)

          })


#' Set Method to store metadata in the SummarySet 
#'
#' @param object A gwasglue2 SummarySet object.
#' @param metadata A dataframe with metadata information
#' @return gwasglue2 SummarySet object with metadata stored
#' @export
#' @docType methods
#' @rdname setMetadata-methods
setGeneric("setMetadata", function(object, metadata) standardGeneric("setMetadata"))
#' @rdname setMetadata-methods
setMethod("setMetadata", "SummarySet",
          function(object,metadata) {
            object@metadata <- as.list(metadata)
            return(object)
          })
          
#' Get Method to retrieve the metadata stored in the SummarySet 
#'
#' @param object A gwasglue2 SummarySet object.
#' @return The gwasglue2 SummarySet metadata
#' @export
#' @docType methods
#' @rdname getMetadata-methods
setGeneric("getMetadata", function(object) standardGeneric("getMetadata"))
#' @rdname getMetadata-methods
setMethod("getMetadata", "SummarySet",
          function(object) {
            return(object@metadata)
          })




setGeneric("getSource", function(object) standardGeneric("getSource"))
setMethod("getSource", "SummarySet",
          function(object) {
            return(object@source)
          })

#' Set Method to create an internal Variant ID for the SummarySet 
#'
#' @param object A gwasglue2 SummarySet object
#' @return A extra '"variantid"' column in the  GWAS summary statistics tibble. The @rdname getSummaryData-methods can be used to retrieve it.
#' @export
#' @docType methods
#' @rdname setVariantid-methods
setGeneric("setVariantid", function(object) standardGeneric("setVariantid"))
#' @rdname setVariantid-methods
setMethod("setVariantid", "SummarySet", function(object) {
  sumstats <- getSummaryData(object)
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

  object@ss[,"variantid"] <- unlist(variantid)
  return(object)
})
 

#' Set Method to store RSID/variants in the SummarySet
#'
#' @param object A gwasglue2 SummarySet object
#' @param variants The RSID/variants associated with the GWAS summary statistics
#' @seealso Similar to [setVariants()]
#' @return The gwasglue2 SummarySet object with RSID/variants stored
#' @export
#' @docType methods
#' @rdname setRSID-methods
setGeneric("setRSID",function(object,variants) standardGeneric("setRSID"))
#' @rdname setRSID-methods
setMethod( "setRSID", "SummarySet",
           function(object,variants) {
             object@variants <- variants
             return(object)
           }
)

#' Get Method to retrieve RSID/variants stored in the SummarySet
#'
#' @param object A gwasglue2 SummarySet object
#' @seealso Similar to [getVariants()]
#' @return The RSID/variants
#' @export
#' @docType methods
#' @rdname getRSID-methods
setGeneric("getRSID",function(object) standardGeneric("getRSID"))
#' @rdname getRSID-methods
setMethod("getRSID","SummarySet",
          function(object) {
            return(object@variants)
          }
)

#' Set Method to store RSID/variants in the SummarySet
#'
#' @param object A gwasglue2 SummarySet object
#' @param variants The RSID/variants associated with the GWAS summary statistics
#' @return The gwasglue2 SummarySet object with RSID/variants stored
#' @seealso Similar to [setRSID()]
#' @export
#' @docType methods
#' @rdname setVariants-methods
setGeneric("setVariants",function(object,variants) standardGeneric("setVariants"))
#' @rdname setVariants-methods
setMethod( "setVariants", "SummarySet",
           function(object,variants) {
             object@variants <- variants
             return(object)
           }
)

#' Get Method to retrieve RSID/variants stored in the SummarySet
#'
#' @param object A gwasglue2 SummarySet object
#' @seealso Similar to [getRSID()]
#' @return The RSID/variants
#' @export
#' @docType methods
#' @rdname getVariants-methods
setGeneric("getVariants",function(object) standardGeneric("getVariants"))
#' @rdname getVariants-methods
setMethod("getVariants","SummarySet",
          function(object) {
            return(object@variants)
          }
)


#' Set Method to store the tools that gwasglue2 is going to convert the SummarySet to
#'
#' @param object A gwasglue2 SummarySet object
#' @param tools The tools
#' @return The gwasglue2 SummarySet object with the tools stored
#' @export
#' @docType methods
#' @rdname setTool-methods
setGeneric("setTool",function(object,tools) standardGeneric("setTool"))
#' @rdname setTool-methods
setMethod( "setTool", "SummarySet",
           function(object, tools) {
             if(!is.null(tools)){
              object@tools <- tools
             message(paste("Gwasglue is going to convert data to ", list(object@tools)))
          }  
             return(object)
           }
)
#' Get Method to retrieve the tools stored in the SummarySet
#'
#' @param object A gwasglue2 SummarySet object
#' @return The tools
#' @export
#' @docType methods
#' @rdname getTool-methods
setGeneric("getTool",function(object) standardGeneric("getTool"))
#' @rdname getTool-methods
setMethod("getTool","SummarySet",
          function(object) {
            return(object@tools)
          }
)


#' Set Method to store the Mendelian randomization (MR) label in the SummarySet
#'
#' @param object A gwasglue2 SummarySet object
#' @param mr_label MR label for SummarySet. It can be either `"exposure"` or `"outcome"`
#' @return The gwasglue2 SummarySet object with the MR labels stored
#' @export
#' @docType methods
#' @rdname setMRlabel-methods
setGeneric("setMRlabel",function(object,mr_label) standardGeneric("setMRlabel"))
#' @rdname setMRlabel-methods
setMethod( "setMRlabel", "SummarySet",
           function(object,mr_label) {
             object@mr_label <- mr_label
             return(object)
           }
)

#' Get Method to retrieve the Mendelian randomization (MR) label linked to the SummarySet
#'
#' @param object A gwasglue2 SummarySet object
#' @return The MR label associated with the SummarySet
#' @export
#' @docType methods
#' @rdname getMRlabel-methods
setGeneric("getMRlabel",function(object) standardGeneric("getMRlabel"))
#' @rdname getMRlabel-methods
setMethod("getMRlabel","SummarySet",
          function(object) {
            return(object@mr_label)
          }
)




# setGeneric("getZscores",function(object) standardGeneric("getZscores"))
# setMethod("getZscores","SummarySet",
#           function(object) {
#             return(object@zscores)
#           }
# )


#' Dimensions of the GWAS Summary Statistics data
#'
#' @param object A gwasglue2 SummarySet object
#' @return The dimensions of the GWAS Summary Statistics data
#' @export
#' @docType methods
#' @rdname dimData-methods
setGeneric("dimData",function(object) standardGeneric("dimData"))
#' @rdname dimData-methods
setMethod("dimData", "SummarySet", function(object) {
    return(dim(object@ss))
})





setMethod(f = "show", signature="SummarySet", definition = function(object) {
  cat("A SummarySet with ", nrow(object@ss), " variants\n")
  print(object@ss)
})
