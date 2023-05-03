
#' An S4 class to represent the Summary Set
#'
#' @slot source Eg.IEUopenGWAS (default NA).
#' @slot ss A tibble with the GWAS summary statistics (default NA).
#' @slot metadata  A data.frame with the metadata associated to ss (default NA).
#' @slot variants The RSID/variants associated with ss (default NA).
#' @slot tools The tools that gwasglue2 is going to convert to (default NA).
#' @slot mr_label Exposure/Outcome (default NA).
setClass("SummarySet",
  slots = c(
    source = "character",
    ss = "tbl_df",
    metadata = "data.frame",
    variants = "character",
    tools = "character",
    mr_label = "character"
  ),
  prototype = prototype(
    source = NA_character_,
    ss = NA_character_,
    metadata = data.frame(NA),
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
#' @param object A SummarySet  S4 object
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
# Set and get methods for Metadata and Source


#' Set Method to store metadata in the SummarySet 
#'
#' @param object A SummarySet  S4 object
#' @param metadata A dataframe with metadata information
#' @return A  SummarySet  S4 object with slot metadata stored
#' @export
#' @docType methods
#' @rdname setMetadata-methods
setGeneric("setMetadata", function(object, metadata) standardGeneric("setMetadata"))
#' @rdname setMetadata-methods
setMethod("setMetadata", "SummarySet",
          function(object,metadata) {
            object@metadata <- as.data.frame(metadata)
            return(object)
          })
          
#' Get Method to retrieve the metadata stored in the SummarySet 
#'
#' @param object A SummarySet  S4 object
#' @return The metadata
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
#' @param object A SummarySet  S4 object
#' @return A extra '"variantid"' column in the  slot @ss
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
 





# Set and get methods for RSID
setGeneric("setRSID",function(object,variants) standardGeneric("setRSID"))
setMethod( "setRSID", "SummarySet",
           function(object,variants) {
             object@variants <- variants
             return(object)
           }
)

setGeneric("getRSID",function(object) standardGeneric("getRSID"))
setMethod("getRSID","SummarySet",
          function(object) {
            return(object@variants)
          }
)

# Set and get methods for tools
setGeneric("setTool",function(object,tools) standardGeneric("setTool"))
setMethod( "setTool", "SummarySet",
           function(object,tools) {
             object@tools <- tools
             message(paste("Gwasglue is going to convert data to ",object@tools))
             return(object)
           }
)

setGeneric("getTool",function(object) standardGeneric("getTool"))
setMethod("getTool","SummarySet",
          function(object) {
            return(object@tools)
          }
)


# Set and get methods for mr_label
setGeneric("setMRlabel",function(object,mr_label) standardGeneric("setMRlabel"))
setMethod( "setMRlabel", "SummarySet",
           function(object,mr_label) {
             object@mr_label <- mr_label
             return(object)
           }
)

setGeneric("getMRlabel",function(object) standardGeneric("getMRlabel"))
setMethod("getMRlabel","SummarySet",
          function(object) {
            return(object@mr_label)
          }
)




setGeneric("getZscores",function(object) standardGeneric("getZscores"))
setMethod("getZscores","SummarySet",
          function(object) {
            return(object@zscores)
          }
)

setGeneric("dimData",function(x) standardGeneric("dimData"))
setMethod("dimData", signature=c(x="SummarySet"), function(x) {
    return(dim(x@ss))
})

setMethod(f = "show", signature="SummarySet", definition = function(object) {
  cat("A SummarySet with ", nrow(object@ss), " variants\n")
  print(object@ss)
})
