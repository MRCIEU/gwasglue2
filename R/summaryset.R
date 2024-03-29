
#' An S4 class to represent the Summary Set
#'
#' @slot ss A tibble with the GWAS summary statistics (default NA).
#' @slot metadata  A list with the metadata associated to ss (default NA).
#' @slot variants The RSID/variants associated with ss (default NA).
#' @slot attributes  Attributes of the SummarySet. Eg. MR label Exposure/Outcome (default NA). 
#' @slot shape The shape of the SummarySet (default NA).
#' @slot shape The shape of the SummarySet (default NA).
#' * `"single"`: single region
#' * `"multiple"`: multiple regions
#' * `"independent"`: independent/scattered variants
#' * `"pruned"`: genome wide - pruned
#' * `"full"`: genome wide - full 
#' @export 
setClass("SummarySet",
  slots = c(
    ss = "tbl_df",
    metadata = "list", 
    variants = "character",
    attributes = "list",
    shape = "character"
  ),
  prototype = prototype(
    ss = NA_character_,
    metadata = list(NA),
    variants = NA_character_,
    attributes = list(NA),
    shape = NA_character_
  ),
  contains = class(dplyr::tibble())
)


#' SummarySet function
#'
#' @param sumstats GWAS summary statistics
#' @keywords internal
#' @importFrom methods new
#' @return A gwasglue2 SummarySet object.
#' @export
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
#' @return gwasglue2 SummarySet object with metadata stored.
#' @export
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

#' Add to metadata in the SummarySet 
#'
#' @param summary_set A gwasglue2 SummarySet object.
#' @param id GWAS study ID.
#' @param sample_size Sample size.
#' @param nsnp Number of variants in the study.
#' @param trait  Phenotype name corresponding the the variant.
#' @param sd Trait standard deviation.
#' @param unit Unit.
#' @param ncontrol Number of controls in study.
#' @param build   genome build version.
#' @param population  Study sample population.
#' @param ncase Number of cases in study.
#' @return gwasglue2 SummarySet object with metadata stored.
#' @export
#' @docType methods
#' @rdname addToMetadata-methods
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
#' @rdname addToMetadata-methods
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
#' @return The gwasglue2 SummarySet metadata.
#' @export
#' @docType methods
#' @rdname getMetadata-methods
setGeneric("getMetadata", function(summary_set) standardGeneric("getMetadata"))
#' @rdname getMetadata-methods
setMethod("getMetadata", "SummarySet",
          function(summary_set) {
            return(summary_set@metadata)
          })


# #' Get Method to retrieve the source information of the GWAS Summary Statistics stored in the SummarySet 
# #'
# #' @param summary_set A gwasglue2 SummarySet object.
# #' @return The source information of the GWAS Summary Statistics (type of file and accession/creation date). 
# #' @export
# #' @docType methods
# #' @rdname getSource-methods
# setGeneric("getSource", function(summary_set) standardGeneric("getSource"))
# #' @rdname getSource-methods
# setMethod("getSource", "SummarySet",
#           function(summary_set) {
#             return(summary_set@source)
#           })

#' Set Method to create an internal Variant ID for the SummarySet 
#'
#' @param summary_set A gwasglue2 SummarySet object
#' @return A extra '"variantid"' column in the  GWAS summary statistics tibble. The [getSummaryData()] can be used to retrieve it.
#' @export
#' @docType methods
#' @rdname setVariantid-methods
setGeneric("setVariantid", function(summary_set) standardGeneric("setVariantid"))
#' @rdname setVariantid-methods
setMethod("setVariantid", "SummarySet", function(summary_set) {
  chr <- getSummaryData(summary_set)$chr
  pos <- getSummaryData(summary_set)$position
  a1 <- getSummaryData(summary_set)$ea
  a2 <- getSummaryData(summary_set)$nea
  
  variantid <- create_variantid(chr,pos,a1,a2)

  summary_set@ss[,"variantid"] <- variantid

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


#' Set Method to store the attributes of the SummarySet
#'
#' @param summary_set A gwasglue2 SummarySet object
#' @param mr_label It can be either `"exposure"` or `"outcome". Default NULL.
#' @param ... Other attributes information
#' @return The gwasglue2 SummarySet object with the attributes stored
#' @export
#' @docType methods
#' @rdname setAttributes-methods
setGeneric("setAttributes",function(summary_set, mr_label = NULL, ...) standardGeneric("setAttributes"))
#' @rdname setAttributes-methods
setMethod( "setAttributes", "SummarySet",
           function(summary_set, mr_label = NULL, ...) {
             summary_set@attributes <- c(summary_set@attributes, mr_label = mr_label, ...)
             return(summary_set)
           }
)

#' Get Method to retrieve the attributes linked to the SummarySet
#'
#' @param summary_set A gwasglue2 SummarySet object
#' @return The attributes associated with the SummarySet
#' @export
#' @docType methods
#' @rdname getAttributes-methods
setGeneric("getAttributes",function(summary_set) standardGeneric("getAttributes"))
#' @rdname getAttributes-methods
setMethod("getAttributes","SummarySet",
          function(summary_set) {
            return(summary_set@attributes)
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


#'  Set the Shape of the gwasglue2 objects 
#' @param object A gwasglue2 SummarySet or DataSet object
#' @param shape The shape of the GWAS data
#' @return The gwasglue2 object with the shape stored
#' @export
#' @docType methods
#' @rdname setShape-methods
setGeneric("setShape",function(object, shape) standardGeneric("setShape"))
#' @rdname setShape-methods
setMethod("setShape", "SummarySet", function(object, shape) {
  
  # check if the shape is allowed
  shapes <- c("single", "multiple", "independent", "pruned", "full")
  if (shape %ni% shapes){
    stop( " This shape is not allowed in gwasglue2. The options available for the 'shape' of the 'SummarySet' object are 'single', 'multiple', 'independent', 'pruned',and 'full'.")
  } else {
    object@shape <- shape
  }
  return(object)
})

#'  Get the Shape of the gwasglue2 objects
#' @param object A gwasglue2 SummarySet or DataSet object.
#' @return The shape of the gwasglue2 object
#' @export
#' @docType methods
#' @rdname getShape-methods
setGeneric("getShape",function(object) standardGeneric("getShape"))
#' @rdname getShape-methods
setMethod("getShape", "SummarySet", function(object) {
  return(object@shape)
})


# show method 
setMethod(f = "show", signature="SummarySet", definition = function(object) {
  
  # set description of SummarySet
  id <- getMetadata(object)$id
  n <-  getMetadata(object)$sample_size
  nvariants <- nrow(getSummaryData(object))
  shape <- getShape(object)

  cat("\nA SummarySet with", nvariants, "variants.\n")
  cat("\nStudy ID:", id, "\n")
  cat("Sample size:", n, "\n")
  
  if (is.na(shape)){
  cat("Shape: No shape defined. Use the setShape() function to add it to the SummarySet\n")
  cat("NOTE: The shape feature is not fully implemented yet. Analyses can continue without defining it.\n")
  } else{
      cat("Shape:", shape, "\n")
  }
  
  cat("\nGWAS summary data: \n")
  print(getSummaryData(object))
  cat("\nTo access the GWAS summary data use the getSummaryData() function.\n")
})
