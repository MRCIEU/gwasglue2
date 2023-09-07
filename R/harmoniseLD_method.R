#  The S4 methods in this file use the functions in `harmonise_ld.R`


#  Build LD matrix
setGeneric("buildLDMatrix", function(dataset, ...) standardGeneric("buildLDMatrix"))
setMethod("buildLDMatrix", "DataSet", function(dataset, bfile = NULL, plink_bin = NULL){
  
  message("Building LD matrix")
  variants<-cbind(dataset@summary_sets[[1]]@ss$chr, dataset@summary_sets[[1]]@ss$position, dataset@summary_sets[[1]]@ss$position, 1:nrow(dataset@summary_sets[[1]]@ss))
    
  dataset@ld_matrix <- ld_matrix_local(variants, bfile=bfile, plink_bin=plink_bin)

  variants_avail <- rownames(dataset@ld_matrix)
  message("\nData available for ", length(variants_avail), " variants")
  dataset@describe$refpop_variants_avail <- length(variants_avail)

  return(dataset)
}
)


#' Get Method to retrieve the Linkage Disequilibrium matrix
#'
#' @param dataset A gwasglue2 DataSet object
#' @return The LD matrix
#' @export
#' @docType methods
#' @rdname getLDMatrix-methods
setGeneric("getLDMatrix",function(dataset) standardGeneric("getLDMatrix"))
#' @rdname getLDMatrix-methods
setMethod("getLDMatrix", "DataSet",
          function(dataset) {
            return(dataset@ld_matrix)
          })

# Harmonise against LD matrix
setGeneric("harmoniseLDMatrix", function(dataset) standardGeneric("harmoniseLDMatrix"))
setMethod( "harmoniseLDMatrix", "DataSet", function(dataset) {
  
  variants_avail <- do.call(rbind, strsplit(rownames(dataset@ld_matrix), split="\\|"))[,1]

  message("Gwasglue is now harmonising the SummarySets against the LD matrix")
  for (i in seq_along(dataset@summary_sets)){
    # message("Gwasglue is now harmonising ", dataset@summary_sets[[i]]@metadata$id, " the against LD matrix!")
 
    # subsetting
    sub_ss <- subset(dataset@summary_sets[[i]]@ss, dataset@summary_sets[[i]]@ss$variantid %in% variants_avail)
    
    index <- match(dataset@summary_sets[[i]]@ss$variantid, variants_avail)
    ld <- dataset@ld_matrix[index, index]


    h <- harmonise_ld_dat(sub_ss,ld)
    dataset@summary_sets[[i]]@ss <- h[[1]]
    

    message("\nThere are ", dim(h[[1]])[1], " variants remaining after harmonising ", dataset@summary_sets[[i]]@metadata$id, " against the LD matrix.")
  }
  
  dataset@ld_matrix <- h[[2]]
  dataset@is_harmonisedLD <- TRUE

  dataset@describe$variants_after_LDharmonization <- nrow(getLDMatrix(dataset))
  
  message("Done!")
  return(dataset)
}
)


#' Check if the DataSet is harmonised against LD matrix
#' 
#' @param dataset A gwasglue2 DataSet object
#' @return TRUE/FALSE
#' @export 
#' @docType methods
#' @rdname isHarmonisedLD-methods
setGeneric("isHarmonisedLD",function(dataset) standardGeneric("isHarmonisedLD"))
#' @rdname isHarmonisedLD-methods
setMethod("isHarmonisedLD","DataSet",function(dataset) {
  return(dataset@is_harmonisedLD)
})


#' Harmonise data against LD matrix
#' 
#' Function to create a LDmatrix gwasglue2 object and set the @slot ld_matrix u
#' @param dataset The DataSet gwasglue2 object
#' @param bfile It corresponds to the path and prefix of the plink files used to build the LD correlation matrix. 
#' @param plink_bin Path to the plink executable
#' @return  The DataSet gwasglue2 object harmonised
#' @export 
harmonise_ld <- function(dataset, bfile = NULL, plink_bin = NULL){

    ld <- buildLDMatrix(dataset, bfile = bfile, plink_bin = plink_bin) 
    
    ds <- harmoniseLDMatrix(ld)
    return(ds)
}
