setGeneric("buildLDMatrix", function(object, ...) standardGeneric("buildLDMatrix"))
setMethod("buildLDMatrix", "DataSet", function(object, ld_ref = NULL, pop = FALSE, bfile = FALSE, plink_bin = NULL){

     message("Building LD matrix")
      if (pop == TRUE){
        object@ld_matrix <- suppressWarnings(ieugwasr::ld_matrix(object@summary_sets[[1]]@ss$rsid, pop = ld_ref, with_alleles = TRUE, bfile = FALSE,plink_bin = plink_bin))
      }
      if (bfile == TRUE){
        object@ld_matrix <- suppressWarnings(ieugwasr::ld_matrix(object@summary_sets[[1]]@ss$rsid, pop=FALSE, with_alleles=TRUE, bfile=ld_ref, plink_bin=plink_bin))
      }

      rsid_avail <- rownames(object@ld_matrix)
      message("\nData available for ", length(rsid_avail), " variants")

  return(object)
}
)


# Get methods for LDMatrix
setGeneric("getLDMatrix",function(object) standardGeneric("getLDMatrix"))
setMethod("getLDMatrix", "DataSet",
          function(object) {
            return(object@ld_matrix)
          })




# Set and get methods for harmonise data against LD matrix
setGeneric("isHarmonisedLD",function(object) standardGeneric("isHarmonisedLD"))
setMethod("isHarmonisedLD","DataSet",function(object) {
  return(object@is_harmonisedLD)
}
)

setGeneric("harmoniseLDMatrix", function(object) standardGeneric("harmoniseLDMatrix"))
setMethod( "harmoniseLDMatrix", "DataSet", function(object) {
  
  rsid_avail <- do.call(rbind, strsplit(rownames(object@ld_matrix), split="_"))[,1]
  message("Gwasglue is now harmonising the SummarySets against the LD matrix")
  for (i in seq_along(object@summary_sets)){
    # message("Gwasglue is now harmonising ", object@summary_sets[[i]]@metadata$id, " the against LD matrix!")
 
    # subseting
    sub_ss <- subset(object@summary_sets[[i]]@ss, object@summary_sets[[i]]@ss$rsid %in% rsid_avail)
    
    index <- match(object@summary_sets[[i]]@ss$rsid, rsid_avail)
    ld <- object@ld_matrix[index, index]


    h <- harmonise_ld_dat(sub_ss,ld)
    object@summary_sets[[i]]@ss <- h[[1]]
    

    message("\nThere are ", dim(h[[1]])[1], " variants remaining after harmonising ", object@summary_sets[[i]]@metadata$id, " against the LD matrix.")
    }
  
  object@ld_matrix <- h[[2]]
  object@is_harmonisedLD <- TRUE
  
  message("Done!")
  return(object)
}
)



#' Harmonise data against LD matrix
#' 
#' Function to create a LDmatrix gwasglue2 object and set the @slot ld_matrix using ieugwasr::ld_matrix()
#' @param dataset The DataSet gwasglue2 object
#' @param ld_ref If bfile = TRUE, corresponds to the path and prefix of the plink files used to build the LD correlation matrix. If pop = TRUE, it corresponds to the population code in ieugwasr (eg. EUR) instead (default NULL). 
#' @param pop logical (default FALSE)
#' @param bfile logical (default FALSE)
#' @param plink_bin Path to the plink executable
#' @return  The DataSet gwasglue2 object harmonised
#' @export 
harmonise_ld <- function(dataset, ld_ref = NULL, pop = FALSE, bfile = FALSE, plink_bin = NULL){

    ld <- buildLDMatrix(dataset, ld_ref = ld_ref, pop = pop, bfile = bfile, plink_bin = plink_bin) 
    
    ds <- harmoniseLDMatrix(ld)
    return(ds)
}