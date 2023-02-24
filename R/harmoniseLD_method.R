

# build LD matrices
setGeneric("buildLDMatrix", function(object, ...) standardGeneric("buildLDMatrix"))
setMethod("buildLDMatrix", "DataSet", function(object, pop = FALSE, bfile = FALSE, plink_bin = NA){
  ldr <- array()
  for (i in seq_along(object@summary_sets)){
    # LD reference dataset name
    ldr[i] <- object@summary_sets[[i]]@ld_ref
      message("Building LD matrix")
      if (pop == TRUE){
        object@ld_matrices[[i]] <- suppressWarnings(ieugwasr::ld_matrix(object@summary_sets[[i]]@variants, pop = ldr[i], with_alleles = TRUE, bfile = FALSE,plink_bin = plink_bin))
      }
      if (bfile == TRUE){
        object@ld_matrices[[i]] <- suppressWarnings(ieugwasr::ld_matrix(object@summary_sets[[i]]@variants, pop=FALSE, with_alleles=TRUE, bfile=ldr[i],plink_bin=plink_bin))
      }
      names(object@ld_matrices)[i] <- ldr[i]
      rsid_avail <- rownames(object@ld_matrices[[i]])
      message("\nData available for ", length(rsid_avail), " variants")

  }
  return(object)
}
)


# Get methods for LDMatrix
setGeneric("getLDMatrix",function(object,...) standardGeneric("getLDMatrix"))
setMethod("getLDMatrix", "DataSet",
          function(object,index) {
            return(object@ld_matrices[[index]])
          })




# Set and get methods for harmonise data against LD matrix
setGeneric("isHarmonisedLD",function(object) standardGeneric("isHarmonisedLD"))
setMethod("isHarmonisedLD","DataSet",function(object) {
  return(object@is_harmonisedLD)
}
)

setGeneric("harmoniseLDMatrix", function(object, args) standardGeneric("harmoniseLDMatrix"))
setMethod( "harmoniseLDMatrix", "DataSet", function(object) {
  message("Gwasglue is now harmonising the data against LD matrices!")
  for (i in seq_along(object@summary_sets)){
    rsid_avail <- do.call(rbind, strsplit(rownames(object@ld_matrices[[i]]), split="_"))[,1]
    # subseting
    sub_ss <- subset(object@summary_sets[[i]]@ss, object@summary_sets[[i]]@ss$rsid %in% rsid_avail)
    
    index <- match(object@summary_sets[[i]]@ss$rsid, rsid_avail)
    ld <- object@ld_matrices[[i]][index, index]


h <- harmonise_ld_dat(sub_ss,ld)
object@summary_sets[[i]]@ss <- h[[1]]
object@ld_matrices[[i]] <- h[[2]]

 message("\nThere are ", dim(h[[1]])[1], " variants remaining after harmonising against a LD matrix.")

object@is_harmonisedLD <- TRUE
 message("Data is harmonised!")

}
return(object)
}
)

