# WIP

# Set and get methods for harmonise data
setGeneric("isHarmonised",function(object) standardGeneric("isHarmonised"))
setMethod("isHarmonised","DataSet",function(object) {
  return(object@is_harmonised)
}
)

setGeneric("harmonise", function(object, args) standardGeneric("harmonise"))
setMethod( "harmonise", "DataSet", function(object) {
  message("Gwasglue is now harmonising the data!")
  for (i in seq_along(object@sumset)){

  }




  object@is_harmonised <- TRUE
  return(object)
}
)



path_to_plink <- "../poc/plink"
# build LD matrices
setGeneric("buildLDMatrix", function(object, args) standardGeneric("buildLDMatrix"))
setMethod( "buildLDMatrix", "DataSet", function(object, pop, bfile, plink_bin) {
  message("Building LD matrix!")
  ld_ref <- array()
  count <- 0
  ldr[1] <- object@sumset[[1]]@ld_ref

  for (i in seq_along(object@sumset)){
    # LD reference dataset name

    if(i == 1 || object@sumset[[i]]@ld_ref != ldr[count]) {
      ldr[count + 1] <- object@sumset[[i]]@ld_ref
      count <- count + 1
      ldr <- ld_ref[count]
      message("Building LD matrix")
      if (pop == TRUE) {
        object@ld_matrices[count] <- suppressWarnings(ieugwasr::ld_matrix(object@sumset[[i]]@variants, pop = ldr, with_alleles = TRUE, bfile = FALSE,plink_bin = plink_bin))
      }
      if (bfile == TRUE){
        object@ld_matrices[count] <- suppressWarnings(ieugwasr::ld_matrix(object@sumset[[i]]@variants, pop=FALSE, with_alleles=TRUE, bfile=ldr,plink_bin=plink_bin))
      }
      names(object@ld_matrices)[count] <- ldr
      rsid_avail <- rownames(object@ld_matrices)
      message("Data available for ", length(rsid_avail), " variants")
    }
    else {
      ldr <- ld_ref[count]
    }

  }

  return(object)
}
)




# Set and get methods for harmonise data against LD matrix
setGeneric("isHarmonisedLD",function(object) standardGeneric("isHarmonisedLD"))
setMethod("isHarmonisedLD","DataSet",function(object) {
  return(object@is_harmonisedLD)
}
)

setGeneric("harmoniseLDMatrix", function(object, args) standardGeneric("harmoniseLDMatrix"))
setMethod( "harmoniseLDMatrix", "DataSet", function(object) {
  message("Gwasglue is now harmonising the data against LD matrices!")
  for (i in seq_along(object@sumset)){

  }
  object@sumset[[i]]@ld_ref



  object@is_harmonisedLD <- TRUE
  return(object)
}
)




