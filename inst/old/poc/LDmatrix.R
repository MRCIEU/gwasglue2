# LDmatrix class

#' An S4 class to represent the LD matrix
############################
#'
#' @slot LDmatrix 
#' @export 
setClass("LDmatrix",
  slots = c(
    ld_matrix = "matrix"
   ),
  prototype = prototype(
    ld_matrix = matrix(NA_real_)
  )

)


#' LDmatrix function
#'
#' @param ... Method to fill the @slot ld_matrix (e.g. use ieugwasr::ld_matrix())
#' @importFrom methods new
#' @return  A LDmatrix gwasglue2 object
#' @export 
LDmatrix <- function(...) {
  new("LDmatrix", 
   ld_matrix = ...)
}

# build LD matrices
setGeneric("buildLDMatrix", function(ld, dataset, ...) standardGeneric("buildLDMatrix"))
setMethod("buildLDMatrix", c("LDmatrix"), function(ld, dataset, ld_ref = NULL, pop = FALSE, bfile = FALSE, plink_bin = NULL){
  

      message("Building LD matrix")
      if (pop == TRUE){
        ld@ld_matrix <- suppressWarnings(ieugwasr::ld_matrix(dataset@summary_sets[[1]]@variants, pop = ld_ref, with_alleles = TRUE, bfile = FALSE,plink_bin = plink_bin))
      }
      if (bfile == TRUE){
        ld@ld_matrix <- suppressWarnings(ieugwasr::ld_matrix(dataset@summary_sets[[1]]@variants, pop=FALSE, with_alleles=TRUE, bfile=ld_ref, plink_bin=plink_bin))
      }

      rsid_avail <- rownames(ld@ld_matrix)
      message("\nData available for ", length(rsid_avail), " variants")

  return(ld)
}
)






# Get methods for LDMatrix
setGeneric("getLDMatrix",function(ld) standardGeneric("getLDMatrix"))
setMethod("getLDMatrix", "LDmatrix",
          function(ld) {
            return(ld@ld_matrix)
          })





# Set and get methods for harmonise data against LD matrix
setGeneric("isHarmonisedLD",function(dataset) standardGeneric("isHarmonisedLD"))
setMethod("isHarmonisedLD","DataSet",function(dataset) {
  return(dataset@is_harmonisedLD)
}
)

setGeneric("harmoniseLDMatrix", function(dataset,ld) standardGeneric("harmoniseLDMatrix"))
setMethod( "harmoniseLDMatrix", "DataSet", function(dataset,ld) {
  message("Gwasglue is now harmonising the data against LD matrix!")
  for (i in seq_along(dataset@summary_sets)){
    rsid_avail <- do.call(rbind, strsplit(rownames(ld@ld_matrix), split="_"))[,1]
    # subseting
    sub_ss <- subset(dataset@summary_sets[[i]]@ss, dataset@summary_sets[[i]]@ss$rsid %in% rsid_avail)
    
    index <- match(dataset@summary_sets[[i]]@ss$rsid, rsid_avail)
    ldm <- ld@ld_matrix[index, index]


h <- harmonise_ld_dat(sub_ss,ldm)
dataset@summary_sets[[i]]@ss <- h[[1]]
ld@ld_matrix <- h[[2]]

 message("\nThere are ", dim(h[[1]])[1], " variants remaining after harmonising against a LD matrix.")

dataset@is_harmonisedLD <- TRUE
 message("Data is harmonised!")

}
return(dataset)
}
)


create_ld_matrix  <- function(dataset, ld_ref = NULL, pop = FALSE, bfile = FALSE, plink_bin = NULL){
  
    # LD reference dataset name
    ld <- LDmatrix()
      message("Building LD matrix")
      if (pop == TRUE){
        ld@ld_matrix <- suppressWarnings(ieugwasr::ld_matrix(dataset@summary_sets[[1]]@variants, pop = ld_ref, with_alleles = TRUE, bfile = FALSE,plink_bin = plink_bin))
      }
      if (bfile == TRUE){
        ld@ld_matrix <- suppressWarnings(ieugwasr::ld_matrix(dataset@summary_sets[[1]]@variants, pop=FALSE, with_alleles=TRUE, bfile=ld_ref,plink_bin=plink_bin))
      }

      rsid_avail <- rownames(ld@ld_matrix)
      message("\nData available for ", length(rsid_avail), " variants")

  return(ld)
}


 
 harmonise_ld <- function(dataset, ld_ref = NULL, pop = FALSE, bfile = FALSE, plink_bin = NULL){

  ld <- create_ld_matrix(dataset = dataset, ld_ref = ld_ref, pop = pop, bfile = bfile, plink_bin = plink_bin)

message("Gwasglue is now harmonising the data against LD matrix!")
  for (i in seq_along(dataset@summary_sets)){
    rsid_avail <- do.call(rbind, strsplit(rownames(ld@ld_matrix), split="_"))[,1]
    # subseting
    sub_ss <- subset(dataset@summary_sets[[i]]@ss, dataset@summary_sets[[i]]@ss$rsid %in% rsid_avail)
    
    index <- match(dataset@summary_sets[[i]]@ss$rsid, rsid_avail)
    ldm <- ld@ld_matrix[index, index]


h <- harmonise_ld_dat(sub_ss,ldm)
dataset@summary_sets[[i]]@ss <- h[[1]]
ld@ld_matrix <- h[[2]]

 message("\nThere are ", dim(h[[1]])[1], " variants remaining after harmonising against a LD matrix.")
 }
dataset@is_harmonisedLD <- TRUE
 message("Data is harmonised!")

 
  return(dataset)
}
 

# Set and get methods for ld_ref
# setGeneric("setLDref",function(object,ld_ref) standardGeneric("setLDref"))
# setMethod( "setLDref", "SummarySet",
#            function(object,ld_ref) {
#              object@ld_ref <- ld_ref
#              message(paste("Gwasglue is going to harmonise data against ", object@ld_ref, "LD correlation matrix"))
#              return(object)
#            }
# )

# setGeneric("getLDref",function(object) standardGeneric("getLDref"))
# setMethod("getLDref","SummarySet",
#           function(object) {
#             return(object@ld_ref)
#           }
# )
