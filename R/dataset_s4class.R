# Functions

clumped_hits <- function(traits,  clump = TRUE, source = ieugwasr::check_access_token()) {
  # ch  <- ieugwasr::tophits(traits[1], clump=clump, source=source)

  #use this line code instead if there is an "unused argument (source = source)" error. It happens if the access was already opened
  ch  <- ieugwasr::tophits(traits[1], clump = clump)
  snps <- ch$rsid
  return(snps)
}

create_dataset_1 <- function(trait, variants) {
  data_set<-  ieugwasr::associations(variants = variants, id = trait)
  return(data_set)
}


x<-clumped_hits(traits="ieu-a-2")

library(tibble)
library(methods)

# First, the list of all member data and their types


setOldClass("data.frame")
setClassUnion("tibbleORlist",c("tbl_df","tbl","dat.frame"))

#setOldClass("tbl","tbl_df")
setClass("Dataset_gwas",
         slots = c(
           trait = "character",
           variants = "character",
           dt_tibble = "tbl_df"
           #,
           #source = "character",
           #factors = "numeric"
          #,annotation = "list"
          ,exposure = "character"
          , methods="vector"

          )
         ,
         contains = class(tibble())
         #contains = "tibbleORlist"
         )

setMethod(
  f = "initialize", "Dataset_gwas",
  function(.Object, ...) {
    .Object <- callNextMethod()
    .Object@dt_tibble <- tibble()
    .Object
  }
)

# Next, the list of all members functions and their signatures
setGeneric("createDataset",function(object) standardGeneric("createDataset"))

# Next, the implementation of all of the member functions
setMethod("createDataset", "Dataset_gwas", function(object){

  object@dt_tibble <- create_dataset_1(trait=object@trait,variants=object@variants)
  return(object)

})


# Finally, we define the constructor of the class - this will likely
# have millions of arguments (eventually)
Dataset <- function(trait,
                    variants,
                    exposure,
                    methods) {
  new("Dataset_gwas", trait=trait, variants=variants,exposure=exposure,methods=methods)
}

dataset <- Dataset(trait="ieu-a-2", variants=x,exposure="ieu-a-2",methods="mr")

str(dataset)
isS4(dataset)

showMethods(class = "Dataset_gwas")


# more than one dataset

setClass("Datasets",
         contains = "Dataset"
         slots = c(
           datasets = "list(DataSet)"
           is_harmonised = "numeric",
         )
)



# Next, the list of all members functions and their signatures
setGeneric("harmonise", function(object) standardGeneric("harmonise"))
setGeneric("harmoniseAgainstLDMatrix", function(object, args) standardGeneric("harmoniseAgainstLDMatrix"))
setGeneric("convertForTwoSampleMR", function(object) standardGeneric("convertForTwoSampleMR"))


