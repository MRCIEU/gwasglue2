

#' An S4 class to represent the ListDataSet object
############################
#'
#' @slot datasets A list of DataSet objects (default NA).
#' @export 
#' @rdname ListDataSets
setClass("ListDataSets",
  slots = c(
    datasets = "list"
    ),
  prototype = prototype(
    datasets = list(NA_character_)
  ),
  contains = c(class(dplyr::tibble()))
)


#' ListDataSets function
#'
#' @param ... Array of gwasglue2 DataSett object names.
#' @importFrom methods new
#' @return  A gwasglue2 ListDataSet object
#' @export 
#' @rdname ListDataSets
ListDataSets <- function(...) {
  new("ListDataSets", datasets = as.list(...))
}



#' Get Method to retrieve the gwasglue2 DataSet object
#'
#' @param list_datasets A gwasglue2 ListDataSet object
#' @param index Index of gwasglue2 SummarySet objects within DataSet
#' @return A gwasglue2 DataSet object
#' @export
#' @docType methods
#' @rdname getDataSet-methods
setGeneric("getDataSet", function(list_datasets,index) standardGeneric("getDataSet"))

#' @rdname getDataSet-methods
setMethod("getDataSet", "ListDataSets",
          function(list_datasets,index) {
            return(list_datasets@datasets[[index]])
          })



# show method 
setMethod(f = "show", signature="ListDataSets", definition = function(object) {
  # set description of DataSet
  length <- getLength(object)

  
  # write
  cat("A 'ListDataSets' object with", length, "DataSets.\n")
 
  cat("\nTo access the gwasglue2 `DataSets` in each of the 'ListDataSets' use getDataSet(list_datasets, index).\n")

})
