

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
#' @rdname geDataSet-methods
setGeneric("getDataSet", function(dataset,index) standardGeneric("getDataSet"))

#' @rdname geDataSet-methods
setMethod("getDataSet", "ListDataSets",
          function(dataset,index) {
            return(dataset@datasets[[index]])
          })
