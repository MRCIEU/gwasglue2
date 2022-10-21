############################
# DataSet class

#' An S4 class to represent the Data Set
#'
#' @slot sumset A list of SummarySet objects (default NA).
#' @slot ld_matrices A list of LD correlation matrices (default NA).
#' @slot is_harmonised logical (default FALSE).
#' @slot is_harmonisedLD logical (default FALSE).
#' @slot is_converted logical (default FALSE).
setClass("DataSet",
  slots = c(
    sumset = "list",
    ld_matrices = "list",
    is_harmonised = "logical",
    is_harmonisedLD = "logical",
    is_converted = "logical"
  ),
  prototype = prototype(
    sumset = list(NA_character_),
    ld_matrices = list(NA_character_),
    is_harmonised = FALSE,
    is_harmonisedLD = FALSE,
    is_converted = FALSE
  ),
  contains = c(class(tibble()))
)


#' DataSet function
#'
#' @param ...
#'
#' @return
DataSet <- function(...) {
  new("DataSet", sumset = list(...))
}

