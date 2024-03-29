% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dataset.R
\docType{methods}
\name{getData}
\alias{getData}
\alias{getData,DataSet-method}
\title{Get Method to retrieve the GWAS Summary Statistics}
\usage{
getData(dataset, index)

\S4method{getData}{DataSet}(dataset, index)
}
\arguments{
\item{dataset}{A gwasglue2 DataSet object}

\item{index}{Index of gwasglue2 SummarySet objects within DataSet}
}
\value{
A tibble with GWAS summary statistics
}
\description{
Get Method to retrieve the GWAS Summary Statistics
}
\seealso{
Similar to \code{\link[=getSummaryData]{getSummaryData()}}
}
