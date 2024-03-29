% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructors.R
\name{add_summaryset}
\alias{add_summaryset}
\title{Add a SummarySet to a DataSet}
\usage{
add_summaryset(
  summary_sets,
  dataset,
  harmonise = TRUE,
  tolerance = 0.08,
  action = 1
)
}
\arguments{
\item{summary_sets}{one or more gwasglue2 Summarysets objects to add to an existent DataSet object. If more than one it should be a list}

\item{dataset}{The gwasglue2 DataSet object to add to}

\item{harmonise}{logical (default TRUE). It harmonises the summary sets in the DataSet against each other.}

\item{tolerance}{Inherited from harmoniseData() (default 0.08)}

\item{action}{Inherited from harmoniseData() (Default 1)
\itemize{
\item \code{action = 1}: Assume all alleles are coded on the forward strand, i.e. do not attempt to flip alleles
\item \code{action = 2}: Try to infer positive strand alleles, using allele frequencies for palindromes (default, conservative);
\item \code{action = 3}: Correct strand for non-palindromic SNPs, and drop all palindromic SNPs from the analysis (more conservative).
}}
}
\value{
A harmonised gwasglue2 DataSet object with input SummarySets added
}
\description{
Add a SummarySet to a DataSet
}
