% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/harmonise_method.R
\docType{methods}
\name{overlapVariants}
\alias{overlapVariants}
\alias{overlapVariants,DataSet-method}
\title{Look for overlapped variants between SummarySets in the DataSet and Resize}
\usage{
overlapVariants(dataset, action)

\S4method{overlapVariants}{DataSet}(dataset, action)
}
\arguments{
\item{dataset}{The gwasglue2 DataSet object}

\item{action}{Level of strictness in dealing with SNPs during harmonisation.
\itemize{
\item \code{action = 1}: Assume all alleles are coded on the forward strand, i.e. do not attempt to flip alleles
\item \code{action = 2}: Try to infer positive strand alleles, using allele frequencies for palindromes (default, conservative);
\item \code{action = 3}: Correct strand for non-palindromic SNPs, and drop all palindromic SNPs from the analysis (more conservative).
}}
}
\value{
The gwasglue2 DataSet object resized
}
\description{
Look for overlapped variants between SummarySets in the DataSet and Resize
}
