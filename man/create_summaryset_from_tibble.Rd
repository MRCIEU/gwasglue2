% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructors.R
\name{create_summaryset_from_tibble}
\alias{create_summaryset_from_tibble}
\title{A function to create a gwasglue2 SummarySet object from a tibble}
\usage{
create_summaryset_from_tibble(
  data = tibble(),
  metadata = NULL,
  qc = FALSE,
  beta_col = "beta",
  se_col = "se",
  samplesize_col = "n",
  pvalue_col = "p",
  chr_col = "chr",
  position_col = "position",
  rsid_col = "rsid",
  effect_allele_col = "ea",
  other_allele_col = "nea",
  eaf_col = "eaf",
  id_col = "id",
  trait_col = "trait",
  build = NULL
)
}
\arguments{
\item{data}{GWAS summary statistics. A tibble}

\item{metadata}{A list with metadata information. If NULL, it creates metadata with information retrieved from the dataset}

\item{qc}{Quality control. It checks the @param data and look for problems that can stop gwasglue2 from runing. If TRUE gwasglue will try to solve the problems.  Default is FALSE}

\item{beta_col}{Name of column with effect sizes. The default is \code{"beta"} for @param type \code{"tibble"} and \code{"ES"}for @param type \code{"vcf"}..}

\item{se_col}{Name of column with standard errors. The default is \code{"se"} for @param type \code{"tibble"} and \code{"SE"}for @param type \code{"vcf"}.}

\item{samplesize_col}{Column name for sample size. The default is \code{"n"} for @param type \code{"tibble"} and \code{"SS"}for @param type \code{"vcf"}.}

\item{pvalue_col}{Name of column with p-value. The default is \code{"p"}.}

\item{chr_col}{Column name for chromosome . The default is \code{"chr"} for @param type \code{"tibble"} and \code{"seqnames"}for @param type \code{"vcf"}.}

\item{position_col}{Column name for the position. Together, with @param chr gives the physical coordinates of the variant. The default is \code{"position"} for @param type \code{"tibble"} and \code{"start"}for @param type \code{"vcf"}.}

\item{rsid_col}{Required name of column with variants rs IDs. The default is \code{"rsid"} for @param type \code{"tibble"} and \code{"ID"}for @param type \code{"vcf"}.}

\item{effect_allele_col}{Name of column with effect allele. Must contain only the characters "A", "C", "T" or "G". The default is \code{"ea"} for @param type \code{"tibble"} and \code{"ALT"}for @param type \code{"vcf"}.}

\item{other_allele_col}{Name of column with non effect allele. Must contain only the characters "A", "C", "T" or "G". The default is \verb{"nea} for @param type \code{"tibble"} and \code{"REF"}for @param type \code{"vcf"}.}

\item{eaf_col}{Name of column with effect allele frequency. The default is \code{"eaf"} for @param type \code{"tibble"} and \code{"AF"}for @param type \code{"vcf"}.}

\item{id_col}{GWAS study ID column. The default is \code{"id"}.}

\item{trait_col}{Column name for the column with phenotype name corresponding the the variant. The default is \code{"trait"}}

\item{build}{Reference genome assembly to generate the genomic data. Default is NULL.
\itemize{
\item Options are \code{"NCBI34"}, \code{"NCBI35"}, \code{"NCBI36"}, \code{"GRCh37"} or "GRCh38".
}}
}
\value{
A gwasglue2 SummarySet object
}
\description{
A function to create a gwasglue2 SummarySet object from a tibble
}
\seealso{
\code{\link[=create_metadata]{create_metadata()}} to create a metadata object
}
