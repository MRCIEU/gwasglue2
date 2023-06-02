
#  Fixed-effects meta-analysis
#  @param beta 
# @param se
# importFrom("stats", "pchisq")

meta.F <- function(beta, se){
  #returns inverse-variance weighted meta-analysis estimate, SE and P-value.
  beta.F = sum(beta / se^2) / sum(1 / se^2)
  se.F = 1 / sqrt(sum(1 / se^2))
  p.F = pchisq( (beta.F / se.F)^2, df = 1, lower.tail = F)
  return(c(beta = beta.F, se = se.F, p = p.F))
}


meta_analyse <- function(dataset, method = "fixed"){
length_dt <- getLength(dataset)
nsnps <- dim(getData(dataset,1))[1]
beta <- matrix(ncol = length_dt, nrow = nsnps)
se <- matrix(ncol = length_dt, nrow = nsnps)

for (i in 1:length_dt){
beta[,i] <- getData(dataset,i)$beta
se[,i] <- getData(dataset,i)$se
}

if(method == "fixed"){
  meta <- matrix(ncol = 3, nrow = nsnps)
  meta <- sapply(1:nsnps, function(i) {
  meta[i,] <- meta.F(beta[i,],se[i,])})
}

return(as_tibble(t(meta)))
}

#' Meta analysis
#' @description Statistical combination of the
#' results from two or more separate studies. 
#' It uses the fixed-effect model assuming that one true 
#' effect size underlies all the studies in the meta-analysis. 
#' @param dataset gwasglue2 DataSet object
#' @importFrom stats pchisq
#' @return gwasglue2 SummarySet object

create_meta_summaryset <- function(dataset) {

  length_dt <- getLength(dataset)
  nsnps <- dim(getData(dataset,1))[1]

  lapply(1:length_dt, function(i) {
    if (is.null(getMetadata(getSummarySet(dataset,i))$sample_size) || is.na(getMetadata(getSummarySet(dataset,i))$sample_size)){
      stop("No sample size information in metadata for at least one of the SummarySets. More details on how to add to metadata in 'help(gwasglue2::create_metadata)' and 'help(gwasglue2::getMetadata)'." )
    }
    })


  # meta data id
  ids <- unlist(lapply(1:length_dt, function(i) {
    id <- getMetadata(getSummarySet(dataset,i))$id }))
  id <- paste(ids, collapse = "||")

  # meta data sample size
  n <- lapply(1:length_dt, function(i) {
    n <- getMetadata(getSummarySet(dataset,i))$sample_size})
  n_meta <- sum(unlist(n))

  # meta data trait
  traits <- unlist(lapply(1:length_dt, function(i) {
    t <- getMetadata(getSummarySet(dataset,i))$trait}))
  trait <- paste(unique(trait), collapse = "||")

  # eaf weighted mean
  eaf <- lapply(1:length_dt, function(i) {
    m <- getData(dataset,i)$eaf * getMetadata(getSummarySet(dataset,i))$sample_size
    })
  eaf_mean <- Reduce("+", eaf) / n_meta

  # create metadata object
  metadata <- create_metadata(metadata = NULL,
                            id = id,
                            sample_size = n_meta,
                            nsnp = nsnps,
                            trait = trait,
                            sd = NA,
                            unit = NA,
                            ncontrol = NA, 
                            build = NA,
                            population = NA,
                            ncase = NA,
                            meta_analysis = TRUE
                            )
  # create the summary statistics tibble
  meta_dt <- meta_analyse(dataset) %>% mutate(n = n_meta, 
                                    chr = getData(dataset,1)$chr, 
                                    position = getData(dataset,1)$position,
                                    rsid = getData(dataset,1)$rsid,
                                    ea = getData(dataset,1)$ea, 
                                    nea = getData(dataset,1)$nea,
                                    eaf = eaf_mean, 
                                    id = id,
                                    trait = trait)
  # create the SummarySet gwasglue2 object
  s <- create_summaryset(data = meta_dt, metadata = metadata)
  return(s)
}


