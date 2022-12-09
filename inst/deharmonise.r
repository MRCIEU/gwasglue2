# This function takes an output from ieugwasr lookup (dat)
# And randomly decides whether to
# - flip strand
# - switch effect allele
# - both
# - neither
# the switched and flipped columns indicate this
# if an allele is flipped, palindromic and intermediate allele frequency then it would be dropped
deharmonise <- function(dat, tolerance=0.08) {
  flip_allele <- function(dat, i) {
    if(dat$ea[i] == "A") dat$ea[i] == "T"
    if(dat$ea[i] == "T") dat$ea[i] == "A"
    if(dat$ea[i] == "C") dat$ea[i] == "G"
    if(dat$ea[i] == "G") dat$ea[i] == "C"
    if(dat$nea[i] == "A") dat$nea[i] == "T"
    if(dat$nea[i] == "T") dat$nea[i] == "A"
    if(dat$nea[i] == "C") dat$nea[i] == "G"
    if(dat$nea[i] == "G") dat$nea[i] == "C"
    return(dat[i,])
  }
  switch_allele <- function(dat, i) {
    temp <- dat$ea[i]
    dat$ea[i] <- dat$nea[i]
    dat$nea[i] <- temp
    dat$beta[i] <- dat$beta[i] * -1
    dat$eaf[i] <- 1 - dat$eaf[i]
    return(dat[i,])
  }
  # choose SNPs to flip
  dat$flipped <- dat$switched <- FALSE
  for(i in 1:nrow(dat)) {
    if(sample(1:2,1) == 1) dat[i,] <- flip_allele(dat, i) %>% mutate(flipped=TRUE)
    if(sample(1:2,1) == 2) dat[i,] <- switch_allele(dat, i) %>% mutate(switched=TRUE)
  }
  # identify SNPs to drop
  dat$palindromic <- paste(dat$ea, dat$nea) %in% c("A T", "T A", "G C", "C G")
  dat$drop <- dat$flipped & dat$palindromic & dat$eaf > (0.5-tolerance) & dat$eaf < (0.5+tolerance)
  return(dat)
}





# Example
dat1 <- ieugwasr::tophits("ieu-a-89")
dat2 <- deharmonise(dat1)
table(dat2$flipped, dat2$switched)
table(dat2$palindromic)
table(dat2$drop)

# need to harmonise dat2 against dat
# should get back the same thing except for instances when palindromic SNPs

