#  this file contains some code used during development

library(devtools)
load_all()



# Checks
devtools::check()
devtools::document()
pkgdown::build_site()



roxygen2::update_collate(".")



library(covr)
report()
usethis::use_github_action("test-coverage")

usethis::use_build_ignore("dev/", escape = TRUE)
usethis::use_build_ignore("docker/", escape = TRUE)
usethis::use_build_ignore(".dockerignore/", escape = TRUE)
usethis::use_build_ignore(".tests/", escape = TRUE)
usethis::use_build_ignore(".devcontainer/", escape = TRUE)
usethis::use_build_ignore("data/", escape = TRUE)

usethis::use_git_ignore("data/")


usethis::use_news_md()


# ********************************************************
# Test harmonise functions
# ********************************************************
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

# test harmonise()
dat1 <- ieugwasr::tophits("ieu-a-89")
dat2 <- deharmonise(dat1)
table(dat2$flipped, dat2$switched)
table(dat2$palindromic)
table(dat2$drop)
rsid=dat1$rsid
A1 <- dat1$ea
A2 <- dat1$nea
B1 <- dat2$ea
B2 <- dat2$nea
betaA <- dat1$beta
betaB <- dat2$beta
fA <- dat1$eaf
fB <- dat2$eaf

tolerance=0.08
action=2
h<- harmonise(rsid, A1, A2, B1, B2, betaA, betaB, fA, fB, tolerance, action)
dim(h)
dim(dat1)

####################
# test  harmoniseData()
x<- clumpTophits(traits = "ieu-a-89")
sumset1 <- SummarySet(traits = "ieu-a-89", variants = x, tools = "mr")
sumset1

#randomly reduce number of snps for deharmonised dataset : test different number of snps
x1 <- x[-round(runif(10,1,length(x)))]
sumset2 <- SummarySet(traits = "ieu-a-89", variants = x1, tools = "mr")
sumset2@ss <-deharmonise(sumset2@ss)

sumset2
table(sumset2@ss$flipped, sumset2@ss$switched)
table(sumset2@ss$palindromic)
table(sumset2@ss$drop)


#randomly reduce number of snps for deharmonised dataset2 : test different number of snps and more than 2 sumset
x3 <- x[-round(runif(20,1,length(x)))]
sumset3 <- SummarySet(traits = "ieu-a-89", variants = x3, tools = "mr")
sumset3@ss <-deharmonise(sumset3@ss)

sumset3
table(sumset3@ss$flipped, sumset3@ss$switched)
table(sumset3@ss$palindromic)
table(sumset3@ss$drop)

# Fill the slots
sumset1 <- setMetadata(sumset1,source = "IEUopenGWAS", traits = "ieu-a-89")
sumset1 <-setMRlabel(sumset1, mr_label = "exposure")
getMetadata(sumset1)
getRSID(sumset1)


sumset2 <- setMetadata(sumset2,source = "IEUopenGWAS", traits = "ieu-a-89")
sumset2 <-setMRlabel(sumset2, mr_label = "outcome")
getMetadata(sumset2)
getRSID(sumset2)


sumset3 <- setMetadata(sumset3,source = "IEUopenGWAS", traits = "ieu-a-89")
sumset3 <-setMRlabel(sumset3, mr_label = "exposure")
getMetadata(sumset3)
getRSID(sumset3)


# create S4 DataSet object
dataset <- DataSet(sumset1,sumset2,sumset3) %>% 
	overlapSNP(object) %>% 
  	harmoniseData(.,tolerance = 0.08,action = 2)
# ********************************************************