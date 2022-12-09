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
sumset1

getMetadata(sumset1)
getRSID(sumset1)


sumset2 <- setMetadata(sumset2,source = "IEUopenGWAS", traits = "ieu-a-89")
sumset2 <-setMRlabel(sumset2, mr_label = "outcome")
sumset2

getMetadata(sumset2)
getRSID(sumset2)


sumset3 <- setMetadata(sumset3,source = "IEUopenGWAS", traits = "ieu-a-89")
sumset3 <-setMRlabel(sumset3, mr_label = "exposure")
sumset3

getMetadata(sumset3)
getRSID(sumset3)


# create S4 DataSet object
dataset <- DataSet(sumset1,sumset2,sumset3)
dataset


object = dataset
source("../R/harmonise.R")
source("../R/harmonise_method.R")
object1 <- overlapSNP(object) %>%
  
  harmoniseData(.,tolerance = 0.08,action = 2)
