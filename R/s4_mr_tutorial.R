
library(tibble)
library(methods)
library(dplyr)

source("ieugwas_utils.R")
source("summaryset.R")
source("dataset.R")
source("convertTo.R")


# create S4 SummarySet objects
# TODO Gib: Here, should we call ids insteadof traits?
x<- clumpTophits(traits = "ieu-a-2")
sumset1 <- SummarySet(traits = "ieu-a-2", variants = x, tools = "mr")
sumset1
sumset2 <- SummarySet(traits="ieu-a-7", variants=x,tools ="mr")
sumset2

# Fill the slots
sumset1 <- setMetadata(sumset1,source = "IEUopenGWAS", traits = "ieu-a-2")
sumset1 <-setMRlabel(sumset1, mr_label = "exposure")
sumset1

getMetadata(sumset1)
getRSID(sumset1)


sumset2 <- setMetadata(sumset2,source = "IEUopenGWAS", traits = "ieu-a-7")
sumset2 <-setMRlabel(sumset2, mr_label = "outcome")
sumset2

getMetadata(sumset2)
getRSID(sumset2)

# create S4 DataSet object
dataset <- DataSet(sumset1,sumset2)
dataset

intersect(getRSID(sumset1),getRSID(sumset2))

# Convert dataset to TwoSampleMR format
dataset_mr <- convertForTwoSampleMR(dataset)
dataset_mr


# Perform the MR analysis (using TwoSampleMR::harmonise_data() function)
mr_result<- TwoSampleMR::harmonise_data(dataset_mr@sumset[[1]]@ss,dataset_mr@sumset[[2]]@ss)  %>%
  TwoSampleMR::mr(., method_list="mr_ivw")

mr_result

