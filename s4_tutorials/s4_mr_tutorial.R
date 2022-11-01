
library(tibble)
library(methods)
library(dplyr)

source("../R/ieugwas_utils.R")
source("../R/summaryset.R")
source("../R/dataset.R")
source("../R/harmonise.R")
source("../R/harmonise_method.R")
source("../R/convertTo.R")


# create S4 SummarySet objects
# TODO Gib: Here, should we call ids instead of traits?
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




# Harmonising inside gwasglue

# create S4 DataSet object

dataset <- DataSet(sumset1,sumset2) %>%
  overlapSNP(.) %>%
  harmoniseData(.,tolerance = 0.08,action = 2)

# Convert dataset to TwoSampleMR format
dataset_mr <- convertForTwoSampleMR(dataset)
dataset_mr

# Perform the MR analysis
mr_result<- merge(dataset_mr@summary_sets[[1]]@ss,dataset_mr@summary_sets[[2]]@ss, by = c("SNP", "mr_keep"))  %>%
  TwoSampleMR::mr(., method_list="mr_ivw")

mr_result


##############################################
#Using TwoSampleMR harmonising function
# create S4 DataSet object
dataset <- DataSet(sumset1,sumset2)
dataset

# Convert dataset to TwoSampleMR format
dataset_mr <- convertForTwoSampleMR(dataset)
dataset_mr


# Perform the MR analysis (using TwoSampleMR::harmonise_data() function)
mr_result<- TwoSampleMR::harmonise_data(dataset_mr@sumset[[1]]@ss,dataset_mr@sumset[[2]]@ss)  %>%
  TwoSampleMR::mr(., method_list="mr_ivw")

mr_result
