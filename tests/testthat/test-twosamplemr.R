# Unit testing for MR analyses
library(dplyr)


d1 <- as_tibble(read.table(system.file("tests", "ieu-a-2_TopHits_sumdata.txt", package="gwasglue2")))
d2 <- as_tibble(read.table(system.file("tests", "ieu-a-7_sumdata_ieu-0-7TopHits.txt", package="gwasglue2")))
data <- list(d1,d2)

# get metadata and create metadata objects
m1 <- read.table(system.file("tests", "ieu-a-2_metadata.txt", package="gwasglue2"))
m2 <- read.table(system.file("tests", "ieu-a-7_metadata.txt", package="gwasglue2"))
metadata <- list(create_metadata(m1), create_metadata(m2))

# mr labels
mr_labels <- c("exposure","outcome")

test_that("twosamplemr gives the same as gwasglue2", {

  #  create dataset and convert it to mr format
  dataset <- lapply(seq_along(data), function(i){
     # create summarysets
    dt <- create_summaryset(data[[i]], metadata=metadata[[i]]) %>% 
      setAttributes(., mr_label = mr_labels[i])
  }) %>%
    # create dataset
    create_dataset(., harmonise = TRUE, tolerance = 0.08, action = 1)  %>%
    # Convert dataset to TwoSampleMR format
   convertForTwoSampleMR(.)


  # Perform the MR analysis 
  mr_result1 <-  merge(getData(dataset,1),getData(dataset,2), by = c("SNP", "mr_keep"))  %>%
    TwoSampleMR::mr(., method_list="mr_ivw") %>% 
    select(id.exposure, id.outcome, method, nsnp, b, se, pval)

 mr_result2 <- TwoSampleMR::make_dat("ieu-a-2","ieu-a-7") %>%
    TwoSampleMR::mr(., method_list="mr_ivw") %>% 
    select(id.exposure, id.outcome, method, nsnp, b, se, pval)

  expect_equal(mr_result1, mr_result2)

})


# # compare with TwoSampleMR

# # TODO we need to standardise TWOsampleMR results (alleles by alphabetical order)


# skip("Not ready yet")
# test_that("simple 2 trait harmonisation matches TwoSampleMR", {

#   dataset <- create_dataset(data, metadata = meta, tools = c("mr"), harmonise = TRUE, tolerance = 0.08, action = 1)
#   # Do the same in TwoSampleMR
#   dat <- TwoSampleMR::make_dat("ieu-a-2","ieu-a-7")

#   expect_equal(
#     cor(dataset@summary_sets[[1]]@ss$beta, dataset@summary_sets[[2]]@ss$beta),
#     cor(dat$beta.exposure, dat$beta.outcome)
#   )
# })
