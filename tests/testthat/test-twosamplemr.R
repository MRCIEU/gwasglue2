




test_that("twosamplemr gives the same as gwasglue2", {

x <- ieugwasr::tophits("ieu-a-2")$rsid
d1 <- ieugwasr::associations(variants = x, id = "ieu-a-2")
d2 <- ieugwasr::associations(variants = x, id = "ieu-a-7")
 
meta1 <-create_metadata(ieugwasr::gwasinfo( "ieu-a-2"))
meta2 <-create_metadata(ieugwasr::gwasinfo( "ieu-a-7"))
meta <- list(meta1,meta2)


  ##############################################
  #Using TwoSampleMR harmonising function
  # create S4 DataSet object
  sumset1 <- create_summaryset(d1, metadata=meta1, tools = "mr")
  sumset2 <- create_summaryset(d2, metadata=meta2, tools ="mr")
  sumset1 <-setMRlabel(sumset1, mr_label = "exposure")
  sumset2 <-setMRlabel(sumset2, mr_label = "outcome")
  dataset <- DataSet(sumset1,sumset2) %>%
  overlapVariants(., strand = "forward") %>%
  harmoniseData(.,tolerance = 0.08, action = 1)


  # Convert dataset to TwoSampleMR format
  dataset_mr <- convertForTwoSampleMR(dataset)


  # Perform the MR analysis 
  mr_result1 <-  merge(getData(dataset_mr,1),getData(dataset_mr,2), by = c("SNP", "mr_keep"))  %>%
    TwoSampleMR::mr(., method_list="mr_ivw") %>% 
    select(id.exposure, id.outcome, method, nsnp, b, se, pval)

  # TODO: convert this to a test
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
