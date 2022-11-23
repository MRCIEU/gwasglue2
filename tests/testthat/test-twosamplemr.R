test_that("twosamplemr gives the same as gwasglue2", {
  skip("Not ready yet")
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
  # TODO: convert this to a test

})
