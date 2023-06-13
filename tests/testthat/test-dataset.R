library(dplyr)
library(ieugwasr)


test_that("compare against 2samplemr", { 

  # lookk for the tophit variants
  x <- ieugwasr::tophits("ieu-a-2")$rsid 
  ids <- c("ieu-a-2",  "ieu-a-7")

  # get metadata and create metadata objects
  metadata <- lapply(seq_along(ids), function(i){
    m <- create_metadata(ieugwasr::gwasinfo(ids[i])) 
  })

# create dataset using create_summaryset() and create_dataset()
  dataset <- lapply(seq_along(ids), function(i){
     # create summarysets
    dt <- create_summaryset(ieugwasr::associations(variants = x, id =ids[i]), metadata=metadata[i])
  }) %>%
    # create dataset
    create_dataset(., harmonise = TRUE, tolerance = 0.08, action = 1)
      
  expect_equal(length(dataset@summary_sets), 2)
})


