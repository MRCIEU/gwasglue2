library(dplyr)
library(ieugwasr)

x <- ieugwasr::tophits("ieu-a-2")$rsid
d1 <- ieugwasr::associations(variants = x, id = "ieu-a-2")
d2 <- ieugwasr::associations(variants = x, id = "ieu-a-7")
data <-list (d1,d2)
 
meta1 <-create_metadata(ieugwasr::gwasinfo( "ieu-a-2"))
meta2 <-create_metadata(ieugwasr::gwasinfo( "ieu-a-7"))
meta <- list(meta1,meta2)





test_that("compare against 2samplemr", { 
dataset <- create_dataset(data, metadata = meta, tools = c("mr"), harmonise = TRUE, tolerance = 0.08, action = 1)
  dataset
  expect_equal(length(dataset@summary_sets), 2)
})


