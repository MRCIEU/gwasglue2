context("creating a dataset")
library(dplyr)

test_that("", {
  dataset <- DataSet(sumset1,sumset2) %>%
    overlapSNP(.) %>%
    harmoniseData(.,tolerance = 0.08,action = 2)
})
