
# GWAS summary data
d1 <- dplyr::as_tibble(read.table(system.file("tests", "ieu-a-2_TopHits_sumdata.txt", package="gwasglue2")))
d2 <- dplyr::as_tibble(read.table(system.file("tests", "ieu-a-7_sumdata_ieu-0-7TopHits.txt", package="gwasglue2")))
data <- list(d1,d2)

# get metadata and create metadata objects
m1 <- read.table(system.file("tests", "ieu-a-2_metadata.txt", package="gwasglue2"))
m2 <- read.table(system.file("tests", "ieu-a-7_metadata.txt", package="gwasglue2"))
metadata <- list(create_metadata(m1), create_metadata(m2))

test_that("compare against 2samplemr", { 

# create dataset using create_summaryset() and create_dataset()
  dataset <- lapply(seq_along(data), function(i){
     # create summarysets
    dt <- create_summaryset(data[[i]], metadata=metadata[[i]])
  }) %>%
    # create dataset
    create_dataset(., harmonise = TRUE, tolerance = 0.08, action = 1)
      
  expect_equal(length(dataset@summary_sets), 2)
})


