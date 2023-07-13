# unit testing for SummarySet object

# GWAS summary data
d1 <- dplyr::as_tibble(read.table(system.file("tests", "ieu-a-2_TopHits_sumdata.txt", package="gwasglue2")))
m1 <- read.table(system.file("tests", "ieu-a-2_metadata.txt", package="gwasglue2"))
meta1 <-create_metadata(m1)

sumset1 <- create_summaryset(d1, metadata = meta1)  %>%
        suppressMessages()


test_that("create summaryset", {
    expect_true(nrow(sumset1@ss) == dim(d1)[1])
})


test_that("get metadata", {
    expect_equal(getMetadata(sumset1), sumset1@metadata)
})

test_that("set mr_label", {
    sumset1 <- setAttributes(sumset1, mr_label = "exposure")
    expect_equal(getAttributes(sumset1)$mr_label, "exposure" )
})

# this test does not work because we arrage the data by chr pos
# test_that("getRSID", {
#     sumset1 <- setRSID(sumset1,sumset1@ss$rsid)
#     expect_equal(getRSID(sumset1), x)
# })


test_that("dim works", {
    expect_equal(dim(sumset1@ss), dimData(sumset1))
})

  

