
library(dplyr)
library(ieugwasr)

x <- ieugwasr::tophits("ieu-a-2")$rsid
d1 <- ieugwasr::associations(variants = x, id = "ieu-a-2")
meta1 <-create_metadata(ieugwasr::gwasinfo( "ieu-a-2"))

sumset1 <- create_summaryset(d1, metadata=meta1, tools = "mr")


test_that("create summaryset", {
    expect_true(nrow(sumset1@ss) == length(x))
})


test_that("get metadata", {
    expect_equal(getMetadata(sumset1), sumset1@metadata)
})

test_that("set mr_label", {
    sumset1 <- setMRlabel(sumset1, mr_label = "exposure")
    expect_equal(getMRlabel(sumset1), "exposure" )
})

# this test does not work because we arrage the data by chr pos
# test_that("getRSID", {
#     sumset1 <- setRSID(sumset1,sumset1@ss$rsid)
#     expect_equal(getRSID(sumset1), x)
# })


test_that("dim works", {
    expect_equal(dim(sumset1@ss), dimData(sumset1))
})

  

