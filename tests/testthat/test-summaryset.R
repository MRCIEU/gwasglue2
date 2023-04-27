
library(dplyr)
library(ieugwasr)

x <- ieugwasr::tophits("ieu-a-2")$rsid
d1 <- ieugwasr::associations(variants = x, id = "ieu-a-2")

sumset1 <- SummarySet(d1)

test_that("create summaryset", {
    expect_true(nrow(sumset1@ss) == length(x))
})

test_that("set metadata", {
    sumset1 <- setMetadata(sumset1,source = "IEUopenGWAS", id = "ieu-a-2")
    expect_true(is.list(sumset1@metadata))
})


test_that("get metadata", {
    expect_equal(getMetadata(sumset1), sumset1@metadata)
})

test_that("set mr_label", {
    sumset1 <- setMRlabel(sumset1, mr_label = "exposure")
    expect_equal(getMRlabel(sumset1), "exposure" )
})


test_that("getRSID", {
    sumset1 <- setRSID(sumset1,sumset1@ss$rsid)
    expect_equal(getRSID(sumset1), x)
})


test_that("dim works", {
    expect_equal(dim(sumset1@ss), dimData(sumset1))
})

  
