skip("problems with server")
library(dplyr)
library(ieugwasr)

x <- ieugwasr::tophits("ieu-a-2")$rsid
sumset1 <- SummarySet(traits = "ieu-a-2", variants = x, tools = "mr")

test_that("create summaryset", {
    expect_true(nrow(sumset1@ss) == length(x))
})

test_that("set metadata", {
    sumset1 <- setMetadata(sumset1,source = "IEUopenGWAS", traits = "ieu-a-2")
    expect_true(is.list(sumset1@metadata))
})


test_that("get metadata", {
    expect_equal(getMetadata(sumset1), sumset1@metadata)
})

test_that("set mr_label", {
    sumset1 <- setMRlabel(sumset1, mr_label = "exposure")
    expect_equal
})


test_that("getRSID", {
    expect_equal(getRSID(sumset1), x)
})


test_that("dim works", {
    expect_equal(dim(sumset1@ss), dimData(sumset1))
})

  
