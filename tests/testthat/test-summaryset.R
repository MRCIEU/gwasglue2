test_that("create summaryset", {
    x <- ieugwasr::tophits("ieu-a-2")$rsid
    sumset1 <- SummarySet(traits = "ieu-a-2", variants = x, tools = "mr")
    sumset1
    expect(nrow(sumset1) == length(x))
})

# create S4 SummarySet objects
# TODO Gib: Here, should we call ids instead of traits?



