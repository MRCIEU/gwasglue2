# Need extensive testing of different harmonisation scenarios
x <- ieugwasr::tophits("ieu-a-2")$rsid
sumset1 <- SummarySet(traits = "ieu-a-2", variants = x, tools = "mr")
sumset2 <- SummarySet(traits="ieu-a-7", variants=x,tools ="mr")

test_that("simple 2 trait harmonisation matches TwoSampleMR", {

  dataset <- DataSet(sumset1,sumset2) %>%
    overlapSNP(.) %>%
    harmoniseData(.,tolerance = 0.08, action = 1)
  # Do the same in TwoSampleMR
  dat <- TwoSampleMR::make_dat(2,7)

  expect_equal(
    cor(dataset@summary_sets[[1]]@ss$beta, dataset@summary_sets[[2]]@ss$beta),
    cor(dat$beta.exposure, dat$beta.outcome)
  )
})
