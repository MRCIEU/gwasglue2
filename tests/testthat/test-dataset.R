
x <- ieugwasr::tophits("ieu-a-2")$rsid
d1 <- ieugwasr::associations(variants = x, id = "ieu-a-2")
d2 <- ieugwasr::associations(variants = x, id = "ieu-a-7")
sumset1 <- constructSummarySet(d1, tools = "mr", source = "IEUopenGWAS", id = "ieu-a-2")
sumset2 <- constructSummarySet(d2, tools ="mr", source = "IEUopenGWAS", id = "ieu-a-7")




test_that("compare against 2samplemr", {
 
  dataset <- DataSet(sumset1,sumset2) %>%
    overlapVariants(.) %>%
    harmoniseData(.,tolerance = 0.08, action = 1,strand = "forward")
  dataset
  expect_equal(length(dataset@summary_sets), 2)
})


