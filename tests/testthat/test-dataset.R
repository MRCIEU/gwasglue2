skip("pipeline changed")
x <- ieugwasr::tophits("ieu-a-2")$rsid
sumset1 <- SummarySet(traits = "ieu-a-2", variants = x, tools = "mr")
sumset2 <- SummarySet(traits="ieu-a-7", variants=x,tools ="mr")

test_that("compare against 2samplemr", {
 
  dataset <- DataSet(sumset1,sumset2) %>%
    overlapSNP(.) %>%
    harmoniseData(.,tolerance = 0.08, action = 1)
  dataset
  expect_equal(length(dataset@summary_sets), 2)
})


