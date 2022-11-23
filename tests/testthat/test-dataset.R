test_that("create a dataset from summary sets", {
  sumset1 <- SummarySet(traits = "ieu-a-2", variants = x, tools = "mr")
  sumset2 <- SummarySet(traits="ieu-a-7", variants=x,tools ="mr")
  sumset3 <- SummarySet(traits="ieu-a-7", variants=x,tools ="mr")

  dataset <- DataSet(sumset1,sumset2,sumset3) %>%
    overlapSNP(.) %>%
    harmoniseData(.,tolerance = 0.08, action = 1)
  
})
