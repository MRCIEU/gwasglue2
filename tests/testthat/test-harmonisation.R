# Need extensive testing of different harmonisation scenarios
library(dplyr)


scenarios <- as_tibble(read.csv(system.file("tests", "harmonisation_scenarios.csv", package="gwasglue2")))

# Get all scenarios, remote flip2 for now
scen <- unique(scenarios$scenario)
scen <- scen[!scen %in% c("flip2")]

# test each scenario 1 by 1
# - Create dataset of input
# - Extract harmonised data
# - Compare against truth

lapply(scen, function(s)
{
  print(s)
  test_that(s, {

  # input 
  x <- subset(scenarios, scenario == s)
    df <- lapply(unique(x$id), \(i) {
      scenarios %>% filter(scenario == s, id == i, version == "input")
    })

  # for "palindrome_flip" action=2, other scen action=1
  action <- ifelse(s %in% c("palindrome_flip"), 2, 1)
  
  # create dataset
  dataset <- lapply(seq_along(data), function(i){
    dt <- create_summaryset(df[[i]]) }) %>%
    create_dataset(., harmonise = TRUE, tolerance = 0.08, action = action) %>% 
        suppressMessages()

  result <- lapply(1:getLength(dataset), function(i) {
    d <- getSummarySet(dataset, i) %>%
         getSummaryData(.) %>%
         select(chr, ea, nea, eaf, beta)})

  # Truth
  truth <- lapply(unique(x$id), \(i) {
    d <- scenarios %>% filter(scenario == s, id == i, version == "truth")  %>% 
      select(chr, ea, nea, eaf, beta)
      })
  expect_equal(truth, result)
  })
})

# TODO flip2