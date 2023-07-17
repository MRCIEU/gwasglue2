# Need extensive testing of different harmonisation scenarios

scenarios <- dplyr::as_tibble(read.csv(system.file("tests", "harmonisation_scenarios.csv", package="gwasglue2")))

# Get all scenarios, remote flip2 for now
scen_names <- unique(scenarios$scenario)
scen_names <- scen_names[!scen_names %in% c("flip2")]

# test each scenario 1 by 1
# - Create dataset of input
# - Extract harmonised data
# - Compare against truth

lapply(scen_names, function(s) {
  test_that(s, {
  # input
  x <- subset(scenarios, scenario == s)
  input <- lapply(unique(x$id), function(i) {
      scenarios %>% dplyr::filter(scenario == s, id == i, version == "input")
    })

  # for "palindrome_flip" action=2, other scen action=1
  action <- ifelse(s %in% c("palindrome_flip"), 2, 1)
  
  # create dataset
  dataset <- lapply(seq_along(input), function(i) {
    dt <- create_summaryset(input[[i]]) }) %>%
    create_dataset(., harmonise = TRUE, tolerance = 0.08, action = action) %>%
        suppressMessages()

  result <- lapply(1:getLength(dataset), function(i) {
    d <- getSummarySet(dataset, i) %>%
        getSummaryData(.) %>%
        dplyr::select(chr, ea, nea, eaf, beta) %>%
        suppressMessages()
    })
  
  # Truth
   truth <- lapply(unique(x$id), function(i) {
    d <- scenarios %>% 
      dplyr::filter(scenario == s, id == i, version == "truth")  %>%
      dplyr::select(chr, ea, nea, eaf, beta)})
  
  expect_equal(result, truth)
  })
})

# TODO flip2