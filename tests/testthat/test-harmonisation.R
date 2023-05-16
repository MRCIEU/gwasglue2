# Need extensive testing of different harmonisation scenarios
library(dplyr)




test_that("harmonisation scenario: easy", {
test <- "easy"
ids <- c(1,2)
action  <- 1
scenarios <- read.csv("harmonisation_scenarios.csv")
scenarios <- as_tibble(scenarios)
tolerance <- 0.08

df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset(data=df, harmonise = TRUE, tolerance = tolerance, action = action)


result <- lapply(seq_along(dataset@summary_sets), function(i) {
  d <- dataset@summary_sets[[i]]@ss %>% select(chr, ea, nea, eaf, beta)})

# Truth
truth<- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "truth")  %>% 
    select(chr, ea, nea, eaf, beta)})

lapply(seq_along(ids), function(i) {
  expect_equal(truth[[i]], result[[i]])})

})

test_that("harmonisation scenario: alphabetical", {
test <- "alphabetical"
ids <- c(1,2)
action  <- 1
scenarios <- read.csv("harmonisation_scenarios.csv")
scenarios <- as_tibble(scenarios)
tolerance <- 0.08

df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset(data=df, harmonise = TRUE, tolerance = tolerance, action = action)


result <- lapply(seq_along(dataset@summary_sets), function(i) {
  d <- dataset@summary_sets[[i]]@ss %>% select(chr, ea, nea, eaf, beta)})

# Truth
truth<- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "truth")  %>% 
    select(chr, ea, nea, eaf, beta)})

lapply(seq_along(ids), function(i) {
  expect_equal(truth[[i]], result[[i]])})

})



test_that("harmonisation scenario: flip1", {
test <- "flip1"
ids <- c(1,2)
action  <- 1
scenarios <- read.csv("harmonisation_scenarios.csv")
scenarios <- as_tibble(scenarios)
tolerance <- 0.08

df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset(data=df, harmonise = TRUE, tolerance = tolerance, action = action)


result <- lapply(seq_along(dataset@summary_sets), function(i) {
  d <- dataset@summary_sets[[i]]@ss %>% select(chr, ea, nea, eaf, beta)})

# Truth
truth<- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "truth")  %>% 
    select(chr, ea, nea, eaf, beta)})

lapply(seq_along(ids), function(i) {
  expect_equal(truth[[i]], result[[i]])})

})


# TODO


# test_that("harmonisation scenario: flip2", {
# test <- "flip2"

# ids <- c(1,2)
# action <- 2

# df <- lapply(seq_along(ids), function(i) {
#   d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

# dataset <- create_dataset(data=df, harmonise = TRUE, tolerance = tolerance, action = action)


# result <- lapply(seq_along(dataset@summary_sets), function(i) {
#   d <- dataset@summary_sets[[i]]@ss %>% select(chr, ea, nea, eaf, beta)})

# # Truth
# truth<- lapply(seq_along(ids), function(i) {
#   d <- scenarios %>% filter(scenario == test, id == ids[i],version == "truth")  %>% 
#     select(chr, ea, nea, eaf, beta)})


# lapply(seq_along(ids), function(i) {
#   expect_equal(truth[i], result[i])})

# })

test_that("harmonisation scenario:  palindrome_flip", {
test <-  "palindrome_flip"
scenarios <- read.csv("harmonisation_scenarios.csv")
scenarios <- as_tibble(scenarios)
tolerance <- 0.08
ids <- c(1,2)
action <- 2
df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset(data=df, harmonise = TRUE, tolerance = tolerance, action = action)


result <- lapply(seq_along(dataset@summary_sets), function(i) {
  d <- dataset@summary_sets[[i]]@ss %>% select(chr, ea, nea, eaf, beta)})

# Truth
truth<- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "truth")  %>% 
    select(chr, ea, nea, eaf, beta)})

lapply(seq_along(ids), function(i) {
  expect_equal(truth[[i]], result[[i]])})
})

test_that("harmonisation scenario:  palindrome_noflip", {
test <-  "palindrome_noflip"
scenarios <- read.csv("harmonisation_scenarios.csv")
scenarios <- as_tibble(scenarios)
tolerance <- 0.08
ids <- c(1,2)
action <- 1
df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset(data=df, harmonise = TRUE, tolerance = tolerance, action = action)


result <- lapply(seq_along(dataset@summary_sets), function(i) {
  d <- dataset@summary_sets[[i]]@ss %>% select(chr, ea, nea, eaf, beta)})

# Truth
truth<- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "truth")  %>% 
    select(chr, ea, nea, eaf, beta)})

lapply(seq_along(ids), function(i) {
  expect_equal(truth[[i]], result[[i]])})
})

test_that("harmonisation scenario:  easy_indels", {
scenarios <- read.csv("harmonisation_scenarios.csv")
scenarios <- as_tibble(scenarios)
tolerance <- 0.08
test <-  "easy_indels"
ids <- c(1,2)
action  <- 1

df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset(data=df, harmonise = TRUE, tolerance = tolerance, action = action)


result <- lapply(seq_along(dataset@summary_sets), function(i) {
  d <- dataset@summary_sets[[i]]@ss %>% select(chr, ea, nea, eaf, beta)})

# Truth
truth<- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "truth")  %>% 
    select(chr, ea, nea, eaf, beta)})

lapply(seq_along(ids), function(i) {
  expect_equal(truth[[i]], result[[i]])})
})

test_that("harmonisation scenario:  indels_flip1", {
scenarios <- read.csv("harmonisation_scenarios.csv")
scenarios <- as_tibble(scenarios)
tolerance <- 0.08
test <-  "indels_flip1"
ids <- c(1,2)
action  <- 1

df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset(data=df, harmonise = TRUE, tolerance = tolerance, action = action)


result <- lapply(seq_along(dataset@summary_sets), function(i) {
  d <- dataset@summary_sets[[i]]@ss %>% select(chr, ea, nea, eaf, beta)})

# Truth
truth<- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "truth")  %>% 
    select(chr, ea, nea, eaf, beta)})

lapply(seq_along(ids), function(i) {
  expect_equal(truth[[i]], result[[i]])})
})

test_that("harmonisation scenario:  indels_drop1", {
scenarios <- read.csv("harmonisation_scenarios.csv")
scenarios <- as_tibble(scenarios)
tolerance <- 0.08
test <-  "indels_drop1"
ids <- c(1,2)
action  <- 1

df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset(data=df, harmonise = TRUE, tolerance = tolerance, action = action)


result <- lapply(seq_along(dataset@summary_sets), function(i) {
  d <- dataset@summary_sets[[i]]@ss %>% select(chr, ea, nea, eaf, beta)})

# Truth
truth<- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "truth")  %>% 
    select(chr, ea, nea, eaf, beta)})

lapply(seq_along(ids), function(i) {
  expect_equal(truth[[i]], result[[i]])})
})

test_that("harmonisation scenario:  multiallele1", {
scenarios <- read.csv("harmonisation_scenarios.csv")
scenarios <- as_tibble(scenarios)
tolerance <- 0.08
test <-  "multiallele1"
ids <- c(1,2)
action  <- 1

df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset(data=df, harmonise = TRUE, tolerance = tolerance, action = action)


result <- lapply(seq_along(dataset@summary_sets), function(i) {
  d <- dataset@summary_sets[[i]]@ss %>% select(chr, ea, nea, eaf, beta)})

# Truth
truth<- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "truth")  %>% 
    select(chr, ea, nea, eaf, beta)})

lapply(seq_along(ids), function(i) {
  expect_equal(truth[[i]], result[[i]])})
})

test_that("harmonisation scenario:  multiallele2", {
scenarios <- read.csv("harmonisation_scenarios.csv")
scenarios <- as_tibble(scenarios)
tolerance <- 0.08
test <-  "multiallele2"
ids <- c(1,2)
action  <- 1

df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset(data=df, harmonise = TRUE, tolerance = tolerance, action = action)


result <- lapply(seq_along(dataset@summary_sets), function(i) {
  d <- dataset@summary_sets[[i]]@ss %>% select(chr, ea, nea, eaf, beta)})

# Truth
truth<- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "truth")  %>% 
    select(chr, ea, nea, eaf, beta)})

lapply(seq_along(ids), function(i) {
  expect_equal(truth[[i]], result[[i]])})
})


test_that("harmonisation scenario:  multiallele3", {
scenarios <- read.csv("harmonisation_scenarios.csv")
scenarios <- as_tibble(scenarios)
tolerance <- 0.08
test <-  "multiallele3"
ids <- c(1,2)
action  <- 1

df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset(data=df, harmonise = TRUE, tolerance = tolerance, action = action)


result <- lapply(seq_along(dataset@summary_sets), function(i) {
  d <- dataset@summary_sets[[i]]@ss %>% select(chr, ea, nea, eaf, beta)})

# Truth
truth<- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "truth")  %>% 
    select(chr, ea, nea, eaf, beta)})

lapply(seq_along(ids), function(i) {
  expect_equal(truth[[i]], result[[i]])})
})



test_that("harmonisation scenario:  multiallele4", {
scenarios <- read.csv("harmonisation_scenarios.csv")
scenarios <- as_tibble(scenarios)
tolerance <- 0.08
test <-  "multiallele4"
ids <- c(1,2)
action  <- 1

df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset(data=df, harmonise = TRUE, tolerance = tolerance, action = action)


result <- lapply(seq_along(dataset@summary_sets), function(i) {
  d <- dataset@summary_sets[[i]]@ss %>% select(chr, ea, nea, eaf, beta)})

# Truth
truth<- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "truth")  %>% 
    select(chr, ea, nea, eaf, beta)})

lapply(seq_along(ids), function(i) {
  expect_equal(truth[[i]], result[[i]])})
})

