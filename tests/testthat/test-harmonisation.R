# Need extensive testing of different harmonisation scenarios
library(dplyr)

scenarios <- as_tibble(read.csv(system.file("tests", "harmonisation_scenarios.csv", package="gwasglue2")))

test_that("harmonisation scenario: easy", {
test <- "easy"

ids <- 1:2

df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset_from_tibble(data=df, harmonise = TRUE, tolerance = 0.08, action = 1)


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
# scenarios <- as_tibble(read.csv("harmonisation_scenarios.csv"))
ids <- 1:2


df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset_from_tibble(data=df, harmonise = TRUE, tolerance = 0.08, action = 1)


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
# scenarios <- as_tibble(read.csv("harmonisation_scenarios.csv"))
ids <- 1:2


df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset_from_tibble(data=df, harmonise = TRUE, tolerance = 0.08, action = 1)


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

# scenarios <- as_tibble(read.csv("harmonisation_scenarios.csv"))
# ids <- 1:2


# df <- lapply(seq_along(ids), function(i) {
#   d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

# dataset <- create_dataset_from_tibble(data=df, harmonise = TRUE, tolerance = 0.08, action = 2)


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
# scenarios <- as_tibble(read.csv("harmonisation_scenarios.csv"))
ids <- 1:2


df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset_from_tibble(data=df, harmonise = TRUE, tolerance = 0.08, action = 2)

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
# scenarios <- as_tibble(read.csv("harmonisation_scenarios.csv"))
ids <- 1:2


df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset_from_tibble(data=df, harmonise = TRUE, tolerance = 0.08, action = 1)


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
test <-  "easy_indels"
# scenarios <- as_tibble(read.csv("harmonisation_scenarios.csv"))
ids <- 1:2


df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset_from_tibble(data=df, harmonise = TRUE, tolerance = 0.08, action = 1)


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
test <-  "indels_flip1"
# scenarios <- as_tibble(read.csv("harmonisation_scenarios.csv"))
ids <- 1:2


df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset_from_tibble(data=df, harmonise = TRUE, tolerance = 0.08, action = 1)


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
test <-  "indels_drop1"
# scenarios <- as_tibble(read.csv("harmonisation_scenarios.csv"))
ids <- 1:2


df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset_from_tibble(data=df, harmonise = TRUE, tolerance = 0.08, action = 1)


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
test <-  "multiallele1"
# scenarios <- as_tibble(read.csv("harmonisation_scenarios.csv"))
ids <- 1:2


df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset_from_tibble(data=df, harmonise = TRUE, tolerance = 0.08, action = 1)


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
test <-  "multiallele2"
# scenarios <- as_tibble(read.csv("harmonisation_scenarios.csv"))
ids <- 1:2


df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset_from_tibble(data=df, harmonise = TRUE, tolerance = 0.08, action = 1)

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
test <-  "multiallele3"
# scenarios <- as_tibble(read.csv("harmonisation_scenarios.csv"))
ids <- 1:2


df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset_from_tibble(data=df, harmonise = TRUE, tolerance = 0.08, action = 1)


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
test <-  "multiallele4"
# scenarios <- as_tibble(read.csv("harmonisation_scenarios.csv"))
ids <- 1:2


df <- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "input")})

dataset <- create_dataset_from_tibble(data=df, harmonise = TRUE, tolerance = 0.08, action = 1)


result <- lapply(seq_along(dataset@summary_sets), function(i) {
  d <- dataset@summary_sets[[i]]@ss %>% select(chr, ea, nea, eaf, beta)})

# Truth
truth<- lapply(seq_along(ids), function(i) {
  d <- scenarios %>% filter(scenario == test, id == ids[i],version == "truth")  %>% 
    select(chr, ea, nea, eaf, beta)})

lapply(seq_along(ids), function(i) {
  expect_equal(truth[[i]], result[[i]])})
})

