# Need extensive testing of different harmonisation scenarios
  skip("Not ready yet")

scenarios <- read.csv("inst/testdata/harmonisation_scenarios.csv")




test_that("harmonisation scenario: easy", {
tolerance <- 0.08
action  <- 1
scenario <- subset(scenarios, scenarios[,1] ==  "easy")
rsid <- unique(scenario$rsid)
A1 <- scenario$a1[1]
A2 <- scenario$a2[1]
B1 <- scenario$a1[2]
B2 <- scenario$a2[2]
betaA <- scenario$beta[1]
betaB <- scenario$beta[2]
fA <- scenario$a1freq[1]
fB <- scenario$a1freq[2]

h <- harmonise(rsid, A1, A2, B1, B2, betaA, betaB, fA, fB, tolerance=tolerance, action=action)

# Truth
truth <- data.frame(rsid=rsid, A1=scenario$a1[3], A2=scenario$a2[3], B1=scenario$a1[4], B2=scenario$a2[4], betaA=scenario$beta[3], betaB=scenario$beta[4], fA=scenario$a1freq[3], fB=scenario$a1freq[4])

expect_equal(truth, h[,1:9])

})

test_that("harmonisation scenario: flip1", {
tolerance <- 0.08
action  <- 1
scenario <- subset(scenarios, scenarios[,1] ==  "flip1")
rsid <- unique(scenario$rsid)
A1 <- scenario$a1[1]
A2 <- scenario$a2[1]
B1 <- scenario$a1[2]
B2 <- scenario$a2[2]
betaA <- scenario$beta[1]
betaB <- scenario$beta[2]
fA <- scenario$a1freq[1]
fB <- scenario$a1freq[2]

h <- harmonise(rsid, A1, A2, B1, B2, betaA, betaB, fA, fB, tolerance=tolerance, action=action)

# Truth
truth <- data.frame(rsid=rsid, A1=scenario$a1[3], A2=scenario$a2[3], B1=scenario$a1[4], B2=scenario$a2[4], betaA=scenario$beta[3], betaB=scenario$beta[4], fA=scenario$a1freq[3], fB=scenario$a1freq[4])

expect_equal(truth, h[,1:9])

})

test_that("harmonisation scenario: flip2", {
tolerance <- 0.08
action  <- 1
scenario <- subset(scenarios, scenarios[,1] ==  "flip2")
rsid <- unique(scenario$rsid)
A1 <- scenario$a1[1]
A2 <- scenario$a2[1]
B1 <- scenario$a1[2]
B2 <- scenario$a2[2]
betaA <- scenario$beta[1]
betaB <- scenario$beta[2]
fA <- scenario$a1freq[1]
fB <- scenario$a1freq[2]

h <- harmonise(rsid, A1, A2, B1, B2, betaA, betaB, fA, fB, tolerance=tolerance, action=action)

# Truth
truth <- data.frame(rsid=rsid, A1=scenario$a1[3], A2=scenario$a2[3], B1=scenario$a1[4], B2=scenario$a2[4], betaA=scenario$beta[3], betaB=scenario$beta[4], fA=scenario$a1freq[3], fB=scenario$a1freq[4])

expect_equal(truth, h[,1:9])
})

test_that("harmonisation scenario: palindrome_flip", {
tolerance <- 0.08
action  <- 2
scenario <- subset(scenarios, scenarios[,1] ==  "palindrome_flip")
rsid <- unique(scenario$rsid)
A1 <- scenario$a1[1]
A2 <- scenario$a2[1]
B1 <- scenario$a1[2]
B2 <- scenario$a2[2]
betaA <- scenario$beta[1]
betaB <- scenario$beta[2]
fA <- scenario$a1freq[1]
fB <- scenario$a1freq[2]

h <- harmonise(rsid, A1, A2, B1, B2, betaA, betaB, fA, fB, tolerance=tolerance, action=action)

# Truth
truth <- data.frame(rsid=rsid, A1=scenario$a1[3], A2=scenario$a2[3], B1=scenario$a1[4], B2=scenario$a2[4], betaA=scenario$beta[3], betaB=scenario$beta[4], fA=scenario$a1freq[3], fB=scenario$a1freq[4])

expect_equal(truth, h[,1:9])

})

test_that("harmonisation scenario: palindrome_noflip", {
tolerance <- 0.08
action  <- 2
scenario <- subset(scenarios, scenarios[,1] ==  "palindrome_noflip")
rsid <- unique(scenario$rsid)
A1 <- scenario$a1[1]
A2 <- scenario$a2[1]
B1 <- scenario$a1[2]
B2 <- scenario$a2[2]
betaA <- scenario$beta[1]
betaB <- scenario$beta[2]
fA <- scenario$a1freq[1]
fB <- scenario$a1freq[2]

h <- harmonise(rsid, A1, A2, B1, B2, betaA, betaB, fA, fB, tolerance=tolerance, action=action)

# Truth
truth <- data.frame(rsid=rsid, A1=scenario$a1[3], A2=scenario$a2[3], B1=scenario$a1[4], B2=scenario$a2[4], betaA=scenario$beta[3], betaB=scenario$beta[4], fA=scenario$a1freq[3], fB=scenario$a1freq[4])

expect_equal(truth, h[,1:9])

})

test_that("harmonisation scenario: easy_indels", {
tolerance <- 0.08
action  <- 1
scenario <- subset(scenarios, scenarios[,1] ==  "easy_indels")
rsid <- unique(scenario$rsid)
A1 <- scenario$a1[1]
A2 <- scenario$a2[1]
B1 <- scenario$a1[2]
B2 <- scenario$a2[2]
betaA <- scenario$beta[1]
betaB <- scenario$beta[2]
fA <- scenario$a1freq[1]
fB <- scenario$a1freq[2]

h <- harmonise(rsid, A1, A2, B1, B2, betaA, betaB, fA, fB, tolerance=tolerance, action=action)

# Truth
truth <- data.frame(rsid=rsid, A1=scenario$a1[3], A2=scenario$a2[3], B1=scenario$a1[4], B2=scenario$a2[4], betaA=scenario$beta[3], betaB=scenario$beta[4], fA=scenario$a1freq[3], fB=scenario$a1freq[4])

expect_equal(truth, h[,1:9])

})


# compare with TwoSampleMR
x <- ieugwasr::tophits("ieu-a-2")$rsid
sumset1 <- SummarySet(traits = "ieu-a-2", variants = x, tools = "mr")
sumset2 <- SummarySet(traits="ieu-a-7", variants = x,tools ="mr")

test_that("simple 2 trait harmonisation matches TwoSampleMR", {

  dataset <- DataSet(sumset1,sumset2) %>%
    overlapVariants(.) %>%
    harmoniseData(.,tolerance = 0.08, action = 1, strand = "forward")
  # Do the same in TwoSampleMR
  dat <- TwoSampleMR::make_dat(2,7)

  expect_equal(
    cor(dataset@summary_sets[[1]]@ss$beta, dataset@summary_sets[[2]]@ss$beta),
    cor(dat$beta.exposure, dat$beta.outcome)
  )
})
