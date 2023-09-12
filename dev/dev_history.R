#  this file contains some code used during development

library(devtools)
load_all()



# Checks
devtools::check()
devtools::document() # update namespace, etc
pkgdown::build_site()
usethis::use_tidy_description()

roxygen2::update_collate(".")

# unit testing
devtools::test()

# create r files
usethis::use_r("liftover.R")

# create a vignette
usethis::use_vignette("SummarySet_DataSet")
usethis::use_vignette("liftover")

# coverage
library(covr)
report()
usethis::use_github_action("test-coverage")

# add files to .Rbuildignore
usethis::use_build_ignore("dev/", escape = TRUE)
usethis::use_build_ignore("docker/", escape = TRUE)
usethis::use_build_ignore(".dockerignore/", escape = TRUE)
usethis::use_build_ignore(".tests/", escape = TRUE)
usethis::use_build_ignore(".inst/tests/", escape = TRUE)
usethis::use_build_ignore(".devcontainer/", escape = TRUE)
usethis::use_build_ignore("data/", escape = TRUE)

# add files to .gitignore
usethis::use_git_ignore("inst/old/test.R")

#  add NEWS.md
usethis::use_news_md()

# debug vignettes
build_vignettes()
devtools::build_rmd("vignettes/liftover.Rmd") 
