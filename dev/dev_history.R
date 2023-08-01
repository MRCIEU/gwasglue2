#  this file contains some code used during development

library(devtools)
load_all()



# Checks
devtools::check()
devtools::document() # update namespace, etc
pkgdown::build_site()



 roxygen2::update_collate(".")

# unit testing
devtools::test()


# create a vignette
usethis::use_vignette("SummarySet_DataSet")


library(covr)
report()
usethis::use_github_action("test-coverage")

# ignore R build
usethis::use_build_ignore("dev/", escape = TRUE)
usethis::use_build_ignore("docker/", escape = TRUE)
usethis::use_build_ignore(".dockerignore/", escape = TRUE)
usethis::use_build_ignore(".devcontainer/", escape = TRUE)
usethis::use_build_ignore("data/", escape = TRUE)
usethis::use_build_ignore("ieugwasr_oauth/", escape = TRUE)
usethis::use_build_ignore(".gitignore", escape = TRUE)
usethis::use_build_ignore(".vscode/", escape = TRUE)
usethis::use_build_ignore("tests/", escape = TRUE)

# ignore git
usethis::use_git_ignore("data/")
usethis::use_git_ignore("ieugwasr_oauth/")
usethis::use_git_ignore(".vscode/")
usethis::use_git_ignore(".devcontainer/")


usethis::use_news_md()

# debug vignettes
build_vignettes()
devtools::build_rmd("vignettes/Strategy.Rmd") 
