
# to load all functions into an R session
devtools::load_all()

# to rebuild documentation:
roxygen2::roxygenise()

# to run all unit tests
testthat::test_dir("tests/testthat")

# to build the entire package
devtools::build(cran = TRUE, manual = TRUE, vignettes = TRUE)

# required check
devtools::check()