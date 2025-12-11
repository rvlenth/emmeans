
# to load all functions into an R session
devtools::load_all()

# to rebuild documentation:
roxygen2::roxygenise()

# to run all unit tests
testthat::test_dir("tests/testthat")

# to build the entire package
devtools::build()

# required check
devtools::check()