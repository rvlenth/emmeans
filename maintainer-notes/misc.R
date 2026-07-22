# for removing trailing white space

ws_function <- function(path) {
    # for a single file
  
    # 1. Read the content of the file
    lines <- readLines(path)

    # 2. Use gsub to remove trailing whitespace
    # The regex ' +$' matches one or more spaces at the end of a line
    clean_lines <- gsub(" +$", "", lines)

    # 3. Write the cleaned lines back to the file
    writeLines(clean_lines, path)

    message("Trailing whitespace removed from: ", path)
}

dir_function <- function(dirc) {
  # to apply to a directory
    files <- as.list(list.files(dirc, pattern = "\\.(R|Rmd)$", ignore.case = TRUE))
    lapply(files, function(f) {
      long_path = here::here(dirc, f)
      ws_function(long_path)
        }
      )
    }

### testing
#ws_function("Vignettes/xtending.Rmd")
#dir_function(getwd()) # should do nothing
#dir_function("vignettes") # gets all .Rmd files
#dir_function("R") # gets all .R files

