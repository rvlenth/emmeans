
# build tar.gz file, navigate to its directory and run the following script

reverse_check = function(reverse = list(which = "all"), Ncpus = 3) {
    
  tar = dir(pattern = "*.tar.gz")
  if(length(tar) == 0)
      stop("No .tar.gz files found")
  ans = menu(tar, graphics = FALSE, title = "Choose tarball, or 0 to cancel")
  if (ans == 0) return(NULL)
  tar = tar[ans]
  
  # path to temp directory
  nm = gsub(".tar.gz", "", tar)
  tmp = strsplit(tempdir(), "\\\\")[[1]]
  tmp = paste(tmp[-length(tmp)], collapse = "\\")
  dir = paste0(tmp, "\\checkdir.for.", nm)
  dir.create(dir)
  
  file.copy(tar, dir)
  
  message("Checking ", tar, "\nDirectory for results: ", dir,
          "\n\nThis will take some time...")
  chk = tools::check_packages_in_dir(dir, reverse = reverse, 
                                     Ncpus = Ncpus, clean = TRUE)
  sink(paste(dir, "00summary.txt", sep = "\\"))
  print(summary(chk))
  sink(NULL)
  message("\nDetails are in ", dir,"\n(See also 00summary.txt there)")
  invisible(chk)
}


