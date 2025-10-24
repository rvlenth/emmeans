---
title: "Notes on preparation and submission of **emmeans**"
author: "Russ Lenth"
date: "2025-10-09"
output: emmeans::.emm_vignette
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Getting it ready for CRAN submission

### Vignette index
The vignettes (but not this one) contain markup within HTML comments that compile a
simple index. If any of these are added or changed, you need to manually update
the index, and that is done by running
```
vigindex::vigindex()
```
This creates another vignette named `vigindex.Rmd`.
The package is available via [https://github.com/rvlenth/vigindex](https://github.com/rvlenth/vigindex).

Here are some typical tags:
```
<!-- @index Transformations; 
            Transformations!Bias correction;
            Examples!`pigs`   -->
```
The comment close *must* be on the same line as the last entry. The second and
third examples create sub-entries. In those cases it's a good idea to ensure
that there is not already a main entry that is similar, but not identical, to
the first part.

### Other manual checking
Make sure the version numbers in `DESCRIPTION` and `NEWS.md` match.

*And* look at the [emmeans package page](https://cran.r-project.org/web/packages/emmeans/index.html) on CRAN,
and go the "CRAN checks: emmeans results." Make sure there are no issues there.

### Build and check the tarball
In RStudio, this is done via the `build` tab, then under "More `v`", 
choose "Build source package. This will create a file named `emmeans-xxx.tar.gz`
in the container of the project directory (just above the `emmeans` directory).

Then check to make sure there are no errors. I always do this in a DOS window (or shell):
```
R CMD check --as-cran emmeans-xxx.tar.gz
```
Obviously, you need to correct any errors that are found.

*Note:* This creates a directory `emmeans.Rcheck/` at the same level as `emmeans/`.
Sometimes I get errors if this directory already exists, because my system sometimes thinks
I have an open file in there somewhere. So I manually delete that directory before checking.

### Doing a thorough check of related packages
Installs all reverse dependencies and checks functionality of updated package.

This is a pain in the butt and optional, but it can prevent surprises that may
come up when CRAN checks the package. I have a script that I developed recently.
It takes a very long time and I usually start it before bedtime, then go through
the results the next day. Here is the script:
```
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
```
After running it, look at the summary, and for each check log file associated 
with an error, look at what happened. Sometimes it's just an unavailable package
or other nuisance. It's a painful process.

If there is a problem, contact the package developer. These issues need to be resolved (CRAN will ask). 

### build a TarBall


### Check against R-Devel
Go to the [Win-Builder site](https://win-builder.r-project.org/upload.aspx)
where they give you three options. Pick the second one (R-devel), which is what
they check against when the package is submitted. (You don't need to check with
R-release, because presumably you already did that earlier.) Click on "Choose file",
and navigate to the `tar.gz` file; then click on on "Submit file." It'll take
an hour or so, and if you get a return message with "Status: OK," you're ready to
submit it to CRAN.

WinBuilder will respond in ~1 hour, and will provide a binary. It's running `R CMD Check`(?). 
They will provide a check log. The desired outcome is "status: OK". If it's something else,
read the check log attached. 

### Submit the package to CRAN
Go to the [CRAN submissions site](https://cran.r-project.org/submit.html), enter
your name and email, choose the `tar.gz` file, and click on "Upload package."
You're not done yet! As maintainer, you'll get an email confirming that you
really are the one who submitted it. Click on the link they give. You'll need to
verify that you've read all the policies and looked at the CRAN and Winbuilder
checks, etc. Then submit the package and wait to hear from them. You'll get a
routine reply once it's gone through preliminary checks. Then within a few
hours, hopefully you'll get an email that it's on its way to CRAN. Or you'll get
the dreaded email that it failed a check against another package; then you have
to respond to that. Let's just assume that won't happen...