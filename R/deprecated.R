##############################################################################
#    Copyright (c) 2012-2017 Russell V. Lenth                                #
#                                                                            #
#    This file is part of the emmeans package for R (*emmeans*)              #
#                                                                            #
#    *emmeans* is free software: you can redistribute it and/or modify       #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 2 of the License, or       #
#    (at your option) any later version.                                     #
#                                                                            #
#    *emmeans* is distributed in the hope that it will be useful,            #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       #
#    along with R and *emmeans*.  If not, see                                #
#    <https://www.r-project.org/Licenses/> and/or                            #
#    <http://www.gnu.org/licenses/>.                                         #
##############################################################################

### Deprecated code - Structural changes Sep 2017
# Documentation stuff is at the end


## 0. We are NOT deprecating lsmeans(), lstrends(), lsm(), etc.
##    We are preserving these into the future, but they are remapped
##    to emmeans(), etc. See file synonyms.R


## 1. We're changing ref.grid and lsmobj classes to emmGrid
# setClass("ref.grid", contains = "emmGrid")
# setClass("lsmobj", contains = "ref.grid")

# For reasons I don't quite understand, I need a new plot.ref.grid method
plot.ref.grid = function(x, ...) {
    plot(x = as(x, "emmGrid"), ...)
}

## 2. Change "dot" names to "underscore" names...
#' @rdname emmeans-deprecated
#' @param ... Arguments passed to other methods
#' @export
ref.grid = function(...) {
    .Deprecated(new = "ref_grid",  old = "ref.grid", package = "emmeans")
    ref_grid(...)
}

### These are in wrappers.R ###
# #' @rdname emmeans-deprecated
# #' @export
# lsm.options = function(...) {
#     .Deprecated("emm_options")
#     emm_options(...)
# }
# #' @rdname emmeans-deprecated
# #' @export
# get.lsm.option = function(...) {
#     .Deprecated("get_emm_option")
#     get_emm_option(...)
# }



## 3. Temporary support for old recover.data and lsm.basis methods

# Keep old generics for now, and add .chk_ routines (which we call in ref_grid and emtrends)
#' @rdname emmeans-deprecated
#' @param object An \code{emmGrid} object
#' @export
recover.data = function(object, ...)
    UseMethod("recover.data")

#' @method recover.data call
#' @export
recover.data.call = function(...)
    recover_data.call(...)


.chk_recover_data = function(object, ...) {  # default method returns char string if it fails
    has.old.meth = paste0("recover.data.", 
                          class(object)[1]) %in% methods("recover.data")
    if (has.old.meth) {
 ## keep it silent      .Deprecated(new = "recover_data", old = "recover.data", package = "emmeans")
        recover.data(object, ...)
    }
    else
        recover_data(object, ...)
}

#' @rdname emmeans-deprecated
#' @export
lsm.basis = function(object, ...)
    UseMethod("lsm.basis")

.chk_emm_basis = function(object, ...) {
    has.old.meth = paste0("lsm.basis.", 
                          class(object)[1]) %in% methods("lsm.basis")
    if (has.old.meth) {
## keep it silent       .Deprecated(new = "emm_basis", old = "lsm.basis", package = "emmeans")
        lsm.basis(object, ...)
    }
    else
        emm_basis(object, ...)
}

#' Deprecated functions and arguments from \pkg{lsmeans}
#' 
#' Over time, some functions and internal structures have been revised or
#' expanded, and as a consequence this necessitates renaming or
#' reconceptualization of the functions it exports. This is a quick reference on
#' how to get what you want if what worked in the past no longer does.
#' 
#' @section List of deprecated functions and arguments:
#' \newcommand{\FCN}{\item{\code{#1()}}}
#' \newcommand{\ARG}{\item{\code{#1}}}
#' \newcommand{\cls}{\dQuote{\code{#1}} class}
#' \newcommand{\CLS}{\item{\cls{#1}}}
#' \describe{
#' \FCN{as.stanfit}{We suggest using \code{as.mcmc()} instead, and plotting
#' the results using functions in \pkg{bayesplot}.}
#' .
#' \CLS{lsmobj}{Both this  and the \cls{ref.grid} have been replaced by the
#' \cls{emmGrid}.}
#' 
#' \FCN{ref.grid}{This has been replaced by \code{ref_grid()}, in hopes of
#' reducing the chance that \code{ref.grid} will be mistaken as an S3 method for
#' class \code{grid}.}
#' 
#' \CLS{ref.grid}{Both this  and the \cls{lsmobj} have been replaced by the
#' \cls{emmGrid}.}
#' 
#' \ARG{trend}{The \code{trend} argument in \code{lsmeans} (now \code{emmeans})
#' is now deprecated. Use \code{emtrends()} instead.}
#' 
#' }% list of deprecated stuff
#' @name emmeans-deprecated
NULL

## Here is a utility that we won't export, but can help clean out lsmeans
## stuff from one's workspace, and unload unnecessary junk
convert_workspace = function(envir = .GlobalEnv) {
    if (exists(".Last.ref.grid", envir = envir)) {
        cat("Deleted .Last.ref.grid\n")
        remove(".Last.ref.grid", envir = envir)
    }
    for (nm in names(envir)) {
        obj <- get(nm)
        if (is(obj, "ref.grid")) {
            cat(paste("Converted", nm, "to class 'emmGrid'\n"))
            assign(nm, as.emmGrid(obj), envir = envir)
        }
    }
    if ("package:lsmeans" %in% search())
        detach("package:lsmeans")
    if ("lsmeans" %in% loadedNamespaces())
        unloadNamespace("lsmeans")
    message("The environment has been converted and lsmeans's namespace is unloaded.\n",
            "Now you probably should save it.")
}


## Here is a non-exported utility to convert .R and .Rmd files
## It's entirely menu-driven.
convert_scripts = function() {
    infiles = utils::choose.files(
        caption = "Select R script(s) or markdown file(s) to be converted",
        multi = TRUE)
    lsm.to.emmGrid = utils::menu(c("yes", "no"), graphics = TRUE, 
                      "lsmxxx() -> emmxxx()?") == 1
    pmm.to.emmGrid = utils::menu(c("yes", "no"), graphics = TRUE, 
                      "pmmxxx() -> emmxxx()?") == 1
    
    for (infile in infiles) {
        buffer = scan(infile, what = character(0), sep = "\n", 
                      blank.lines.skip = FALSE)
        
        buffer = gsub("library *\\(\"*'*lsmeans\"*'*\\)", "library(\"emmeans\")", buffer)
        buffer = gsub("require *\\(\"*'*lsmeans\"*'*\\)", "require(\"emmeans\")", buffer)
        buffer = gsub("lsmeans::", "emmeans::", buffer)
        buffer = gsub("ref\\.grid *\\(", "ref_grid(", buffer)
        opt.idx = grep("lsm\\.option", buffer)
        if (length(opt.idx) > 0) {
            buffer[opt.idx] = gsub("ref\\.grid", "ref_grid", buffer[opt.idx])
            buffer[opt.idx] = gsub("lsmeans", "emmeans", buffer[opt.idx])
            buffer[opt.idx] = gsub("lsm\\.options *\\(", "emm_options(", buffer[opt.idx])
            buffer[opt.idx] = gsub("get\\.lsm\\.option *\\(", "get_emm_option(", buffer[opt.idx])
        }
        buffer = gsub("\\.lsmc", ".emmc", buffer)
        
        if (lsm.to.emmGrid) {
            buffer = gsub("lsmeans *\\(", "emmeans(", buffer)
            buffer = gsub("lsmip *\\(", "emmip(", buffer)
            buffer = gsub("lstrends *\\(", "emtrends(", buffer)
            buffer = gsub("lsm *\\(", "emmGrid(", buffer)
            buffer = gsub("lsmobj *\\(", "emmobj(", buffer)
        }
        if (pmm.to.emmGrid) {
            buffer = gsub("pmmeans *\\(", "emmeans(", buffer)
            buffer = gsub("pmmip *\\(", "emmip(", buffer)
            buffer = gsub("pmtrends *\\(", "emtrends(", buffer)
            buffer = gsub("pmm *\\(", "emmGrid(", buffer)
            buffer = gsub("pmmobj *\\(", "emmobj(", buffer)
        }
        outfile = file.path(dirname(infile), sub("\\.", "-emmGrid.", basename(infile)))
        write(buffer, outfile)
        cat(paste(infile, "\n\twas converted to\n", outfile, "\n"))
    }
}