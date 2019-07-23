##############################################################################
#    Copyright (c) 2012-2019 Russell V. Lenth                                #
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



# Just define the function for now. When we get to R version 3.6 or so
# maybe we can we require R >= 3.4 (first that has hasName())
# and add utils::hasName to imports (in emmeans-package.R)

hasName = function(x, name)
   match(name, names(x), nomatch = 0L) > 0L



### NOTE: Revised just after version 1.3.1 release to move CSS file to inst/css
###       because devtools and relatives will delete inst/doc without notice!

# NOTE: Excluded from documentation
# Custom Vignette format
# 
# This is used to format HTML vignettes the way its developer wants them.
#
# @param ... Arguments passed to \code{rmarkdown::html_document}
#
# @return R Markdown format used by \code{rmarkdown::render}
#' @export
.emm_vignette = function(css = system.file("css", "clean-simple.css", package = "emmeans"), 
                         highlight = NULL,  ...) {
    rmarkdown::html_document(theme = NULL, highlight = highlight,
                             fig_width = 3, fig_height = 3, 
                             css = css, pandoc_args = "", ...)
###    css = css, pandoc_args = "--strip-comments", ...)
}


### Dynamic registration of S3 methods
# Code borrowed from hms pkg. I omitted some type checks etc. because
# this is only for internal use and I solemnly promise to behave myself.
register_s3_method = function(pkg, generic, class, envir = parent.frame()) {
    fun = get(paste0(generic, ".", class), envir = envir)
    if (isNamespaceLoaded(pkg)) {
        registerS3method(generic, class, fun, envir = asNamespace(pkg))
    }
    # Register hook in case package is later unloaded & reloaded
    setHook(
        packageEvent(pkg, "onLoad"),
        function(...) {
            registerS3method(generic, class, fun, envir = asNamespace(pkg))
        }
    )
}

.onLoad = function(libname, pkgname) {
    if (.requireNS("coda", fail = .nothing)) {
        register_s3_method("coda", "as.mcmc", "emmGrid")
        register_s3_method("coda", "as.mcmc.list", "emmGrid")
    }
    if (.requireNS("multcomp", fail = .nothing)) {
        register_s3_method("multcomp", "glht", "emmlf")
        register_s3_method("multcomp", "glht", "emmGrid")
        register_s3_method("multcomp", "cld", "emmGrid")
        register_s3_method("multcomp", "modelparm", "emmwrap")
    }
}


#' @rdname extending-emmeans
#' @section Registering S3 methods for a model class:
#' The \code{.emm_register} function is provided as a convenience to conditionally 
#' register your
#' S3 methods for a model class, \code{recover_data.foo} and \code{emm_basis.foo},
#' where \code{foo} is the class name. Your package should implement an
#' \code{.onLoad} function and call \code{.emm_register} if \pkg{emmeans} is
#' installed. See the example.
#'
#' @param classes Character names of one or more classes to be registered.
#'   The package must contain the functions \code{recover_data.foo} and
#'   \code{emm_basis.foo} for each class \code{foo} listed in \code{classes}.
#' @param pkgname Character name of package providing the methods (usually
#'    should be the second argument of \code{.onLoad})
#'
#' @export
#'
#' @examples
#' \dontrun{
#' #--- If your package provides recover_data and emm_grid methods for class 'mymod',
#' #--- put something like this in your package code -- say in zzz.R:
#'   .onLoad = function(libname, pkgname) {
#'     if (requireNamespace("emmeans", quietly = TRUE))
#'       emmeans::.emm_register("mymod", pkgname)
#'   }
#' }
.emm_register = function(classes, pkgname) {
    envir = asNamespace(pkgname)
    for (class in classes) {
        register_s3_method("emmeans", "recover_data", class, envir)
        register_s3_method("emmeans", "emm_basis", class, envir)
    }
}

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
        outfile = file.path(dirname(infile), sub("\\.", "-emm.", basename(infile)))
        write(buffer, outfile)
        cat(paste(infile, "\n\twas converted to\n", outfile, "\n"))
    }
}
