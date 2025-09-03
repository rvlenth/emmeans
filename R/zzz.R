##############################################################################
#    Copyright (c) 2012-2024 Russell V. Lenth                                #
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

### Encouragement to wean ourselvs from the titdyverse
#' Dare to be un-"tidy"!
#' 
#' Users who use \pkg{emmeans} functions as part of a pipeline -- or post-process 
#' those results in some other way -- are likely missing some important information. 
#' 
#' Your best bet is to display the actual results without any post-processing.
#' That's because \code{emmeans} and its relatives have their own \code{summary} 
#' and \code{print} methods that display annotations that may be helpful in
#' explaining what you have. If you just pipe the results into the next step,
#' those annotations are stripped away and you never see them. Statistical
#' analysis is not just a workflow; it is a discipline that involves care in
#' interpreting intermediate results, and thinking before moving on. 
#' 
#' @examples
#' neur.glm <- glm(Pain ~ Treatment + Sex + Age, family = binomial(),
#'             data = neuralgia)
#'             
#' ### The actual results with annotations (e.g. ests are on logit scale):
#' emmeans(neur.glm, "Treatment")
#' 
#' ### Post-processed results lose the annotations
#' if(requireNamespace("tibble")) {
#'     emmeans(neur.glm, "Treatment") |> tibble::as_tibble()
#' }
#' 
#' @name untidy
NULL



# NOTE: Excluded from documentation
# Custom Vignette format
# 
# This is used to format HTML vignettes the way its developer wants them.
#
# @param ... Arguments passed to \code{rmarkdown::html_document}
#
# @return R Markdown format used by \code{rmarkdown::render}
#' @rdname extending-emmeans
#' @order 51
#' @param css,package,highlight Arguments for \code{.emm_vignette}, which is
#' a clean and simple alternative to such as \code{html_document} for use
#' as the output style of a Markdown file. All the vignettes in the
#' \pkg{emmeans} package use this output style.
#' @export
.emm_vignette = function(css = system.file("css", "clean-simple.css", package = "emmeans"), 
                         highlight = NULL,  ...) {
    rmarkdown::html_document(theme = NULL, highlight = highlight,
                             fig_width = 3, fig_height = 3, 
                             css = css, ###pandoc_args = "", 
                             ...)
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
        register_s3_method("coda", "as.mcmc", "emm_list")
        register_s3_method("coda", "as.mcmc.list", "emm_list")
    }
    if (.requireNS("multcomp", fail = .nothing)) {
        register_s3_method("multcomp", "glht", "emmlf")
        register_s3_method("multcomp", "glht", "emmGrid")
        register_s3_method("multcomp", "cld", "emmGrid")
        register_s3_method("multcomp", "cld", "emm_list")
        register_s3_method("multcomp", "modelparm", "emmwrap")
    }
    if(.requireNS("xtable", fail = .nothing)) {
        register_s3_method("xtable", "xtable", "emmGrid")
        register_s3_method("xtable", "xtable", "summary_emm")
        register_s3_method("xtable", "print", "xtable_emm")
    }
}

.onAttach <- function(libname, pkgname) {
    packageStartupMessage(
        "Welcome to emmeans.\n",
        "Caution: You lose important information if you filter this package's results.\n",
        "See '? untidy'")
}



#' @rdname extending-emmeans
#' @order 29
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
#' @param qdrg Logical value. If \code{FALSE}, the \code{recover_data} and
#'     \code{emm_basis} methods are registered. If \code{TRUE}, the \code{qdrg}
#'     method for each class is registered instead.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' #--- If your package provides recover_data and emm_grid methods for class 'mymod',
#' #--- put something like this in your package code -- say in zzz.R:
#'   .onLoad <- function(libname, pkgname) {
#'     if (requireNamespace("emmeans", quietly = TRUE))
#'       emmeans::.emm_register("mymod", pkgname)
#'   }
#' }
.emm_register = function(classes, pkgname, qdrg = FALSE) {
    envir = asNamespace(pkgname)
    for (class in classes) {
        if(qdrg) 
            register_s3_method("emmeans", "recover_data", class, envir)
        else {
            register_s3_method("emmeans", "recover_data", class, envir)
            register_s3_method("emmeans", "emm_basis", class, envir)
        }
    }
}

