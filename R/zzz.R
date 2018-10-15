##############################################################################
#    Copyright (c) 2012-2016 Russell V. Lenth                                #
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

# NOTE: Excluded from documentation
# Custom Vignette format
# 
# This is used to format HTML vignettes the way its developer wants them.
#
# @param ... Arguments passed to \code{rmarkdown::html_document}
#
# @return R Markdown format used by \code{rmarkdown::render}
#' @export
.emm_vignette = function(css = system.file("doc", "clean-simple.css", package = "emmeans"), 
                         highlight = NULL,  ...) {
    rmarkdown::html_document(theme = NULL, highlight = highlight,
                             fig_width = 3, fig_height = 3, 
                             css = css, ...)
}


### Dynamic registration of S3 methods
# Code borrowed from hms pkg. I omitted some type checks etc. because
# this is only for internal use and I solemnly promise to behave myself.
register_s3_method = function(pkg, generic, class) {
    fun = get(paste0(generic, ".", class), envir = parent.frame())
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
    if (requireNamespace("coda", quietly = TRUE)) {
        register_s3_method("coda", "as.mcmc", "emmGrid")
        register_s3_method("coda", "as.mcmc.list", "emmGrid")
    }
    if (requireNamespace("multcomp", quietly = TRUE)) {
        register_s3_method("multcomp", "glht", "emmlf")
        register_s3_method("multcomp", "glht", "emmGrid")
        register_s3_method("multcomp", "cld", "emmGrid")
        register_s3_method("multcomp", "modelparm", "emmwrap")
    }
    # regardless...
    register_s3_method("emmeans", "cld", "emmGrid")
}

# .onAttach = function(libname, pkgname) {
#     packageStartupMessage (
#         "NOTE: As of emmeans versions > 1.2.3,\n",
#         "      The 'cld' function will be deprecated in favor of 'CLD'.\n",
#         "      You may use 'cld' only if you have package:multcomp attached."
#     )
# }

