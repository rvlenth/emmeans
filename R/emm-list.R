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

# Methods for emm_list objects

# First, ehere is documentation for the emm_list class


#' The \code{emm_list} class
#' 
#' An \code{emm_list} object is simply a list of
#' \code{\link[=emmGrid-class]{emmGrid}} objects. Such a list is returned,
#' for example, by \code{\link{emmeans}} with a two-sided formula or a list as its
#' \code{specs} argument.
#' 
#' Methods for \code{emm_list} objects include \code{summary}, \code{CLD},
#' \code{coef}, \code{confint}, \code{contrast}, \code{pairs}, \code{plot},
#' \code{print}, and
#' \code{test}. These are all the same as those methods for \code{emmGrid}
#' objects, with an additional \code{which} argument (integer) to specify which 
#' members of the list to use. The default is \code{which = seq_along(object)};
#' i.e., the method is applied to every member of the \code{emm_list} object.
#' The exception is \code{plot}, where only the \code{which[1]}th element is 
#' plotted.
#' 
#' As an example,
#' to summarize a single member -- say the second one -- of an \code{emm_list}, 
#' one may use \code{summary(object, which = 2)}, but it is probably preferable 
#' to directly summarize it using \code{summary(object[[2]])}.
#'
#' @rdname emm_list-object
#' @name emm_list
#' @aliases CLD.emm_list
NULL


#' @export
#' @method str emm_list
str.emm_list = function(object, ...) {
    for(nm in names(object)) {
        cat(paste("$", nm, "\n", sep=""))
        str(object[[nm]])
        cat("\n")
    }
}


# summary.emm_list et all take an argument which that allows doing a subset
# Each returns a regular 'list'


#' @export
#' @method summary emm_list
summary.emm_list <- function(object, ..., which = seq_along(object))
    lapply(object[which], function(x) {
        if (inherits(x, "summary.emmGrid"))  x
        else summary.emmGrid(x, ...)
    })

#' @export
#' @method print emm_list
print.emm_list = function(x, ...) {
    print(summary(x, ...))
}

#' @export
#' @method contrast emm_list
contrast.emm_list = function(object, ... , which = seq_along(object)) {
    lapply(object[which], contrast, ...)
}

#' @export
#' @method pairs emm_list
pairs.emm_list = function(x, ..., which = seq_along(x)) {
    lapply(x[which], pairs, ...)
}

#' @export
#' @method test emm_list
test.emm_list = function(object, ..., which = seq_along(object)) {
    lapply(object[which], test, ...)
}

#' @export
#' @method confint emm_list
confint.emm_list = function(object, ..., which = seq_along(object)) {
    lapply(object[which], confint, ...)
}

#' @export
#' @method CLD emm_list
CLD.emm_list = function(object, ..., which = seq_along(object)) {
    if (length(which) > 1)
        warning("`CLD()` called with a list of ", length(which), " objects. ",
             "Only the first one was used.")
    CLD(object[[which[1]]], ...)
}

#' @export
#' @method coef emm_list
coef.emm_list = function(object, ..., which = seq_along(object)) {
    lapply(object[which], coef, ...)
}

# plot just plots one

#' @export
#' @method plot emm_list
plot.emm_list = function(x, ..., which = 1) {
    plot.emmGrid(x[[which[1]]], ...)
}

#' @export
#' @method rbind emm_list
rbind.emm_list = function(...) {
    do.call(rbind, ...)
}

#' @export
#' @method as.data.frame emm_list
as.data.frame.emm_list = function(x, ...) {
    as.data.frame(rbind(x, ...))
}

#' @export
#' @method as.list emm_list
as.list.emm_list = function(x, ...) {
    rtn = list()
    for (nm in names(x))
        rtn[[nm]] = as.list.emmGrid(x[[nm]])
    attr(rtn, "emm_list") = TRUE
    rtn
}

#' @export
#' @return \code{as.emm_list} returns an object of class \code{emm_list}.
#' 
#' @rdname as.emmGrid
#' @order 3
as.emm_list = function(object, ...) {
    if (is.null(attr(object, "emm_list")))
        as.emmGrid(object, ...)
    else
        lapply(object, as.emmGrid, ...)
}


