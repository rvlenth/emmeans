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


# First, here is documentation for the emm_list class


#' The \code{emm_list} class
#' 
#' An \code{emm_list} object is simply a list of
#' \code{\link[=emmGrid-class]{emmGrid}} objects. Such a list is returned,
#' for example, by \code{\link{emmeans}} with a two-sided formula or a list as its
#' \code{specs} argument. Several methods for this class are provided, as detailed below.
#' 
#' @param object,x an object of class \code{emm_list}
#' @param ... additional arguments passed to corresponding \code{emmGrid} method
#' @param which integer vector specifying which elements to select.
#' 
#' @return a \code{list} of objects returned by the corresponding \code{emmGrid}
#'   method (thus, often, another \code{emm_list} object). However, if
#'   \code{which} has length 1, the one result is not wrapped in a list.
#' 
#' @rdname emm_list-object
#' @name emm_list
#' @order 1
NULL

#' @rdname emm_list-object
#' @order 99
#' 
#' @note \code{I_bet()} provides a default value for \code{which} noisily.
#' It is used in cases where users
#' most likely intended to call a method for an \code{emmGrid} object
#' rather than an \code{emm_list} object. 
#' An informative message is displayed
#' and then the \code{which}th element is returned. 
#'
#' @export
I_bet = function(which) {
    message("I bet you wanted to call this with just object[[", which, "]]")
    which
}

# Internal utility to noisily return one of an emm_list
.chk.list = function(object, which = 1) {
    if (inherits(object, "emm_list"))
        object = object[[I_bet(which)]]
    object
}

# My own lapply() function that drops when the dimension is 1
.lapply = function(...) {
    rtn = lapply(...)
    if (length(rtn) == 1)   rtn[[1]]
    else                    rtn
}


#' @export
#' @method str emm_list
#' @rdname emm_list-object
#' @order 2
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
#' @rdname emm_list-object
#' @order 3
summary.emm_list <- function(object, ..., which = seq_along(object))
    .lapply(object[which], function(x) {
        if (inherits(x, "summary.emmGrid"))  x
        else summary.emmGrid(x, ...)
    })

#' @export
#' @method print emm_list
#' @rdname emm_list-object
#' @order 4
#' @note No \code{export} option is provided for printing an \code{emm_list}
#' (see \code{\link{print.emmGrid}}). If you wish to export these objects, you 
#' must do so separately for each element in the list.
#'
print.emm_list = function(x, ...) {
    print(summary(x, ...))
}

#' @export
#' @method contrast emm_list
#' @rdname emm_list-object
#' @order 11
contrast.emm_list = function(object, ... , which = I_bet(1)) {
    .lapply(object[which], contrast, ...)
}

#' @export
#' @method pairs emm_list
#' @rdname emm_list-object
#' @order 12
pairs.emm_list = function(x, ..., which = I_bet(1)) {
    .lapply(x[which], pairs, ...)
}

#' @export
#' @method test emm_list
#' @rdname emm_list-object
#' @order 6
test.emm_list = function(object, ..., which = seq_along(object)) {
    .lapply(object[which], test, ...)
}

#' @export
#' @method confint emm_list
#' @rdname emm_list-object
#' @order 5
confint.emm_list = function(object, ..., which = seq_along(object)) {
    .lapply(object[which], confint, ...)
}

#' @export
#' @method coef emm_list
#' @rdname emm_list-object
#' @order 7
coef.emm_list = function(object, ..., which = I_bet(1)) {
    .lapply(object[which], coef, ...)
}

# plot just plots one

#' @export
#' @method plot emm_list
#' @rdname emm_list-object
#' @order 15
#' @note The \code{plot} method uses only the first element of \code{which}; the others are ignored.
plot.emm_list = function(x, ..., which = I_bet(1)) {
    plot.emmGrid(x[[which[1]]], ...)
}

#' @rdname rbind.emmGrid
#' @order 3
#' @param which Integer vector of subset of elements to use; if missing, all are combined
#' @return The \code{rbind} method for \code{emm_list} objects simply combines 
#' the \code{emmGrid} objects comprising the first element of \code{...}.
#' @export
#' @method rbind emm_list
#' @examples
#' 
#' ### Working with 'emm_list' objects
#' mod <- lm(conc ~ source + factor(percent), data = pigs)
#' all <- emmeans(mod, list(src = pairwise ~ source, pct = consec ~ percent))
#' rbind(all, which = c(2, 4), adjust = "mvt")
rbind.emm_list = function(..., which, adjust = "bonferroni") {
    elobj = list(...)[[1]]
    if(!missing(which))
        elobj = elobj[which]
    class(elobj) = c("emm_list", "list")
    update(do.call(rbind.emmGrid, elobj), adjust = adjust)
}

#' @export
#' @method as.data.frame emm_list
as.data.frame.emm_list = function(x, ...) {
    data.frame(rbind(x, ..., check.names = FALSE))
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


### Others we won't document

update.emm_list = function(object, ...)
    update.emmGrid(object[[I_bet(1)]])

predict.emm_list = function(object, ...)
    predict.emmGrid(object[[I_bet(1)]], ...)

vcov.emm_list = function(object, ...)
    vcov.emmGrid(object[[I_bet(1)]], ...)

xtable.emm_list = function(x, ...)
    xtable.emmGrid(x[[I_bet(1)]], ...)

