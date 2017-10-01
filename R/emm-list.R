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
        if (inherits(x, "summary.emm"))  x
        else summary.emm(x, ...)
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
#' @method confint emm_list
cld.emm_list = function(object, ..., which = seq_along(object)) {
    lapply(object[which], cld, ...)
}

#' @export
#' @method coef emm_list
coef.emm_list = function(object, ..., which = seq_along(object)) {
    lapply(object[which], coef, ...)
}


