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

### Wrappers for those who want other (familiar) terminology: lsmeans 
                                                    #xxx  and pmmeans xxx#

### general-purpose wrapper for creating pmxxxxx and lsxxxxx functions
### use subst arg to specify e.g. "ls" #xxx or "pm"
.emwrap = function(emmfcn, subst, ...) {
    result = emmfcn(...)

        if (inherits(result, "emmGrid"))
        result = .sub.em(result, subst)
    else if(inherits(result, "emm_list")) {
        for (i in seq_along(result))
            result[[i]] = .sub.em(result[[i]], subst)
        names(result) = gsub("^em", subst, names(result))
    }
    result
}

# returns an updated emmGrid object with estName "em..." replaced by "xx..."
.sub.em = function(object, subst) {
    nm = object@misc$estName
    update(object, estName = gsub("^em", subst, nm))
}



### Exported implementations

# lsmeans family
#' Wrappers for alternative naming of EMMs
#' 
#' These are wrappers for \code{\link{emmeans}} and related functions to provide
#' backward compatibility, or for users who may prefer to
#' use other terminology than \dQuote{estimated marginal means} -- namely 
#' \dQuote{least-squares means}. These functions also provide the functionality
#' formerly provided by the \pkg{lsmeans} package, which is now just a front-end
#' for \pkg{emmeans}.
#' 
#' For each function with \code{ls}\emph{xxxx} in its name,
#' the same function named \code{em}\emph{xxxx} is called. Any estimator names or 
#' list items beginning with \dQuote{em} are replaced with \dQuote{ls} 
#' before the results are returned
#' 
#' @param ... Arguments passed to the corresponding \code{em}\emph{xxxx} function
#' 
#' @return The result of the call to \code{em}\emph{xxxx}, suitably modified.
#' @rdname wrappers
#' @aliases wrappers
#' @seealso \code{\link{emmeans}}, \code{\link{emtrends}}, \code{\link{emmip}},
#'          \code{\link{emm}}, \code{\link{emmobj}}, \code{\link{emm_options}},
#'          \code{\link{get_emm_option}}
#' @export
#' @examples
#' pigs.lm <- lm(log(conc) ~ source + factor(percent), data = pigs)
#' lsmeans(pigs.lm, "source")
lsmeans = function(...)
    .emwrap(emmeans, subst = "ls", ...)


#' @rdname wrappers
#' @export
lstrends = function(...)
    .emwrap(emtrends, subst = "ls", ...)


#' @rdname wrappers
#' @export
lsmip = function(...)
    emmip(...)


#' @rdname wrappers
#' @export
lsm = function(...)
    emm(...)


#' @rdname wrappers
#' @export
lsmobj = function(...)
    .emwrap(emmobj, subst = "ls", ...)

#' @rdname wrappers
#' @export
lsm.options = function(...) {
    .Deprecated("emm_options")
    args = list(...)
    nms = names(args)
    nms = gsub("ref.grid", "ref_grid", nms)
    nms = gsub("lsmeans", "emmeans", nms)
    names(args) = nms
    do.call(emm_options, args)
}

#' @rdname wrappers
#' @param x Character name of desired option
#' @param default default value to return if \code{x} not found
#' 
#' @return \code{get.lsm.option} and \code{lsm.options} remap options from
#'   and to corresponding options in the \pkg{emmeans} options system.
#' @export
get.lsm.option = function(x, default = emm_defaults[[x]]) {
    .Deprecated("get_emm_option")
    if(x == "ref.grid") x = "ref_grid"
    if(x == "lsmeans") x = "emmeans"
    get_emm_option(x, default = default)
}
    