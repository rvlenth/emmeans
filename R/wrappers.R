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

### Wrappers for those who want other (familiar) terminology: lsmeans and pmmeans

### general-purpose wrapper for creating pmxxxxx and lsxxxxx functions
### use subst arg to specify e.g. "ls" or "pm"
.emwrap = function(emmfcn, subst, ...) {
    result = emmfcn(...)

        if (inherits(result, "emm"))
        result = .sub.em(result, subst)
    else if(inherits(result, "emm_list")) {
        for (i in seq_along(result))
            result[[i]] = .sub.em(result[[i]], subst)
        names(result) = gsub("^em", subst, names(result))
    }
    result
}

# returns an updated emm object with estName "em..." replaced by "xx..."
.sub.em = function(object, subst) {
    nm = object@misc$estName
    update(object, estName = gsub("^em", subst, nm))
}



### Exported implementations

# lsmeans family
#' Wrappers for alternative naming of EMMs
#' 
#' These are wrappers for \code{\link{emmeans}} and realated functions to provide
#' backward compatibility, or for users who may prefer to
#' use other terminology than \dQuote{estimated marginal means} -- namely 
#' \dQuote{least-squares means} or \dQuote{predicted marginal means}.
#' 
#' For each function with \code{ls}\emph{xxxx} or \code{pm}\emph{xxxx} in its name,
#' the same function named \code{em}\emph{xxxx} is called. Any estimator names or 
#' list items beginning with \dQuote{em} are replaced with \dQuote{ls} or 
#' \dQuote{pm} before the results are returned
#' 
#' @param ... Arguments passed to the corresponding \code{em}\emph{xxxx} function
#' 
#' @return The result of the call to \code{em}\emph{xxxx}, suitably modified
#' @rdname wrappers
#' @aliases wrappers
#' @export
#' @examples
#' pigs.lm <- lm(log(conc) ~ source + factor(percent), data = pigs)
#' lsmeans(pigs.lm, "source")
lsmeans = function(...)
    .emwrap(emmeans, subst = "ls", ...)

#' @rdname wrappers
#' @export
pmmeans = function(...)
    .emwrap(emmeans, subst = "pm", ...)



#' @rdname wrappers
#' @export
lstrends = function(...)
    .emwrap(emtrends, subst = "ls", ...)

#' @rdname wrappers
#' @export
pmtrends = function(...)
    .emwrap(emtrends, subst = "pm", ...)



#' @rdname wrappers
#' @export
lsmip = function(...)
    emmip(...)

#' @rdname wrappers
#' @export
pmmip = function(...)
    emmip(...)



#' @rdname wrappers
#' @export
lsm = function(...)
    emm(...)

#' @rdname wrappers
#' @export
pmm = function(...)
    emm(...)



#' @rdname wrappers
#' @export
lsmobj = function(...)
    .emwrap(emmobj, subst = "ls", ...)

#' @rdname wrappers
#' @export
pmmobj = function(...)
    .emwrap(emmobj, subst = "pm", ...)

