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

# Methods for emm_list objects


# First, here is documentation for the emm_list class


#' The \code{emm_list} class
#' 
#' An \code{emm_list} object is simply a list of
#' \code{\link[=emmGrid-class]{emmGrid}} objects. Such a list is returned,
#' for example, by \code{\link{emmeans}} with a two-sided formula or a list as its
#' \code{specs} argument. Several methods for this class are provided, as detailed below.
#' Typically, these methods just quietly do the same thing as their \code{emmGrid}
#' methods, using the first element of the list. You can specify \code{which}
#' to select a different element, or just run the corresponding \code{emmGrid}
#' method on \code{object[[k]]}.
#' 
#' @param object,x an object of class \code{emm_list}
#' @param ... additional arguments passed to corresponding \code{emmGrid}
#'   method. In addition, the user may include a logical argument \code{drop}
#'   that is akin to \code{\link{drop}} and the argument of the same name in
#'   subscripting matrices and data frames. When \code{drop} is \code{TRUE} (the
#'   default), then when the result is a \code{list} of length 1, the
#'   \code{list} structure is removed.
#' @param which integer vector specifying which elements to select;
#'   if \code{NULL},
#'   we try to guess which elements make sense. Usually, this is all elements having 
#'   names that start with \sQuote{em} or \sQuote{ls},
#'   or the first element if no matches are found. However, in \code{coef.emm_list},
#'   these are the ones we \emph{exclude}.
#' 
#' @return a \code{list} of objects returned by the corresponding \code{emmGrid}
#'   method (thus, often, another \code{emm_list} object). However, if
#'   \code{which} has length 1, the one result is not wrapped in a list.
#' 
#' @rdname emm_list-object
#' @name emm_list
#' @order 1
#' @examples
#' mod <- lm(conc ~ source, data = pigs)
#' obj <- emmeans(mod, pairwise ~ source)
#' 
#' linfct(obj)
#' 
#' coef(obj)     # done only for the contrasts
#' 
#' contrast(obj, "consec")  # done only for the means
#' 
#' contrast(obj, "eff", drop = FALSE)   # kept as a list
NULL

# Internal utility to noisily return one of an emm_list
# Call with ... argument so that message is suppressed if we specify 'which'
.chk.list = function(object, which, ...) {
    if (inherits(object, "emm_list")) {
        if (missing(which))
            which = 1
        object = object[[which]]
    }
    object
}

# My own lapply() function that drops when the dimension is 1
.lapply = function(X, ..., drop = TRUE) {
    oldClass(X) = "list"
    rtn = lapply(X, ...)
    if (drop && (length(rtn) == 1))   {
        rtn = rtn[[1]]
    }
    else if(length(rtn) != 0)                   {
        cls = ifelse(class(rtn[[1]])[1] == "emmGrid", "emm_list", "summary_eml")
        class(rtn) = c(cls, "list")
    }
    rtn
}


#' @export
#' @method str emm_list
#' @rdname emm_list-object
#' @order 12
str.emm_list = function(object, ...) {
    for(nm in names(object)) {
        cat(paste("$", nm, "\n", sep=""))
        str(object[[nm]])
        cat("\n")
    }
}


# summary.emm_list et al take an argument 'which' that allows doing a subset
# Each returns a regular 'list'


#' @export
#' @method summary emm_list
#' @return \code{summary.emm_list} returns an object
#' of class \code{summary_eml}, which is a list of \code{summary_emm}
#' objects.
#' @rdname emm_list-object
#' @order 13
summary.emm_list <- function(object, ..., which = seq_along(object)) {
    # .lapply(object[which], function(x) {
    #     if (inherits(x, "summary_emm"))  x
    #     else summary.emmGrid(x, ...)
    # })
    if(length(which) == 1)
        summary.emmGrid(object[[which]], ...)
    else
        .lapply(object[which], \(x) summary.emmGrid(x, ...))
}

#' @export
summary.summary_eml = function(object, ...) object

#' @export
#' @rdname emm_list-object
#' @order 25
#' @param row.names,optional Required arguments of \code{as.data.frame}, ignored
#' @method as.data.frame summary_eml
as.data.frame.summary_eml = function(x, row.names = NULL, optional = FALSE, which, ...) {
    rbind(x, which = which)
}

#' @export
print.summary_eml = function(x, ...) {
    attr(x, "class") = NULL
    print(x)
}

#' @export
#' @method print emm_list
#' @rdname emm_list-object
#' @order 14
#' @note No \code{export} option is provided for printing an \code{emm_list}
#' (see \code{\link{print.emmGrid}}). If you wish to export these objects, you 
#' must do so separately for each element in the list.
#'
print.emm_list = function(x, ...) {
    print(summary(x, ...))
}

### Utility fcn to identify which elements' names start with "em" or "lm"
### this only operates if which == NULL
.guess.which = guess.which = function(object, which) {
    if(!is.null(which))
        return(which)
    rtn = which(substr(paste0(names(object), "xx"), 1, 2) %in% c("em", "ls"))
    if(length(rtn) == 0)
        rtn = 1
    rtn
}

#' @export
#' @method contrast emm_list
#' @rdname emm_list-object
#' @order 3
contrast.emm_list = function(object, ... , which = NULL) {
    which = .guess.which(object, which)
    rtn = .lapply(object[which], contrast, ...)
    if(is.list(rtn))
        names(rtn) = paste("contrasts of", names(rtn))
    rtn
}

#' @export
#' @method pairs emm_list
#' @rdname emm_list-object
#' @order 4
pairs.emm_list = function(x, ..., which = NULL) {
    which = .guess.which(x, which)
    rtn = .lapply(x[which], pairs, ...)
    names(rtn) = paste("comparisons of", names(rtn))
    rtn
}

#' @export
#' @method test emm_list
#' @rdname emm_list-object
#' @order 6
test.emm_list = function(object, ..., which = seq_along(object)) {
    which = .guess.which(object, which)
    .lapply(object[which], test, ...)
}

#' @export
#' @method confint emm_list
#' @rdname emm_list-object
#' @order 7
confint.emm_list = function(object, ..., which = seq_along(object)) {
    which = .guess.which(object, which)
    .lapply(object[which], confint, ...)
}

#' @export
#' @method coef emm_list
#' @rdname emm_list-object
#' @order 9
coef.emm_list = function(object, ..., which = NULL) {
    which = - .guess.which(object, which)
    .lapply(object[which], coef, ...)
}

#' @export
#' @method linfct emm_list
#' @rdname emm_list-object
#' @order 11
linfct.emm_list = function(object, ..., which = seq_along(object)) {
    .lapply(object[which], \(x) attr(x, "linfct"), ...)
}


#' @export
#' @method plot emm_list
#' @rdname emm_list-object
#' @order 8
#' @note The \code{plot} method uses only the first element of \code{which}; the others are ignored.
plot.emm_list = function(x, ..., which = 1) {
    which = .guess.which(x, which)
    plot.emmGrid(x[[which[1]]], ...)
}

#' @rdname rbind.emmGrid
#' @order 23
#' @param which Integer vector of subset of elements to use;
#' if missing, we use all elements. 
#' @return The \code{rbind} method for \code{emm_list} objects simply combines 
#' the \code{emmGrid} objects comprising the first element of \code{...}.
#' Note that the returned object is not yet summarized, so any \code{adjust}
#' parameters apply to the combined \code{emmGrid}.
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

#' @rdname rbind.emmGrid
#' @order 33
#' @return The \code{rbind} method for \code{summary_emm} objects (or a list thereof)
#' returns a single \code{summary_emm} object. This combined object
#' \emph{preserves} any adjusted P values or confidence limits in the
#' original summaries, since those quantities have already been computed.
#' @export
#' @method rbind summary_emm
rbind.summary_emm = function(..., which) {
    slobj = list(...)
    if(!all(sapply(slobj, \(z) inherits(z, "data.frame")))) {
        # workaround to make tern.gee::lsmeans() work
        slobj = lapply(slobj, \(z) if(inherits(z, "data.frame")) data.frame(z) else z)
       return(do.call("rbind", slobj))
    }
    rbind.summary_eml(slobj, which = which)
}

#' 
#' @export
#' @method rbind summary_eml
rbind.summary_eml = function(..., which) {
    x = list(...)[[1]]
    if(!missing(which))
        x = x[which]
    nms.lst = lapply(x, names)
    bys = unique(do.call(c, lapply(x, \(z) attr(z, "by.vars"))))
    pris = unique(do.call(c, lapply(x, \(z) attr(z, "pri.vars"))))
    if (length(x) == 1) {
        attr(x[[1]], "pri.vars") = c(pris, bys)
        attr(x[[1]], "by.vars") = NULL
        return (x[[1]])
    }
    nms = pris = union(bys, pris)
    for (n in nms.lst)
        nms = union(nms, n)
    nums = setdiff(nms, pris)  # numeric columns
    xx = lapply(x, function(df) {
        d = data.frame(matrix(".", nrow = nrow(df), ncol = length(nms),
                       dimnames = list(NULL, nms)))
        d[, nums] = NA
        d[, names(df)] = df
        d
    })
    rtn = do.call("rbind", xx)
    row.names(rtn) = NULL
    class(rtn) = c("summary_emm", "data.frame")
    attr(rtn, "pri.vars") = pris
    attr(rtn, "estName") = attr(x[[1]], "estName")
    mesg = otr.mesg = attr(x[[1]], "mesg")
    for (i in 2:length(x)) {
        mesg = intersect(mesg, attr(x[[i]], "mesg"))
        otr.mesg = union(otr.mesg, attr(x[[i]], "mesg"))
    }
    otr.mesg = setdiff(otr.mesg, mesg)
    if (length(otr.mesg) > 0) 
        mesg = c(mesg, "The following messages apply only to some rows:",
                 paste("*", otr.mesg))
    attr(rtn, "mesg") = mesg
    rtn
}

#' @export
#' @rdname emm_list-object
#' @order 24
#' @return The \code{as.data.frame} methods return a single data frame via
#' \code{as.data.frame(rbind(x))}.
#' See also \code{\link{rbind.emm_list}} and \code{\link{as.data.frame.emmGrid}}
#' @method as.data.frame emm_list
as.data.frame.emm_list = function(x, ...) {
    if (length(x) > 1)
        warning("Note: 'as.data.frame' has combined your ", length(x), " sets of results into one object,\n",
                "and this affects things like adjusted P values. Refer to the annotations.")
    as.data.frame(rbind(x, ..., check.names = FALSE))
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
#' @exportS3Method update emm_list
update.emm_list = function(object, ...)
    update.emmGrid(object[[1]])

#' @exportS3Method predict emm_list
predict.emm_list = function(object, ...)
    predict.emmGrid(object[[1]], ...)

#' @exportS3Method vcov emm_list
vcov.emm_list = function(object, ...)
    vcov.emmGrid(object[[1]], ...)

xtable.emm_list = function(x, ...)
    xtable.emmGrid(x[[1]], ...)

