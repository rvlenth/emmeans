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

# Functions used only internally and of fairly general use
# More of these with specific uses only within a certain context remain there.

## Alternative to all.vars, but keeps vars like foo$x and foo[[1]] as-is
##   Passes ... to all.vars
#' @export
.all.vars = function(expr, retain = c("\\$", "\\[\\[", "\\]\\]"), ...) {
    if (is.null(expr))
        return(character(0))
    if (!inherits(expr, "formula")) {
        expr = try(eval(expr), silent = TRUE)
        if(inherits(expr, "try-error")) {
            return(character(0))
        }
    }
    repl = paste("_Av", seq_along(retain), "_", sep = "")
    for (i in seq_along(retain))
        expr = gsub(retain[i], repl[i], expr)
    subs = switch(length(expr), 1, c(1,2), c(2,1,3))
    vars = all.vars(as.formula(paste(expr[subs], collapse = "")), ...)
    retain = gsub("\\\\", "", retain)
    for (i in seq_along(retain))
        vars = gsub(repl[i], retain[i], vars)
    if(length(vars) == 0) vars = "1"   # no vars ---> intercept
    vars = trimws(unlist(strsplit(vars, "[*:+]")))
}

### parse a formula of the form lhs ~ rhs | by into a list
### of variable names in each part
### Returns character(0) for any missing pieces
.parse.by.formula = function(form) {
    allv = .all.vars(form)
    ridx = ifelse(length(form) == 2, 2, 3)
    allrhs = as.character(form)[ridx]
    allrhs = gsub("\\|", "+ .by. +", allrhs) # '|' --> '.by.'
    allrhs = .all.vars(stats::reformulate(allrhs))
    bidx = grep(".by.", allrhs, fixed = TRUE)
    if (length(bidx) == 0) { # no 'by' vars
        by = character(0)
        rhs = allrhs
    }
    else {
        rhs = allrhs[seq_len(bidx - 1)]
        by = setdiff(allrhs, c(rhs, ".by."))
    }
    lhs = setdiff(allv, allrhs)
    list(lhs = lhs, rhs = rhs, by = by)
}


# Utility to pick out the args that can be passed to a function
.args.for.fcn = function(fcn, args) {
    oknames = names(as.list(args(fcn)))
    mat = pmatch(names(args), oknames)
    args = args[!is.na(mat)]
    mat = mat[!is.na(mat)]
    names(args) = oknames[mat]
    args
}

# Create a list and give it class class.name
.cls.list <- function(class.name, ...) {
    result <- list(...)
    class(result) <- c(class.name, "list")
    result
}

### Not-so-damn-smart replacement of diag() that will 
### not be so quick to assume I want an identity matrix
### returns matrix(x) when x is a scalar
#' @export
.diag = function(x, nrow, ncol) {
    if(is.matrix(x))
        diag(x)
    else if((length(x) == 1) && missing(nrow) && missing(ncol)) 
        matrix(x)
    else 
        diag(x, nrow, ncol)
}

# Utility that returns TRUE if getOption("emmeans")[[opt]] is TRUE
.emmGrid.is.true = function(opt) {
    x = get_emm_option(opt)
    if (is.logical(x))  x
    else FALSE
}


# return list of row indexes in tbl for each combination of by
# tbl should be a data.frame
.find.by.rows = function(tbl, by) {
    if (is.null(by))
        return(list(seq_len(nrow(tbl))))
    if (any(is.na(match(by, names(tbl)))))
        stop("'by' variables are not all in the grid")    
    bylevs = tbl[ , by, drop = FALSE]
    by.id = do.call("paste", bylevs)
    uids = unique(by.id)
    result = lapply(uids, function(id) which(by.id == id))
    names(result) = uids
    result
}

# calculate the offset for the given grid
#' @export
.get.offset = function(terms, grid) {
    off.idx = attr(terms, "offset")
    offset = rep(0, nrow(grid))
    tvars = attr(terms, "variables")
    for (i in off.idx)
        offset = offset + eval(tvars[[i+1]], grid)
    offset
}

######################################################################
### Contributed by Jonathon Love, https://github.com/jonathon-love ###
######################################################################
# reformulate for us internally in emmeans
# same as stats::reformulate, except it surrounds term labels with backsticks
#
# RVL note: I renamed it .reformulate to avoid certain issues.
#   For example I need reformulate() sometimes to strip off function calls
#   and this .reformulate works quite differently.
#
.reformulate <- function (termlabels, response = NULL, intercept = TRUE)
{
    if (!is.character(termlabels) || !length(termlabels))
        stop("'termlabels' must be a character vector of length at least one")
    has.resp <- !is.null(response)
    termtext <- paste(if (has.resp)
        "response", "~", paste0("`", trimws(termlabels), "`", collapse = "+"), collapse = "")
    if (!intercept)
        termtext <- paste(termtext, "- 1")
    rval <- eval(parse(text = termtext, keep.source = FALSE)[[1L]])
    if (has.resp)
        rval[[2L]] <- if (is.character(response))
            as.symbol(response)
    else response
    environment(rval) <- parent.frame()
    rval
}

