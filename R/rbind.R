##############################################################################
#    Copyright (c) 2012-2025 Russell V. Lenth                                #
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

# rbind method for emmGrid objects
# ... and some other methods

#' Combine or subset \code{emmGrid} objects
#'
#' These functions provide methods for \code{\link[base:cbind]{rbind}} and
#' \code{\link[base:Extract]{[}} that may be used to combine \code{emmGrid} objects
#' together, or to extract a subset of cases. The primary reason for 
#' doing this would be to obtain multiplicity-adjusted results for smaller
#' or larger families of tests or confidence intervals. 
#' 
#' @param ... In \code{rbind}, object(s) of class \code{emmGrid} or \code{summary_emm}.
#'   In others, additional arguments passed to other methods
#' @param deparse.level (required but not used)
#' @param adjust Character value passed to \code{\link{update.emmGrid}}
#' 
#' @note \code{rbind} throws an error if there are incompatibilities in
#'   the objects' coefficients, covariance structures, etc. But they 
#'   are allowed to have different factors; a missing level \code{'.'}
#'   is added to factors as needed.
#' 
#' @note These functions generally reset \code{by.vars} to \code{NULL};
#' so if you want to keep any \dQuote{by} variables, you should follow-up
#' with \code{\link{update.emmGrid}}.
#'
#' @return A revised object of class \code{emmGrid}
#' @order 1
#' @method rbind emmGrid
#' @export
#' @examples
#' warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
#' warp.rg <- ref_grid(warp.lm)
#' 
#' # Do all pairwise comparisons within rows or within columns, 
#' # all considered as one faily of tests:
#' w.t <- pairs(emmeans(warp.rg, ~ wool | tension))
#' t.w <- pairs(emmeans(warp.rg, ~ tension | wool))
#' rbind(w.t, t.w, adjust = "mvt")
#' update(w.t + t.w, adjust = "fdr")  ## same as above except for adjustment
rbind.emmGrid = function(..., deparse.level = 1, adjust = "bonferroni") {
    objs = list(...)
    if (!all(sapply(objs, inherits, "emmGrid")))
        stop("All objects must inherit from 'emmGrid'")
    bhats = lapply(objs, function(o) o@bhat)
    bhat = bhats[[1]]
    if(!all(sapply(bhats, function(b) (length(b) == length(bhat)) 
                   && (sum((b - bhat)^2, na.rm = TRUE) == 0))))
        stop("All objects must have the same fixed effects")
    Vs = lapply(objs, function(o) o@V)
    V = Vs[[1]]
    if(!all(sapply(Vs, function(v) sum((v - V)^2) == 0)))
        stop("All objects must have the same covariances")
    obj = objs[[1]]
    linfcts = lapply(objs, function(o) o@linfct)
    obj@linfct = do.call(rbind, linfcts)
    bnms = unlist(lapply(objs, function(o) o@misc$by.vars))
    grids = lapply(objs, function(o) o@grid)
    gnms = unique(c(bnms, unlist(lapply(grids, names))))
    gnms = setdiff(gnms, c(".wgt.", ".offset.")) # exclude special names
    grid = data.frame(.tmp. = seq_len(n <- nrow(obj@linfct)))
    for (g in gnms)
        grid[[g]] = rep(".", n)
    grid[[".wgt."]] = grid[[".offset."]] = 0
    grid$.tmp. = NULL
    n.before = 0
    for (g in grids) {
        rows = n.before + seq_along(g[[1]])
        n.before = max(rows)
        for (nm in setdiff(names(g), c(".wgt.", ".offset.")))
            grid[rows, nm] = as.character(g[[nm]])
        if (!is.null(g$.wgt.)) grid[rows, ".wgt."] = g$.wgt.
        if (!is.null(g$.offset.)) grid[rows, ".offset."] = g$.offset.
    }
    if (all(grid$.wgt. == 0)) 
        grid$.wgt. = 1
    if (all(grid$.offset. == 0)) 
        grid$.offset. = NULL
    avgd.over = unique(unlist(lapply(objs, function(o) o@misc$avgd.over)))
    attr(avgd.over, "qualifier") = " some or all of"
    obj@grid = grid
    obj@levels = lapply(gnms, function(nm) unique(grid[[nm]]))
    names(obj@levels) = gnms
    obj@roles$predictors = setdiff(names(obj@levels), obj@roles$multresp)
    obj@misc$con.coef = obj@misc$orig.grid = obj@misc$.pairby = NULL
    update(obj, pri.vars = gnms, by.vars = NULL, adjust = adjust, initMesg = NULL,
           estType = "rbind",
           famSize = n,   # instead of round((1 + sqrt(1 + 8*n)) / 2, 3),
           avgd.over = avgd.over)
}
#' @rdname rbind.emmGrid
#' @order 4
#' @param e1,e2,x,object Objects of class \code{emmGrid}
#' @return The result of \code{e1 + e2} is the same as \code{rbind(e1, e2)}
#' @method + emmGrid
#' @export
"+.emmGrid" = function(e1, e2) {
    if(!is(e2, "emmGrid"))
        stop("'+.emmGrid' works only when all objects are class `emmGrid`", call. = FALSE)
    rbind(e1, e2)
}


### Subset a reference grid
# if drop = TRUE, the levels of factors are reduced
#' @rdname rbind.emmGrid
#' @order 5
#' @param i Integer vector of indexes
#' @param drop.levels Logical value. If \code{TRUE}, the \code{"levels"} slot in
#'   the returned object is updated to hold only the predictor levels that actually occur
#'   
#' @method [ emmGrid
#' @export
#' @examples
#'
#' # Show only 3 of the 6 cases
#' summary(warp.rg[c(2, 4, 5)])
"[.emmGrid" = function(x, i, adjust, drop.levels = TRUE, ...) {
    x@linfct = x@linfct[i, , drop = FALSE]
    x@grid = x@grid[i, , drop = FALSE]                  
    x = update(x, pri.vars = names(x@grid), famSize = length(i), estType = "[")
    x@misc$orig.grid = x@misc$con.coef = x@misc$.pairby = NULL
    x@misc$by.vars = NULL
    if(!missing(adjust))
        x@misc$adjust = adjust
    if(!is.null(disp <- x@misc$display))
        x@misc$display = disp[i]
    if (drop.levels) {
        for (nm in names(x@levels))
            x@levels[[nm]] = unique(x@grid[[nm]])
    }
    x
}

#' @rdname rbind.emmGrid
#' @order 9
#' @param subset logical expression indicating which rows of the grid to keep
#' @method subset emmGrid
#' @export
#' @examples
#' 
#' # After-the-fact 'at' specification
#' subset(warp.rg, wool == "A")  ## or warp.rg |> subset(wool == "A")
#' 
subset.emmGrid = function(x, subset, ...) {
    sel = eval(substitute(subset), envir = x@grid, enclos = parent.frame())
    x[sel, ...]
}

#' @rdname rbind.emmGrid
#' @order 6
#' @param n integer number of entries to include (or exclude if negative)
#' @importFrom utils head
#' @method head emmGrid
#' @export
head.emmGrid = function(x, n = 6, ...) {
    s = sign(n[1])
    n = min(abs(n[1]), nrow(x@grid))
    x[s * seq_len(n)]
}

#' @rdname rbind.emmGrid
#' @order 7
#' @importFrom utils tail
#' @method tail emmGrid
#' @export
tail.emmGrid = function(x, n = 6, ...) {
    s = sign(n[1])
    N = nrow(x@grid)
    n = min(abs(n[1]), N)
    x[s * (seq_len(n) - n + N)]
}

