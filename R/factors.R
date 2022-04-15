##############################################################################
#    Copyright (c) 2012-2022 Russell V. Lenth                                #
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



# Combine two or more factors (named in facs) into a new factor named newname
#' Manipulate factors in a reference grid
#' 
#' These functions manipulate the velels of factors comprising a reference
#' grid by combining factor levels, splitting a factor's levels into 
#' combinations of newly-defined factors, or by creating a grouping factor in which 
#' factor(s) levels are nested.
#' 
#' 
#' @param object An object of class \code{emmGrid}
#' @param facs Character vector. The names of the factors to combine
#' @param newname Character value. The name of the new factor
#' @param drop Logical value. If \code{TRUE}, any levels of the new factor that 
#'     are dropped if all occurrences in the newly reconstructed object have 
#'     weight zero. If \code{FALSE}, all levels are retained. 
#'     (This argument is ignored if there is no \code{.wgt.} column
#'     in \code{object@grid}.)
#' @rdname manip-factors
#'
#' @return A modified object of class \code{emmGrid}
#'
#' @section The \code{comb_facs} function:
#' \code{comb_facs} combines the levels of factors into a single factor
#' in the reference grid (similar to \code{\link{interaction}}). This new factor
#' replaces the factors that comprise it.
#' 
#' \emph{Additional note:}
#' The choice of whether to drop levels or not can make a profound difference.
#' If the goal is to combine factors for use in \code{joint_tests}, we advise \emph{against} 
#' \code{drop = TRUE} because that might change the weights used in deriving marginal means.
#' If combining factors in a nested structure, dropping unused cases can considerably reduce 
#' the storage required.
#'  
#'
#' @export
#'
#' @examples
#' mtcars.lm <- lm(mpg ~ factor(vs)+factor(cyl)*factor(gear), data = mtcars)
#' (v.c.g <- ref_grid(mtcars.lm))
#' (v.cg <- comb_facs(v.c.g, c("cyl", "gear")))
#'   
#' # One use is obtaining a single test for the joint contributions of two factors:
#' joint_tests(v.c.g)
#' 
#' joint_tests(v.cg)
#' 
#' # undo the 'comb_facs' operation:
#' split_fac(v.cg, "cyl.gear", list(cyl = c(4, 6, 8), gear = 3:5))
#' 
comb_facs = function(object, facs, newname = paste(facs, collapse = "."),
                     drop = FALSE) {
    if((length(facs)  < 1))
        stop("No factors have been specified to combine")
    levs = object@levels
    if(any(sapply(facs, function(x) !(x %in% names(levs)))))
        stop("Unknown factor(s)")
    if (newname %in% facs)
        stop("newname of '", newname, "' cannot be used. Specify something else")
    grid = object@grid
    idx = sapply(facs, function(x) which(names(levs) == x))
    levs[[newname]] = do.call(paste, c(do.call(expand.grid, levs[idx]), sep = ":"))
    grid[[newname]] = do.call(paste, c(grid[idx], sep = ":"))
    levs[idx] = grid[idx] = NULL
    if (drop && !is.null(w <- grid$.wgt.)) {
        rows = lapply(levs[[newname]], function(nm) which(grid[[newname]] == nm))
        maxw = sapply(rows, function(i) max(w[i]))
        zw = which(maxw == 0)
        if(length(zw) > 0) {
            r = unlist(rows[zw])
            levs[[newname]] = levs[[newname]][-r]
            grid = grid[-r, , drop = FALSE]
            object@linfct = object@linfct[-r, , drop = FALSE]
            if(!is.null(object@misc$display))
                object@misc$display = object@misc$display[-r]
        }
    }
    if(!is.null(nests <- object@model.info$nesting)) {
        # which of facs are nested?
        m = match(facs, names(nests), nomatch = 0)
        if (length(m[m > 0]) == length(facs)) { # all are nested, find common nests
            nst = nests[[m[1]]]
            for (mm in m[-1])
                nst = intersect(nst, nests[[mm]])
            if (length(nst) > 0)
                nests[newname] = nst
        }
        nests[m] = NULL # zap out all facs that are nested
        
        # look for where other factors are nested in facs factors
        for (nm in names(nests)) {
            # factors nested in any facs are also nested in the new factor
            if (any(facs %in% nests[[nm]])) {
                m = match(facs, nests[[nm]], nomatch = 0)
                nests[[nm]] = c(nests[[nm]][-m], newname)
            }
        }
        if (length(nests) == 0) nests = NULL
        object@model.info$nesting = nests
    }
    
    object@levels = levs
    object@grid = grid
    object@roles$predictors = c(object@roles$predictors[-idx], newname)
    
    ord = .std.order(grid, levs)
    object[ord]
}

#' @rdname manip-factors
#' @param fac The name of a factor that is part of the grid in \code{object}
#' @param newfacs A named list with the names of new factors
#'   and their levels. The names must not already exist in the object,
#'   and the product of the lengths of the levels must equal the number
#'   of levels of \code{fac}.
#'   
#' @section The \code{split_fac} function:
#' The levels in \code{newfacs} are expanded via \code{\link{expand.grid}} into
#' combinations of levels, and the factor \code{fac} is replaced by those 
#' factor combinations. Unlike \code{add_grouping}, this creates a crossed, 
#' rather than a nested structure. Note that the order of factor combinations
#' is systematic with the levels of first factor in \code{newfacs} varying 
#' the fastest; and those factor combinations are assigned respectively
#' to the levels of \code{fac} as displayed in \code{str(object)}.
#'
#' @export
#'
#' @examples
#' IS.glm <- glm(count ~ spray, data = InsectSprays, family = poisson)
#' IS.emm <- emmeans(IS.glm, "spray")
#' IS.new <- split_fac(IS.emm, "spray", list(A = 1:2, B = c("low", "med", "hi")))
#' str(IS.new)
#'
split_fac = function(object, fac, newfacs) {
    if (!(fac %in% names(object@grid)))
        stop("The factor '", fac, "' is not in the reference grid")
    newg = expand.grid(newfacs)
    ref = object@levels[[fac]]
    if(nrow(newg) != length(ref))
        stop("Mismatch between 'fac' levels (", length(ref),
             ") and new factor combinations (", nrow(newg), ")")
    iy = as.numeric(factor(object@grid[[fac]], levels = ref))
    repl = newg[iy, , drop = FALSE]
    idx = which(names(object@levels) == fac)
    i = seq_along(object@grid)
    object@grid = cbind(object@grid[i < idx], repl, object@grid[i>idx])
    i = seq_along(object@levels)
    object@levels = c(object@levels[i < idx], newfacs, object@levels[i > idx])
    object@roles$predictors = object@misc$pri.vars = names(object@levels)
    object@misc$by.vars = NULL
    object
}


### Create a grouping factor and add it to a ref grid

#' @rdname manip-factors
#'
#' @param newname Character name of grouping factor to add (different from any
#'   existing factor in the grid)
#' @param refname Character name(s) of the reference factor(s)
#' @param newlevs Character vector or factor of the same length as that of the (combined) levels for 
#'   \code{refname}. The grouping factor \code{newname} will have the unique
#'   values of \code{newlevs} as its levels. The order of levels in \code{newlevs}
#'   is the same as the order of the level combinations produced by 
#'   \code{\link{expand.grid}} applied to the levels of \code{refname} -- that is, the
#'   first factor's levels change the fastest and the last one's vary the slowest.
#'
#' @section The \code{add_grouping} function:
#' This function adds a grouping factor to an existing reference grid or other 
#' \code{emmGrid} object, such that the levels of one or more existing factors (call them the
#' reference factors) are mapped to a smaller number of levels of the new
#' grouping factor. The reference factors are then nested in a 
#' new grouping factor named \code{newname}, and a new nesting structure
#' \code{refname \%in\% newname}.
#' This facilitates obtaining marginal means of the grouping factor, and 
#' contrasts thereof.
#' 
#' 
#' \emph{Additional notes:} By default, the levels of \code{newname} will be ordered
#'   alphabetically. To dictate a different ordering of levels, supply 
#'   \code{newlevs} as a \code{factor} having its levels in the desired order.
#'   
#' When \code{refname} specifies more than one factor, this can
#'   fundamentally (and permanently) change what is meant by the levels of those
#'   individual factors. For instance, in the \code{gwrg} example below, there
#'   are two levels of \code{wool} nested in each \code{prod}; and that implies
#'   that we now regard these as four different kinds of wool. Similarly, there
#'   are five different tensions (L, M, H in prod 1, and L, M in prod 2).
#'   
#' @export
#'
#' @examples
#' fiber.lm <- lm(strength ~ diameter + machine, data = fiber)
#' ( frg <- ref_grid(fiber.lm) )
#' 
#' # Suppose the machines are two different brands
#' brands <- factor(c("FiberPro", "FiberPro", "Acme"), levels = c("FiberPro", "Acme"))
#' ( gfrg <- add_grouping(frg, "brand", "machine", brands) )
#' 
#' emmeans(gfrg, "machine")
#' 
#' emmeans(gfrg, "brand")
#' 
#' ### More than one reference factor
#' warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
#' gwrg <- add_grouping(ref_grid(warp.lm), 
#'     "prod",  c("tension", "wool"),  c(2, 1, 1,  1, 2, 1))
#'         # level combinations:         LA MA HA  LB MB HB
#' 
#' emmeans(gwrg, ~ wool * tension)   # some NAs due to impossible combinations
#' 
#' emmeans(gwrg, "prod")
#' 
add_grouping = function(object, newname, refname, newlevs) {
    if(!is.null(object@model.info$nesting[[refname]]))
        stop("'", refname, "' is already nested in another factor; cannot re-group it")
    if(newname %in% object@roles$predictors)
        stop("'", newname, "' is already the name of an existing predictor")
    rlevs = do.call(paste, do.call(expand.grid, object@levels[refname]))
    if (length(newlevs) != length(rlevs))
        stop("Length of 'newlevs' doesn't match # levels of '", refname, "'")
    newlevs = factor(newlevs)
    glevs = levels(newlevs)
    k = length(glevs)
    
    one = matrix(1, nrow = k, ncol = 1)
    object@linfct = kronecker(one, object@linfct)
    object@levels[[newname]] = glevs
    object@roles$predictors = c(object@roles$predictors, newname)
    
    ref = do.call(paste, object@grid[refname]) # obs levels of rlevs
    wgt = object@grid$.wgt.
    if (is.null(wgt)) wgt = rep(1, nrow(object@grid))
    offset = object@grid$.offset.
    ogrid = object@grid[setdiff(names(object@grid), c(".wgt.", ".offset."))]
    grid = data.frame()
    valid = logical(0) # flag for rows that make sense
    for (i in 1:k) {
        g = ogrid
        g[[newname]] = glevs[i]
        g$.wgt. = wgt
        g$.offset. = offset
        grid = rbind(grid, g)
        alevs = rlevs[newlevs == glevs[i]]
        valid = c(valid, ref %in% alevs)
    }
    # screen out invalid rows
    grid[!valid, ".wgt."] = 0
    object@linfct[!valid, ] = NaN
    object@misc$pri.vars = c(object@misc$pri.vars, newname)
    if(is.null(disp <- object@misc$display))
        object@misc$display = valid
    else
        object@misc$display = disp & valid
    object@grid = grid
    
    # update nesting structure
    nesting = object@model.info$nesting
    if (is.null(nesting))
        nesting = list()
    for (nm in names(nesting))
        if (any(refname %in% nesting[[nm]]))
            nesting[[nm]] = c(nesting[[nm]], newname)
    for (nm in refname)
        nesting[[nm]] = newname   ### ??? should it be c(nesting[[nm]], newname)
    object@model.info$nesting = nesting
    
    object
}

