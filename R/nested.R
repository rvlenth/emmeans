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

# Code supporting nested models

# This code replies on nested structures specified in a named list like
#    list(a = "b", c = c("d", "e"))
# ... to denote a %in% b, c %in% d*e

### Create a grouping factor and add it to a ref grid

#' Add a grouping factor
#' 
#' This function adds a grouping factor to an existing reference grid or other 
#' \code{emmGrid} object, such that the levels of an existing factor (call it the
#' reference factor) are mapped to a smaller number of levels of the new
#' grouping factor. The reference factor is then nested in the grouping factor. 
#' This facilitates obtaining marginal means of the grouping factor, and 
#' contrasts thereof.
#'
#' @param object An \code{emmGrid} object
#' @param newname Character name of grouping factor to add (different from any
#'   existing factor in the grid)
#' @param refname Character name of the reference factor
#' @param newlevs Character vector or factor of the same length as that of the levels for 
#'   \code{refname}. The grouping factor \code{newname} will have the unique
#'   values of \code{newlevs} as its levels.
#'
#' @return A revised \code{emmGrid} object having an additional factor named 
#'   \code{newname}, and a new nesting structure \code{refname \%in\% newname}
#' 
#' @note By default, the levels of \code{newname} will be ordered
#'   alphabetically. To dictate a different ordering of levels, supply 
#'   \code{newlevs} as a \code{factor} having its levels in the required order.
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
add_grouping = function(object, newname, refname, newlevs) {
    if(!is.null(object@model.info$nesting[[refname]]))
        stop("'", refname, "' is already nested in another factor; cannot re-group it")
    if(newname %in% object@roles$predictors)
        stop("'", newname, "' is already the name of an existing predictor")
    rlevs = object@levels[[refname]]
    if (length(newlevs) != length(rlevs))
        stop("Length of 'newlevs' doesn't match # levels of '", refname, "'")
    newlevs = factor(newlevs)
    glevs = levels(newlevs)
    k = length(glevs)
    
    one = matrix(1, nrow = k, ncol = 1)
    object@linfct = kronecker(one, object@linfct)
    object@levels[[newname]] = glevs
    object@roles$predictors = c(object@roles$predictors, newname)
    
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
        valid = c(valid, g[[refname]] %in% alevs)
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
        if (refname %in% nesting[[nm]])
            nesting[[nm]] = c(nesting[[nm]], newname)
    nesting[[refname]] = newname
    object@model.info$nesting = nesting
    
    object
}



### ----- Rest of this file is used only internally ---------

# Internal function to deal with nested structures. 
#   rgobj        -- an emmGrid object
#   specs, ...   -- arguments for emmeans
#   nesting      -- a named list of nesting info
# This function works by subsetting rgobj as needed, and applying emmeans
# to each subsetted object
# This is a servant to emmeans.character.emmGrid, so we can assume specs is character
.nested_emm = function(rgobj, specs, by = NULL, ..., nesting) {
    # # Trap something not supported for these... This doesn't work
    # dots = list(...)
    # if("weights" %in% dots)
    #     if(!is.na(pmatch(dots$weights, "show.levels")))
    #         stop('weights = "show.levels" is not supported for nested models.')

    orig.by = by   # save original 'by' vars
    #### Two issues to worry about....
    # (1) specs contains nested factors. We need to include their grouping factors
    xspecs = intersect(union(specs, by), names(nesting))
    if (length(xspecs) > 0) {
        xgrps = unlist(nesting[xspecs])
        specs = union(union(xspecs, xgrps), specs)  # expanded specs with flagged ones first
        by = setdiff(by, xspecs) # can't use nested factors for grouping
    }
    # (2) If we average over any nested factors, we need to do it separately
    avg.over = setdiff(names(rgobj@levels), union(specs, by))
    afacs = intersect(names(nesting), avg.over) ### DUH!names(nesting)[names(nesting) %in% avg.over]
    rgobj@misc$display = NULL  ## suppress warning messages from emmeans
    
    if (length(afacs) == 0)  { # no nesting issues; just use emmeans
        result = emmeans(rgobj, specs, by = by, ...)
    }
    else { # we need to handle each group separately
        sz = sapply(afacs, function(nm) length(nesting[[nm]]))
        # use highest-order one first: potentially, we end up calling this recursively
        afac = afacs[rev(order(sz))][1] 
        otrs = setdiff(afacs, afac)   # other factors than afac
        grpfacs = union(nesting[[afac]], otrs)
        gspecs = union(specs, union(by, grpfacs))
        grpids = as.character(interaction(rgobj@grid[, grpfacs]))
        grps = do.call(expand.grid, rgobj@levels[grpfacs])  # all combinations of group factors
        result = NULL
        rg = rgobj
        actually.avgd.over = character(0) # keep track of this from emmeans calls
        for (i in seq_len(nrow(grps))) {
            sig = as.character(interaction(grps[i, ]))
            rows = which(grpids == sig)
            grd = rgobj@grid[rows, , drop = FALSE]
            lf = rgobj@linfct[rows, , drop = FALSE]
            # Reduce grid to infacs that actually appear in  this group
            nzg = grd[grd$.wgt. > 0, , drop = FALSE]
            rows = integer(0)
            # focus on levels of afac that exist in this group
            levs = unique(nzg[[afac]])
            rg@levels[[afac]] = levs
            rows = union(rows, which(grd[[afac]] %in% levs))
            rg@grid = grd[rows, , drop = FALSE]
            rg@linfct = lf[rows, , drop = FALSE]
            for (j in seq_along(grpfacs))
                rg@levels[[grpfacs[j]]] = grps[i, j]
            emmGrid = suppressMessages(emmeans(rg, gspecs, ...))
            actually.avgd.over = union(actually.avgd.over, emmGrid@misc$avgd.over)
            if (is.null(result))
                result = emmGrid
            else {
                result@grid = rbind(result@grid, emmGrid@grid)
                result@linfct = rbind(result@linfct, emmGrid@linfct)
            }
        }
        for (j in seq_along(grpfacs))
            result@levels[grpfacs[j]] = rgobj@levels[grpfacs[j]]
        
        result@misc$avgd.over = setdiff(actually.avgd.over, gspecs)
        result@misc$display = NULL
        nkeep = intersect(names(nesting), names(result@levels))
        if (length(nkeep) > 0)
            result@model.info$nesting = nesting[nkeep]
        else
            result@model.info$nesting = NULL
        
        # Note: if any nesting remains, this next call recurs back to this function
        result = emmeans(result, specs, by = by, ...)
    }
    
    if (length(xspecs) > 0)
        result@misc$display = .find.nonempty.nests(result, xspecs, nesting)

    # preserve any nesting that still exists
    nesting = nesting[names(nesting) %in% names(result@levels)]
    result@model.info$nesting =   if (length(nesting) > 0) nesting    else NULL
    
    # resolve 'by'
    by = orig.by
    if (length(xspecs <- intersect(by, names(nesting))))
        by = union(unlist(nesting[xspecs]), by)
    result@misc$by.vars = by
    
    result
}


### contrast function for nested structures
.nested_contrast = function(rgobj, method = "eff", by = NULL, adjust, ...) {
    nesting = rgobj@model.info$nesting
    # Prevent meaningless cases -- if A %in% B, we can't have A in 'by' without B
    # Our remedy will be to EXPAND the by list
    for (nm in intersect(by, names(nesting)))
        if (!all(nesting[[nm]] %in% by)) {
            by = union(by, nesting[[nm]])
            message("Note: Grouping factor(s) for '", nm, "' have been added to the 'by' list.")
        }

    if(!is.character(method))
        stop ("Non-character contrast methods are not supported with nested objects")
    
    testcon = get(paste0(method, ".emmc"))(1:3)
    if(missing(adjust)) 
        adjust = attr(testcon, "adjust")
    estType = attr(testcon, "type")

    wkrg = rgobj # working copy
    facs = setdiff(names(wkrg@levels), by)  # these are the factors we'll combine & contrast
    if (length(facs) == 0)
        stop("There are no factor levels left to contrast. Try taking nested factors out of 'by'.")
    if (!is.null(display <- wkrg@misc$display))
        wkrg = wkrg[which(display), drop.levels = TRUE]
    wkrg@model.info$nesting = wkrg@misc$display = NULL
    by.rows = .find.by.rows(wkrg@grid, by)
    if(length(by.rows) == 1)
        result = contrast(wkrg, method = method, by = by, ...)
    else {
        result = lapply(by.rows, function(rows) {
            contrast.emmGrid(wkrg[rows, drop.levels = TRUE], method = method, 
                              by = by, adjust = adjust, ...)
        })
        # Have to define .wgt. for nested emmGrid. Use average weight - seems most sensible
        for (i in seq_along(by.rows))
            result[[i]]@grid$.wgt. = mean(wkrg@grid[[".wgt."]][by.rows[[i]]])
        result$adjust = ifelse(is.null(adjust), "none", adjust)
        result = do.call(rbind.emmGrid, result)
        result = update(result, by = by, 
                        estType = ifelse(is.null(estType), "contrast", estType))
        cname = setdiff(names(result@levels), by)
        result@model.info$nesting[[cname]] = by
    }
    result@misc$orig.grid = result@misc$con.code = NULL

    for (nm in by) {
        if (nm %in% names(nesting))
            result@model.info$nesting[[nm]] = intersect(nesting[[nm]], by)
    }
    result
}


# Internal function to find nonempty cells in nested structures in rgobj for xfacs
# Returns logical vector, FALSE are rows of the grid we needn't display
.find.nonempty.nests = function(rgobj, xfacs = union(names(nesting), unlist(nesting)), 
                                nesting = rgobj@model.info$nesting) {
    grid = rgobj@grid
    keep = rep(TRUE, nrow(grid))
    for (x in xfacs) {
        facs = union(x, nesting[[x]])
        combs = do.call(expand.grid, rgobj@levels[facs])
        levs = as.character(interaction(combs))
        glevs = as.character(interaction(grid[facs]))
        
        for (lev in levs) {
            idx = which(glevs == lev)
            if (all(grid$.wgt.[idx] == 0)) {
                keep[idx] = FALSE
                levs[levs==lev] = ""
            }
        }
    }
    keep
}


# Internal function to find nesting
# We look at two things:
# (1) structural nesting - i.e., any combinations of
#     factors A and B for which each level of A occurs with one and only one
#     level of B. If so, we deem A %in% B.
# (2) Model-term nesting - cases where a factor appears not as a main effect
#     but only in higher-order terms. This is discovered using the 1s and 2s in 
#     trms$factors
# The function returns a named list, e.g., list(A = "B")
# If none found, an empty list is returned.
.find_nests = function(grid, trms, coerce, levels) {
    result = list()
    
    # only consider cases where levels has length > 1
    lng = sapply(levels, length)
    nms = names(levels[lng > 1])
    if (length(nms) < 2)
        return (result)
    g = grid[grid$.wgt. > 0, nms, drop = FALSE]
    for (nm in nms) {
        x = levels[[nm]]
        # exclude other factors this is already nested in
        excl = sapply(names(result), function(lnm)
            ifelse(nm %in% result[[lnm]], lnm, ""))
        otrs = setdiff(nms[!(nms == nm)], excl)
        max.levs = sapply(otrs, function(n) {
            max(sapply(x, function(lev) length(unique(g[[n]][g[[nm]] == lev]))))
        })
        if (any(max.levs == 1)) 
            result[[nm]] = otrs[max.levs == 1]
    }
    
    # Now look at factors attribute
    fac = attr(trms, "factors")
    if (length(fac) > 0) {
        if (!is.null(coerce)) for (stg in coerce) {
            subst = paste(.all.vars(stats::reformulate(stg)), collapse = ":")
            for (i in 1:2)
                dimnames(fac)[[i]] = gsub(stg, subst, dimnames(fac)[[i]], 
                                          fixed = TRUE)
        }
        fac = fac[intersect(nms, row.names(fac)), , drop = FALSE]
        
        ### new code
        nms = row.names(fac)
        cols = dimnames(fac)[[2]]
        pert = setdiff(nms, intersect(nms, cols)) # pertinent - no main effect in model
        for (nm in pert) {
            pfac = fac[, fac[nm, ] == 1, drop = FALSE]  # cols where nm appears
            if (ncol(pfac) == 0) { # case where there is no 1 in a row
                pfac = fac[, fac[nm, ] == 2, drop = FALSE]
                pfac[nm, ] = 1 # make own entry 1 so it isn't nested in self
            }
            nst = .strip.supersets(apply(pfac, 2, function(col) nms[col == 2]))
            if (length(nst) > 0)
                result[[nm]] = union(result[[nm]], nst)
        }
    }
    # include nesting factors that are themselves nested
    for (nm in names(result))
        result[[nm]] = union(unlist(result[result[[nm]]]), result[[nm]])
    
    result
}

# strip supersets from a list and condense down to a character vector
# e.g., lst = list("a", c("A", "B")) --> "A"
.strip.supersets = function(lst) {
    if (is.list(lst) && (length(lst) > 1)) {
        lst = lst[order(sapply(lst, length))]  # order by length
        for (i in length(lst):2) {
            tst = sapply(lst[1:(i-1)], function(x) all(x %in% lst[[i]]))
            if (any(tst)) lst[[i]] = NULL
        }
    }
    unique(unlist(lst))
}

# internal function to format a list of nested levels
.fmt.nest = function(nlist) {
    if (length(nlist) == 0)
        "none"
    else {
        tmp = lapply(nlist, function(x) 
            if (length(x) == 1) x
            else                paste0("(", paste(x, collapse = "*"), ")")
        )
        paste(sapply(names(nlist), function (nm) paste0(nm, " %in% ", tmp[[nm]])),
              collapse = ", ")
    }
}

# internal function to parse a nesting string & return a list
# spec can be a named list, character vector    ####, or formula
.parse_nest = function(spec) {
    if (is.null(spec))
        return(NULL)
    if (is.list(spec))
        return (spec)
    result = list()
    # break up any comma delimiters
    spec = trimws(unlist(strsplit(spec, ",")))
    for (s in spec) {
        parts = strsplit(s, "[ ]+%in%[ ]+")[[1]]
        grp = .all.vars(stats::reformulate(parts[2]))
        result[[parts[[1]]]] = grp
    }
    if(length(result) == 0)
        result = NULL
    result
}



# ### I'm removing this because I now think it creates more problems than it solves
# #
# # courtesy function to create levels for a nested structure factor %in% nest
# # factor: factor (or interaction() result)
# # ...:    factors in nest
# # SAS:    if (FALSE|TRUE), reference level in each nest is (first|last)
# nested = function(factor, ..., SAS = FALSE) {
#     nfacs = list(...)
#     if (length(nfacs) == 0)
#         return(factor)
#     nfacs$drop = TRUE
#     nest = do.call(interaction, nfacs)
#     result = as.character(interaction(factor, nest, sep = ".in."))
#     ores = unique(sort(result))
#     nlev = levels(nest)
#     flev = levels(factor)
#     refs = lapply(nlev, function(nst) {
#         r = ores[ores %in% paste0(flev, ".in.", nst)]
#         ifelse (SAS, rev(r)[1], r[1])
#     })
#     result[result %in% refs] = ".nref."
#     ores[ores %in% refs] = ".nref."
#     ores = setdiff(ores, ".nref.")
#     if (SAS)
#         factor(result, levels = c(ores, ".nref."))
#     else
#         factor(result, levels = c(".nref.", ores))
# }

