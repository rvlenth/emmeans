##############################################################################
#    Copyright (c) 2012-2019 Russell V. Lenth                                #
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

### functions to implement different families of contrasts
### All return a matrix or data frame whose columns are the desired contrasts coefs
### with appropriate row and column names
### Also they have two attributes: 
###   "desc" is an expanded description of the family,
###   "adjust" is the default multiplicity adjustment (used if adjust="auto" in emmeans)

#' Contrast families
#'
#' Functions with an extension of \code{.emmc} provide for named contrast
#' families. One of the standard ones documented here may be used, or the user
#' may write such a function.
#'
#' Each standard contrast family has a default multiple-testing adjustment as
#' noted below. These adjustments are often only approximate; for a more
#' exacting adjustment, use the interfaces provided to \code{glht} in the
#' \pkg{multcomp} package.
#'
#' \code{pairwise.emmc}, \code{revpairwise.emmc}, and \code{tukey.emmc} generate
#' contrasts for all pairwise comparisons among estimated marginal means at the
#' levels in levs. The distinction is in which direction they are subtracted.
#' For factor levels A, B, C, D, \code{pairwise.emmc} generates the comparisons
#' A-B, A-C, A-D, B-C, B-D, and C-D, whereas \code{revpairwise.emmc} generates
#' B-A, C-A, C-B, D-A, D-B, and D-C. \code{tukey.emmc} invokes
#' \code{pairwise.emmc} or \code{revpairwise.emmc} depending on \code{reverse}.
#' The default multiplicity adjustment method is \code{"tukey"}, which is only
#' approximate when the standard errors differ.
#'
#' \code{poly.emmc} and \code{opoly.emmc} generate orthogonal polynomial contrasts.
#' \code{poly.emmc} uses equally-spaced factor levels; coefficients are derived from the
#' \code{\link[stats]{poly}} function, but an \emph{ad hoc} algorithm is used to
#' scale them to integer coefficients that are (usually) the same as in
#' published tables of orthogonal polynomial contrasts. On the other hand, 
#' \code{opoly.emmc}'s coefficients are always normalized (sum of squares equals 1),
#' but allows the user to choose alternate reference points in \code{scores}, as in the
#' \code{\link{contr.poly}} function.
#' In both cases, the default multiplicity adjustment method is \code{"none"}.
#'
#' \code{trt.vs.ctrl.emmc} and its relatives generate contrasts for comparing
#' one level (or the average over specified levels) with each of the other
#' levels. The argument \code{ref} should be the index(es) (not the labels) of
#' the reference level(s). \code{trt.vs.ctrl1.emmc} is the same as
#' \code{trt.vs.ctrl.emmc} with a reference value of 1, and
#' \code{trt.vs.ctrlk.emmc} is the same as \code{trt.vs.ctrl} with a reference
#' value of \code{length(levs)}. \code{dunnett.emmc} is the same as
#' \code{trt.vs.ctrl}. The default multiplicity adjustment method is
#' \code{"dunnettx"}, a close approximation to the Dunnett adjustment.
#' \emph{Note} in all of these functions, it is illegal to have any overlap
#' between the \code{ref} levels and the \code{exclude} levels. If any is found,
#' an error is thrown.
#'
#' \code{consec.emmc} and \code{mean_chg.emmc} are useful for contrasting
#' treatments that occur in sequence. For a factor with levels A, B, C, D,
#' \code{consec.emmc} generates the comparisons B-A, C-B, and D-C, while
#' \code{mean_chg.emmc} generates the contrasts (B+C+D)/3 - A, (C+D)/2 -
#' (A+B)/2, and D - (A+B+C)/3. With \code{reverse = TRUE}, these differences go
#' in the opposite direction.
#'
#' \code{eff.emmc} and \code{del.eff.emmc} generate contrasts that compare each
#' level with the average over all levels (in \code{eff.emmc}) or over all other
#' levels (in \code{del.eff.emmc}). These differ only in how they are scaled.
#' For a set of k EMMs, \code{del.eff.emmc} gives weight 1 to one EMM and weight
#' -1/(k-1) to the others, while \code{eff.emmc} gives weights (k-1)/k and -1/k
#' respectively, as in subtracting the overall EMM from each EMM. The default
#' multiplicity adjustment method is \code{"fdr"}. This is a Bonferroni-based
#' method and is slightly conservative; see \code{\link[stats]{p.adjust}}.
#' 
#' \code{nrmlz.emmc} is a wrapper that can be used with any other \code{.emmc} 
#' function that will normalize the contrast coefficients so that the sum of its 
#' squares equals 1. Just provide the root name of the function in \code{family},
#' along with any other arguments to pass to it.
#' 
#' \code{wtcon.emmc} generates weighted contrasts based on the function 
#' \code{\link[multcomp]{contrMat}} function in the \pkg{multcomp} package,
#' using the provided \code{type} as documented there. If the user provides 
#' \code{wts}, they have to conform to the length of \code{levs}; however,
#' if \code{wts} is not specified, \code{contrast} will fill-in what is
#' required, and usually this is safer (especially when \code{by != NULL}
#' which usually means that the weights are different in each \code{by} group).
#'
#' \code{identity.emmc} simply returns the identity matrix (as a data frame),
#' minus any columns specified in \code{exclude}. It is potentially useful in
#' cases where a contrast function must be specified, but none is desired.
#'
#' @rdname emmc-functions
#' @aliases emmc-functions
#' @param levs Vector of factor levels
#' @param exclude integer vector of indices, or character vector of levels to
#'   exclude from consideration. These levels will receive weight 0 in all
#'   contrasts. Character levels must exactly match elements of \code{levs}.
#' @param include integer or character vector of levels to include (the
#'   complement of \code{exclude}). An error will result if the user specifies
#'   both \code{exclude} and \code{include}.
#' @param ... Additional arguments, passed to related methods as appropriate
#'
#' @return A data.frame, each column containing contrast coefficients for levs.
#'   The "desc" attribute is used to label the results in emmeans, and the
#'   "adjust" attribute gives the default adjustment method for multiplicity.
#'
#' @note Caution is needed in cases where the user alters the ordering of
#'   results (e.g., using the the \code{"[...]"} operator), because the
#'   contrasts generated depend on the order of the levels provided. For
#'   example, suppose \code{trt.vs.ctrl1} contrasts are applied to two \code{by}
#'   groups with levels ordered (Ctrl, T1, T2) and (T1, T2, Ctrl) respectively,
#'   then the contrasts generated will be for (T1 - Ctrl, T2 - Ctrl) in the
#'   first group and (T2 - T1, Ctrl - T1) in the second group, because the first
#'   level in each group is used as the reference level.
#'
#' @examples
#' warp.lm <- lm(breaks ~ wool*tension, data = warpbreaks)
#' warp.emm <- emmeans(warp.lm, ~ tension | wool)
#' contrast(warp.emm, "poly")
#' contrast(warp.emm, "trt.vs.ctrl", ref = "M")
#' \dontrun{
#' ## Same when enhanced labeling is used:
#' contrast(warp.emm, "trt.vs.ctrl", 
#'          enhance.levels = "tension", ref = "tensionM")}
#' 
#' # Comparisons with grand mean
#' contrast(warp.emm, "eff")
#' # Comparisons with a weighted grand mean
#' contrast(warp.emm, "eff", wts = c(2, 5, 3))
#'
#' # Compare only low and high tensions
#' # Note pairs(emm, ...) calls contrast(emm, "pairwise", ...)
#' pairs(warp.emm, exclude = 2)
#' # (same results using exclude = "M" or include = c("L","H") or include = c(1,3))
#'
#'# Same contrasts as above but with normalized contrast coefficients
#'contrast(warp.emm, "nrmlz", family = "pairwise", include = c(1, 3))
#'
#' ### Setting up a custom contrast function
#' revhelmert.emmc <- function(levs, ...) {
#'     M <- as.data.frame(contr.helmert(levs)[rev(seq_along(levs)), ])
#'     names(M) <- paste(rev(levs)[-1],"vs later")
#'     attr(M, "desc") <- "reverse Helmert contrasts"
#'     M
#' }
#' contrast(warp.emm, "revhelmert")
#' 
#' \dontrun{
#' # See what is used for polynomial contrasts with 6 levels
#' emmeans:::poly.emmc(1:6)
#' }
#' @name contrast-methods
#' @export
pairwise.emmc = function(levs, exclude = integer(0), include, ...) {
    exclude = .get.excl(levs, exclude, include)
    k = length(levs)
    M = data.frame(levs=levs)
    for (i in setdiff(seq_len(k-1), exclude)) {
        for (j in setdiff(i + seq_len(k-i), exclude)) { ###for (j in (i+1):k) {
            con = rep(0,k)
            con[i] = 1
            con[j] = -1
            nm = paste(levs[i], levs[j], sep = " - ")
            M[[nm]] = con
        }
    }
    row.names(M) = levs
    M = M[-1]
    attr(M, "desc") = "pairwise differences"
    attr(M, "adjust") = "tukey"
    attr(M, "type") = "pairs"
    attr(M, "famSize") = k - length(exclude)
    if(length(exclude) > 0)
        attr(M, "famSize") = length(levs) - length(exclude)
    M
}

# all pairwise trt[j] - trt[i], j > i
#' @rdname emmc-functions
#' @export
revpairwise.emmc = function(levs, exclude = integer(0), include, ...) {
    exclude = .get.excl(levs, exclude, include)
    k = length(levs)
    M = data.frame(levs=levs)
    for (i in setdiff(1 + seq_len(k - 1), exclude)) {
        for (j in setdiff(seq_len(i-1), exclude)) {
            con = rep(0,k)
            con[i] = 1
            con[j] = -1
            nm = paste(levs[i], levs[j], sep = " - ")
            M[[nm]] = con
        }
    }
    row.names(M) = levs
    M = M[-1]
    attr(M, "desc") = "pairwise differences"
    attr(M, "adjust") = "tukey"
    attr(M, "type") = "pairs"
    if(length(exclude) > 0)
        attr(M, "famSize") = length(levs) - length(exclude)
    M
}

# pseudonym
#' @rdname emmc-functions
#' @param reverse Logical value to determine the direction of comparisons
#' @export
tukey.emmc = function(levs, reverse = FALSE, ...) {
    if (reverse)
        revpairwise.emmc(levs, ...)
    else
        pairwise.emmc(levs, ...)
}

# Poly contrasts - scaled w/ integer levels like most tables
# ad hoc scaling works for up to 13 levels
#' @rdname emmc-functions
#' @param max.degree Integer specifying the maximum degree of polynomial contrasts
#' @export
poly.emmc = function(levs, max.degree = min(6, k-1), ...) {
    nm = c("linear", "quadratic", "cubic", "quartic", paste("degree",5:20))
    k = length(levs)
    M = as.data.frame(poly(seq_len(k), min(20,max.degree)))
    for (j in seq_len(ncol(M))) {
        con = M[, j]
        pos = which(con > .01)
        con = con / min(con[pos])
        z = max(abs(con - round(con)))
        while (z > .05) {
            con = con / z
            z = max(abs(con - round(con)))
        }
        M[ ,j] = round(con)
    }
    row.names(M) = levs
    names(M) = nm[seq_len(ncol(M))]
    attr(M, "desc") = "polynomial contrasts"
    attr(M, "adjust") = "none"
    M
}

# Orthonormal poly contrasts - not rescaled
#' @rdname emmc-functions
#' @param scores Set of values of length \code{length(levs)} over which
#'   orthogonal polynomials are computed. The default scores are the 
#'   consecutive integers \code{seq_along(levs)}.
#'   (If \code{exclude} or \code{include}
#'   are used, the default scores are subsetted accordingly; and if \code{scores} 
#'   is specified, its length must be the same as that of the subsetted \code{levs}).
#' @export
opoly.emmc = function(levs, max.degree = min(6, k-1), scores, 
                      exclude = integer(0), include, ...) {
    if (MS <- missing(scores))
        scores = seq_along(levs)
    exclude = .get.excl(levs, exclude, include)
    if(length(exclude) > 0) {
        oldlevs = levs
        if (MS)
            scores = scores[-exclude]
        levs = levs[-exclude]
    }
    nm = c("linear", "quadratic", "cubic", "quartic", paste("degree",5:20))
    k = length(levs)
    if (length(scores) != k)
        stop("In opoly.emmc: Lengths of 'scores' must equal ", length(levs), call. = FALSE)
    M = contr.poly(k, contrasts = TRUE, scores = scores, sparse = FALSE)
    if(length(exclude) > 0) {
        MM = M
        levs = oldlevs
        M = matrix(0, nrow = length(oldlevs), ncol = ncol(MM))
        M[-exclude, ] = MM
    }
    M = as.data.frame(M[, seq_len(min(k-1, max.degree)), drop = FALSE])
    row.names(M) = levs
    names(M) = nm[seq_len(ncol(M))]
    attr(M, "desc") = "normalized polynomial contrasts"
    attr(M, "adjust") = "none"
    M
}


# All comparisons with a control; ref = index of control group
# New version -- allows more than one control group (ref is a vector)
#' @rdname emmc-functions
#' @param ref Integer(s) or character(s) specifying which level(s) to use 
#'   as the reference. Character values must exactly match elements of \code{levs}
#'   (including any enhancements -- see examples)
#' @export
trt.vs.ctrl.emmc = function(levs, ref = 1, reverse = FALSE, 
                            exclude = integer(0), include, ...) {
    ref = .num.key(levs, ref)
    exclude = .get.excl(levs, exclude, include)
    if (length(ref) == 0 || (min(ref) < 1) || (max(ref) > length(levs)))
        stop("In trt.vs.ctrl.emmc(), 'ref' levels are out of range", call. = FALSE)
    k = length(levs)
    cnm = ifelse(length(ref)==1, 
                 levs[ref], 
                 paste("avg(", paste(levs[ref], collapse=","), ")", sep=""))
    templ = rep(0, length(levs))
    templ[ref] = -1 / length(ref)
    M = data.frame(levs=levs)
    if (length(intersect(exclude, ref)) > 0)
        stop("'exclude' set cannot overlap with 'ref'")
    skip = c(ref, exclude)
    for (i in seq_len(k)) {
        if (i %in% skip) next
        con = templ
        con[i] = 1
        if (reverse)
            nm = paste(cnm, levs[i], sep = " - ")
        else
            nm = paste(levs[i], cnm, sep = " - ")
        M[[nm]] = con
    }
    row.names(M) = levs
    M = M[-1]
    if (reverse)
        M = -M
    attr(M, "desc") = "differences from control"
    attr(M, "adjust") = "dunnettx"
    if(length(exclude) > 0)
        attr(M, "famSize") = length(levs) - length(exclude)
    M
}

# control is 1st level
#' @rdname emmc-functions
#' @export
trt.vs.ctrl1.emmc = function(levs, ref = 1, ...) {
    trt.vs.ctrl.emmc(levs, ref = ref, ...)
}

# control is last level
#' @rdname emmc-functions
#' @export
trt.vs.ctrlk.emmc = function(levs, ref = length(levs), ...) {
    trt.vs.ctrl.emmc(levs, ref = ref, ...)
}

# pseudonym for trt.vs.ctrl
#' @rdname emmc-functions
#' @export
dunnett.emmc = function(levs, ref = 1, ...) {
    trt.vs.ctrl.emmc(levs, ref = ref, ...)
}

# effects contrasts. Each mean versus the average of all

#' @rdname emmc-functions
#' @param wts Optional weights to use with \code{eff.emmc} and \code{del.eff.emmc} contrasts.
#'   These default to equal weights.
#'   If \code{exclude} or \code{include} are specified, \code{wts} may be
#'   either the same length as \code{levs} or the length of the included levels.
#'   In the former case, weights for any excluded levels are set to zero.
#'   \code{wts} has no impact on the results unless there are at least
#'   three levels included in the contrast.
#' @export
eff.emmc = function(levs, exclude = integer(0), include, wts = rep(1, length(levs)), ...) {
    exclude = .get.excl(levs, exclude, include)
    if ((length(exclude) > 0) && (length(wts) == length(levs) - length(exclude))) {
        tmp = rep(0, length(levs))
        tmp[-exclude] = wts
        wts = tmp
    }
    if (length(wts) != length(levs))
        stop("length of 'wts' must equal the number of levels",
             " or the number of included levels")
    wts[exclude] = 0
    
    M = data.frame(row.names = levs)
    wts = wts / sum(wts)
    for (i in setdiff(seq_along(levs), exclude)) {
        con = -wts
        con[i] = 1 + con[i]
        nm = paste(levs[i], "effect")
        M[[nm]] = con
    }
    attr(M, "desc") = "differences from grand mean"
    attr(M, "adjust") = "fdr"
    if(length(exclude) > 0)
        attr(M, "famSize") = length(levs) - length(exclude)
    M
}

#' @rdname emmc-functions
#' @export
del.eff.emmc = function(levs, exclude = integer(0), include, wts = rep(1, length(levs)), ...) {
    M = eff.emmc(levs, exclude, include, wts, ...)
    use = setdiff(seq_along(levs), .get.excl(levs, exclude, include))
    for(i in seq_along(M))
        M[[i]] = M[[i]] / M[[i]][use[i]]
    attr(M, "desc") = "differences from mean of others"
    M
}

# Contrasts to compare consecutive levels:
# (-1,1,0,0,...), (0,-1,1,0,...), ..., (0,...0,-1,1)
#' @rdname emmc-functions
#' @export
consec.emmc = function(levs, reverse = FALSE, exclude = integer(0), include, ...) {
    exclude = .get.excl(levs, exclude, include)
    sgn = ifelse(reverse, -1, 1)
    tmp = rep(0, length(levs))
    k = length(levs) - length(exclude)
    active.rows = setdiff(seq_along(levs), exclude)
    M = data.frame(levs=levs)
    nms = levs[active.rows]
    for (i in seq_len(k-1)) {
        con = rep(0, k)
        con[i] = -sgn
        con[i+1] = sgn
        tmp[active.rows] = con
        nm = ifelse(reverse,
                    paste(nms[i], "-", nms[i+1]),
                    paste(nms[i+1], "-", nms[i]))
        M[[nm]] = tmp
    }
    row.names(M) = levs
    M = M[-1]
    attr(M, "desc") = "changes between consecutive levels"
    attr(M, "adjust") = "mvt"
    if(length(exclude) > 0)
        attr(M, "famSize") = length(levs) - length(exclude)
    M
}

# Mean after minus mean before
# e.g., (-1, 1/3,1/3,1/3), (-1/2,-1/2, 1/2,1/2), (-1/3,-1/3,-1/3, 1)
#' @rdname emmc-functions
#' @export
mean_chg.emmc = function(levs, reverse = FALSE, exclude = integer(0), include, ...) {
    exclude = .get.excl(levs, exclude, include)
    sgn = ifelse(reverse, -1, 1)
    k = length(levs) - length(exclude)
    tmp = rep(0, length(levs))
    M = data.frame(levs=levs)
    active.rows = setdiff(seq_along(levs), exclude)
    nms = levs[active.rows]
    for (i in seq_len(k-1)) {
        kmi = k - i
        con = rep(c(-sgn/i, sgn/kmi), c(i, kmi)) 
        nm = paste(nms[i], nms[i+1], sep="|")
        tmp[active.rows] = con
        M[[nm]] = tmp
    }
    row.names(M) = levs
    M = M[-1]
    attr(M, "desc") = "mean after minus mean before"
    attr(M, "adjust") = "mvt"
    if(length(exclude) > 0)
        attr(M, "famSize") = length(levs) - length(exclude)
    M
}

#' @rdname emmc-functions
#' @export
helmert.emmc <- function(levs, exclude = integer(0), include, ...) {
    exclude = .get.excl(levs, exclude, include)
    k = length(levs) - length(exclude)
    M = contr.helmert(k)
    if (length(exclude) > 0) {
        M = matrix(0, nrow = length(levs), ncol = k-1)
        M[-exclude, ] = contr.helmert(k)
        lbl = levs[-exclude]
    }
    else
        lbl = levs
    M <- as.data.frame(M)
    names(M) <- paste(lbl[-1],"vs earlier")
    attr(M, "desc") <- "Helmert contrasts"
    attr(M, "adjust") = "none"
    if(length(exclude) > 0)
        attr(M, "famSize") = k
    M
}

#' @rdname emmc-functions
#' @param family name of contrast family to use
#' @export
nrmlz.emmc = function(levs, family, ...) {
    M = get(paste0(family, ".emmc"))(levs, ...)
    for (i in seq_len(ncol(M)))
        M[, i] = M[, i] / sqrt(sum(M[, i]^2))
    attr(M, "desc") = paste("Normalized", attr(M, "desc"))
    M
}

# weighted contrasts
#' @rdname emmc-functions
#' @param cmtype the \code{type} argument passed to \code{\link[multcomp]{contrMat}}
#' @export
wtcon.emmc = function(levs, wts = rep(1, length(levs)), cmtype = "GrandMean", ...) {
    if (!requireNamespace("multcomp"))
        stop("The 'multcomp' package must be installed to use 'wtcon' contrasts", call. = FALSE)
    
    names(wts) = levs
    M = data.frame(t(multcomp::contrMat(wts, type = cmtype, ...)))
    M
}



# Non-contrasts -- just pass thru, possibly excluding some levels
#' @rdname emmc-functions
#' @export
identity.emmc = function(levs, exclude = integer(0), include, ...) {
    exclude = .get.excl(levs, exclude, include)
    k = length(levs) - length(exclude)
    M = as.data.frame(diag(length(levs)))
    names(M) = levs
    if(length(exclude) > 0)
        M = M[ , -exclude, drop = FALSE]
    attr(M, "desc") = "Identity"
    attr(M, "famSize") = k
    attr(M, "adjust") = "none"
    M
}



### utility to translate character keys to index keys
#' @rdname extending-emmeans
#' @order 43
#' @param levs,key The \code{.num.key} function returns the numeric indices of
#' the levels in \code{levs} to the set of all levels in \code{key}
#' @export
.num.key = function(levs, key) {
    if(!is.null(raw <- attr(levs, "raw")))
        levs = raw
    orig.key = key
    if (is.character(key))
        key = match(key, levs)
    # if (any(is.na(key)))
    #     warning("One or more of: '", paste(orig.key, collapse = "','"), "' not found in '",
    #          paste(levs, collapse = "','"), "'", call. = FALSE)
    # if (any(key > length(levs)) || any(key < 1))
    #     stop("Numeric index not in 1 : length(levs)")
    key[key %in% seq_along(levs)] # I think I'll just silently remove unmatched levels
}

### utility to find exclude levels from either exclude or include
### Also returns numeric version
#' @rdname extending-emmeans
#' @order 34
#' @param exc,inc Arguments for \code{.get.excl} which is useful
#' in writing \code{.emmc} functions for generating contrast coefficients,
#' and supports arguments \code{exclude} or \code{include} for excluding
#' or specifying which levels to use.
#' @export
.get.excl = function(levs, exc, inc) {
    if (!missing(inc)) {
        if(length(exc) > 0)
            stop("Cannot specify both 'exclude' and 'include'", call. = FALSE)
        inc = .num.key(levs, inc)
        exc = setdiff(seq_along(levs), inc)
    }
    .num.key(levs, exc)
}