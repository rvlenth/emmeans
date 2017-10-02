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

### functions to implement different families of contrasts
### All return a matrix or data frame whose columns are the desired contrasts coefs
### with appropriate row and column names
### Also they have two attributes: 
###   "desc" is an expanded description of the family,
###   "adjust" is the default multiplicity adjustment (used if adjust="auto" in emmeans)

# all pairwise trt[i] - trt[j], i < j
#' Contrast families
#'
#' @rdname emmc-functions
#' @aliases emmc-functions 
#' 
#' @param levs Vector of factor levels
#' @param ... Additional arguments (these are ignored, but needed to make these functions 
#'   interchangeable)
#' Each contrast family has a default multiple-testing adjustment as noted
#' below. These adjustments are often only approximate; for a more exacting
#' adjustment, use the interfaces provided to \code{\link[multcomp]{glht}}
#' in the \pkg{multcomp} package.
#'
#' \code{pairwise.emmc}, \code{revpairwise.emmc}, and \code{tukey.emmc} generate
#' contrasts for all pairwise comparisons among least-squares means at the
#' levels in levs. The distinction is in which direction they are subtracted.
#' For factor levels A, B, C, D, \code{pairwise.emmc} generates the comparisons
#' A-B, A-C, A-D, B-C, B-D, and C-D, whereas \code{revpairwise.emmc} generates
#' B-A, C-A, C-B, D-A, D-B, and D-C. \code{tukey.emmc} invokes
#' \code{pairwise.emmc} or \code{revpairwise.emmc} depending on \code{reverse}. 
#' The default multiplicity adjustment method is \code{"tukey"}, which is only
#' approximate when the standard errors differ.
#'
#' \code{poly.emmc} generates orthogonal polynomial contrasts, assuming 
#' equally-spaced factor levels. These are derived from the
#' \code{\link[stats]{poly}} function, but an \emph{ad hoc} algorithm is used to
#' scale them to integer coefficients that are (usually) the same as in
#' published tables of orthogonal polynomial contrasts. The default multiplicity
#' adjustment method is \code{"none"}.
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
#'
#' \code{consec.emmc} and \code{mean_chg.emmc} are useful for contrasting
#' treatments that occur in sequence. For a factor with levels A, B, C, D, E,
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
#' @return A data.frame, each column containing contrast coefficients for levs.
#'   The "desc" attribute is used to label the results in emmeans, and the
#'   "adjust" attribute gives the default adjustment method for multiplicity.
#'
#' @examples
#' warp.lm <- lm(breaks ~ wool*tension, data = warpbreaks)
#' warp.emm <- emmeans(warp.lm, ~ tension | wool)
#' contrast(warp.emm, "poly")
#' 
#' ### Setting up a custom contrast function
#' helmert.emmc <- function(levs, ...) {
#'     M <- as.data.frame(contr.helmert(levs))
#'     names(M) <- paste(levs[-1],"vs earlier")
#'     attr(M, "desc") <- "Helmert contrasts"
#'     M
#' }
#' contrast(warp.emm, "helmert")
#' \dontrun{
#' # See what is used for polynomial contrasts with 6 levels
#' emmeans:::poly.emmc(1:6)
#' }
#' @name contrast-methods
pairwise.emmc = function(levs, ...) {
    k = length(levs)
    M = data.frame(levs=levs)
    for (i in seq_len(k-1)) {
        for (j in (i + seq_len(k-i))) { ###for (j in (i+1):k) {
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
    M
}

# all pairwise trt[j] - trt[i], j > i
#' @rdname emmc-functions
revpairwise.emmc = function(levs, ...) {
    k = length(levs)
    M = data.frame(levs=levs)
    for (i in 2:k) {
        for (j in seq_len(i-1)) {
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
    M
}

# pseudonym
#' @rdname emmc-functions
#' @param reverse Logical value to determine the direction of comparisons
tukey.emmc = function(levs, reverse = FALSE) {
    if (reverse)
        revpairwise.emmc(levs)
    else
        pairwise.emmc(levs)
}

# Poly contrasts - scaled w/ integer levels like most tables
# ad hoc scaling works for up to 13 levels
#' @rdname emmc-functions
#' @param max.degree Integer specifying the maximum degree of polynomial contrasts
poly.emmc = function(levs, max.degree = min(6, k-1)) {
    nm = c("linear", "quadratic", "cubic", "quartic", paste("degree",5:20))
    k = length(levs)
    M = as.data.frame(poly(seq_len(k), min(20,max.degree)))
    for (j in seq_len(ncol(M))) {
        con = M[ ,j]
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

# All comparisons with a control; ref = index of control group
# New version -- allows more than one control group (ref is a vector)
#' @rdname emmc-functions
#' @param ref Integer(s) specifying which level(s) to use as the reference
trt.vs.ctrl.emmc = function(levs, ref = 1) {
    if ((min(ref) < 1) || (max(ref) > length(levs)))
        stop("Reference levels are out of range")
    k = length(levs)
    cnm = ifelse(length(ref)==1, 
        levs[ref], 
        paste("avg(", paste(levs[ref], collapse=","), ")", sep=""))
    templ = rep(0, length(levs))
    templ[ref] = -1 / length(ref)
    M = data.frame(levs=levs)
    for (i in seq_len(k)) {
        if (i %in% ref) next
        con = templ
        con[i] = 1
        nm = paste(levs[i], cnm, sep = " - ")
        M[[nm]] = con
    }
    row.names(M) = levs
    M = M[-1]
    attr(M, "desc") = "differences from control"
    attr(M, "adjust") = "dunnettx"
    M
}

# control is 1st level
#' @rdname emmc-functions
trt.vs.ctrl1.emmc = function(levs, ...) {
    trt.vs.ctrl.emmc(levs, ref = 1)
}

# control is last level
#' @rdname emmc-functions
#' @inheritParams pairwise
trt.vs.ctrlk.emmc = function(levs, ...) {
    trt.vs.ctrl.emmc(levs, ref = length(levs))
}

# pseudonym
#' @rdname emmc-functions
dunnett.emmc = function(levs, ref = 1) {
    trt.vs.ctrl.emmc(levs, ref = ref)
}

# effects contrasts. Each mean versus the average of all
#' @rdname emmc-functions
eff.emmc = function(levs, ...) {
    k = length(levs)
    M = data.frame(levs=levs)
    for (i in seq_len(k)) {
        con = rep(-1/k, k)
        con[i] = (k-1)/k
        nm = paste(levs[i], "effect")
        M[[nm]] = con
    }
    row.names(M) = levs
    M = M[-1]
    attr(M, "desc") = "differences from grand mean"
    attr(M, "adjust") = "fdr"
    M
}

# "deleted" effects contrasts. 
# Each mean versus the average of all others
#' @rdname emmc-functions
del.eff.emmc = function(levs, ...) {
    k = length(levs)
    M = as.matrix(eff.emmc(levs,...)) * k / (k-1)
    M = as.data.frame(M)
    attr(M, "desc") = "differences from mean of others"
    attr(M, "adjust") = "fdr"
    M
}

# Contrasts to compare consecutive levels:
# (-1,1,0,0,...), (0,-1,1,0,...), ..., (0,...0,-1,1)
#' @rdname emmc-functions
consec.emmc = function(levs, reverse = FALSE, ...) {
    sgn = ifelse(reverse, -1, 1)
    k = length(levs)
    M = data.frame(levs=levs)
    for (i in seq_len(k-1)) {
        con = rep(0, k)
        con[i] = -sgn
        con[i+1] = sgn
        nm = ifelse(reverse,
                    paste(levs[i], "-", levs[i+1]),
                    paste(levs[i+1], "-", levs[i]))
        M[[nm]] = con
    }
    row.names(M) = levs
    M = M[-1]
    attr(M, "desc") = "changes between consecutive levels"
    attr(M, "adjust") = "mvt"
    M
}

# Mean after minus mean before
# e.g., (-1, 1/3,1/3,1/3), (-1/2,-1/2, 1/2,1/2), (-1/3,-1/3,-1/3, 1)
#' @rdname emmc-functions
mean_chg.emmc = function(levs, reverse = FALSE, ...) {
    sgn = ifelse(reverse, -1, 1)
    k = length(levs)
    M = data.frame(levs=levs)
    for (i in seq_len(k-1)) {
        kmi = k - i
        con = rep(c(-sgn/i, sgn/kmi), c(i, kmi)) 
        nm = paste(levs[i], levs[i+1], sep="|")
        M[[nm]] = con
    }
    row.names(M) = levs
    M = M[-1]
    attr(M, "desc") = "mean after minus mean before"
    attr(M, "adjust") = "mvt"
    M
}


