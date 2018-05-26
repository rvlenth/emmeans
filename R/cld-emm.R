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

# Runs the function multicompLetters from the multcompView package
# returns an error if not installed
.mcletters = function(..., Letters=c("1234567890",LETTERS,letters), reversed=FALSE) {
    if(!requireNamespace("multcompView", quietly = TRUE)) {
        stop("The 'multcompView' package must be installed to use cld methods")
        return (list(monospacedLetters = "?"))
    }
    
    # Expand strings to individual letters
    Letters = as.character(unlist(sapply(Letters, function(stg) {
        sapply(seq_len(nchar(stg)), function(i) substr(stg, i, i))
    })))
    
    result = multcompView::multcompLetters(..., Letters=Letters, reversed=reversed)
    if (is.null(result$monospacedLetters))
        result$monospacedLetters = result$Letters
    result
}

# Hack around providing a generic

#' @export cld
#' @rdname cld.emmGrid
cld = function (object, ...) {
    if(requireNamespace("multcomp", quietly = TRUE))
        multcomp::cld(object, ...)
    else
        UseMethod("cld")
}

# S3 method for emmGrid
#' Extract and display information on all pairwise comparisons of least-squares means.
#'
#' @aliases cld
#' 
#' @param object An object of class \code{emmGrid}
#' @param details Logical value determining whether detailed information on tests of 
#'   pairwise comparisons is displayed
#' @param sort Logical value determining whether the EMMs are sorted before the comparisons 
#'   are produced. When \code{TRUE}, the results are displayed according to 
#'   \code{reversed}.
#' @param by Character value giving the name or names of variables by which separate 
#'   families of comparisons are tested. If NULL, all means are compared. 
#'   If missing, the object's \code{by.vars} setting, if any, is used.
#' @param alpha Numeric value giving the significance level for the comparisons
#' @param Letters Character vector of letters to use in the display. Any strings of 
#'   length greater than 1 are expanded into individual characters
#' @param reversed Logical value (passed to \code{multcompView::multcompLetters}.) 
#'   If \code{TRUE}, the order of use of the letters is reversed. 
#'   In addition, if both \code{sort} and \code{reversed} are TRUE, the sort 
#'   order of results is reversed.
#' @param ... Arguments passed to \code{\link{contrast}} (for example, 
#'   an \code{adjust} method)
#'
#' This function uses the Piepho (2004) algorithm (as implemented in the 
#' \pkg{multcompView} package) to generate a compact letter display of all
#' pairwise comparisons of least-squares means. The function obtains (possibly
#' adjusted) P values for all pairwise comparisons of means, using the
#' \code{\link{contrast}} function with \code{method = "pairwise"}. When a P
#' value exceeds \code{alpha}, then the two means have at least one letter in
#' common.
#' 
#' @return When details == FALSE, an object of class \code{summary.ref_grid}
#'   (which inherits from \code{data.frame}) showing the summary of EMMs with 
#'   an added column named \code{.groups} containing the \code{cld} information. 
#'   When \code{details == TRUE}, a \code{list} with the object just described, 
#'   as well as the summary of the contrast results showing each comparison, 
#'   its estimate, standard error, t ratio, and adjusted P value.
#' 
#' @references Hans-Peter Piepho (2004) An algorithm for a letter-based representation 
#'   of all pairwise comparisons, Journal of Computational and Graphical Statistics, 13(2), 
#'   456-466.
#' 
#' @seealso \code{\link[multcomp]{cld}} in the \pkg{multcomp} package
#' 
#' @export cld.emmGrid
#'
#' @examples
#' warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
#' warp.emm <- emmeans(warp.lm, ~ tension | wool)
#' cld(warp.emm)                  # implicitly uses by = "wool"
#' cld(warp.emm, by = "tension")  # overrides implicit 'by'
#' 
#' # Mimic grouping bars and compare all 6 means
#' cld(warp.emm, by = NULL, Letters = "||||||||", alpha = .01)
cld.emmGrid = function(object, details=FALSE, sort=TRUE, 
                    by, alpha=.05, 
                    Letters = c("1234567890",LETTERS,letters), 
                    reversed=FALSE, ...) {
    emmtbl = summary(object, ...)
    if(missing(by)) 
        by = object@misc$by.vars
    if (sort) {
        args = list()
        for (nm in by) args[[nm]] = emmtbl[[nm]]
        args$.emmGrid. = emmtbl[[attr(emmtbl, "estName")]]
        ord = do.call("order", args)
        emmtbl = emmtbl[ord, , as.df = FALSE]
        object@grid = object@grid[ord, , drop = FALSE]
        object@linfct = object@linfct[ord, , drop = FALSE]
    }
    attr(emmtbl, "by.vars") = by
    object@misc$by.vars = by
    
    prwise = contrast(object, "revpairwise", by=by)    
    pwtbl = test(prwise, ...)
    
    p.boo = (pwtbl$p.value < alpha)
    if(is.null(by)) {
        by.rows = list(seq_len(nrow(pwtbl)))
        by.out = list(seq_len(nrow(emmtbl)))
    }
    else {
        by.rows = .find.by.rows(pwtbl, by)
        by.out = .find.by.rows(emmtbl, by)
    }
    # Create comps matrix reflecting order generated by pairwise.emmc
    icol = jcol = numeric(0)
    # create fake row indexes in revpairwise order for use by .mcletters
    k = length(by.out[[1]])
    for (i in 2:k) {
        icol = c(icol, seq_len(i-1))
        jcol = c(jcol, rep(i, i-1))
    }
    na.p = which(is.na(p.boo))
    # Take care of non-est cases. This is surprisingly complicated,
    # because it's possible we have some emmeans that are non-est
    # but comparisons are est'ble. So cases to exclude must be missing in
    # the table of means, AND appar somewhere in the indexes of NA p values
    # All that said, it still messes up because I didn't track the indexes correctly
    # excl.rows = intersect(which(is.na(emmtbl$SE)), union(icol[na.p], jcol[na.p]))
    # So I'll just go with which est's are missing
    excl.rows = which(is.na(emmtbl$SE))
    p.boo[na.p] = FALSE
    
    labs = paste(icol,jcol,sep="-")
    ltrs = rep("", nrow(emmtbl))
    for (i in seq_len(length(by.rows))) {
        pb = p.boo[by.rows[[i]]]
        names(pb) = labs
        mcl = .mcletters(pb, Letters = Letters, reversed = reversed)$monospacedLetters
        ltrs[by.out[[i]]] = paste(" ", mcl, sep="")
    }
    # any missing estimates get blanks...
    ltrs[excl.rows] = ""
    
    emmtbl[[".group"]] = ltrs
    if(sort && reversed) for (i in seq_len(length(by.out))) {
        r = by.out[[i]]
        emmtbl[r, ] = emmtbl[rev(r), ]
    }
    
    attr(emmtbl, "mesg") = c(attr(emmtbl,"mesg"), attr(pwtbl, "mesg"), 
                             paste("significance level used: alpha =", alpha))
        
    if (details)
        list(emmeans = emmtbl, comparisons = pwtbl)
    else
        emmtbl
}
