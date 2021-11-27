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

# Runs the function multcompLetters from the multcompView package
# returns an error if not installed
.mcletters = function(..., Letters=c("1234567890",LETTERS,letters), reversed=FALSE) {
    .requireNS("multcompView", 
               "The 'multcompView' package must be installed to use CLD methods")

    # Expand strings to individual letters
    Letters = as.character(unlist(sapply(Letters, function(stg) {
        sapply(seq_len(nchar(stg)), function(i) substr(stg, i, i))
    })))
    
    result = multcompView::multcompLetters(..., Letters=Letters, reversed=reversed)
    if (is.null(result$monospacedLetters))
        result$monospacedLetters = result$Letters
    result
}


### Lingering support for multcomp::cld -- registered dynamically in zzz.R
### NOTE: MUST KEEP the rdname of CLD.emmGrid 
###       because it's referenced by augmentedRCBD package
#' Compact letter displays
#' 
#' A method for \code{multcomp::cld()} is provided for users desiring to produce 
#' compact-letter displays (CLDs). 
#' This method uses the Piepho (2004) algorithm (as implemented in the
#' \pkg{multcompView} package) to generate a compact letter display of all
#' pairwise comparisons of estimated marginal means. The function obtains (possibly
#' adjusted) P values for all pairwise comparisons of means, using the
#' \code{\link{contrast}} function with \code{method = "pairwise"}. When a P
#' value exceeds \code{alpha}, then the two means have at least one letter in
#' common.
#' 
#' @rdname CLD.emmGrid
#' @order 1
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
#' @references Piepho, Hans-Peter (2004) An algorithm for a letter-based 
#'   representation of all pairwise comparisons, 
#'   Journal of Computational and Graphical Statistics, 
#'   13(2), 456-466.
#'   
#' @note
#' We warn that such displays encourage a poor
#' practice in interpreting significance tests. CLDs are misleading because they
#' visually group means with comparisons \emph{P} > \code{alpha} as though they 
#' are equal, when in fact we have only failed to prove that they differ.
#' As alternatives, consider \code{\link{pwpp}} (graphical display of \emph{P} 
#' values) or \code{\link{pwpm}} (matrix display).
#'
#' @method cld emmGrid
#' @examples
#' if(requireNamespace("multcomp")) {
#'     pigs.lm <- lm(log(conc) ~ source + factor(percent), data = pigs)
#'     pigs.emm <- emmeans(pigs.lm, "percent", type = "response")
#'     multcomp::cld(pigs.emm, alpha = 0.10, Letters = LETTERS)
#' }
cld.emmGrid = function(object, details=FALSE, sort=TRUE, 
                       by, alpha=.05, 
                       Letters = c("1234567890",LETTERS,letters), 
                       reversed=FALSE, ...) {
    if (!is.na(object@post.beta)[1]) {
        message("NOTE: Summary and groupings are based on frequentist results")
        object@post.beta = matrix(NA)
    }
    emmtbl = summary(object, ...)
    if(missing(by)) 
        by = object@misc$by.vars
    if (sort) {
        args = list()
        for (nm in by) args[[nm]] = emmtbl[[nm]]
        args$.emmGrid. = emmtbl[[attr(emmtbl, "estName")]]
        ord = do.call("order", unname(args))
        emmtbl = emmtbl[ord, , as.df = FALSE]
        if (!is.null(object@misc$display)) {
            use = which(object@misc$display)
            object@linfct = object@linfct[use, , drop = FALSE]
            object@grid = object@grid[use, , drop = FALSE]
            object@misc$display = NULL
        }
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
    
    ### This code moved to inside the loop
    # Create comps matrix reflecting order generated by pairwise.emmc
    # icol = jcol = numeric(0)
    # create fake row indexes in revpairwise order for use by .mcletters
    # k = length(by.out[[1]])
    # for (i in 2:k) {
    #     icol = c(icol, seq_len(i-1))
    #     jcol = c(jcol, rep(i, i-1))
    # }
    
    ltrs = rep("", nrow(emmtbl))
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
    
    for (i in seq_len(length(by.rows))) {
        # Create comps matrix reflecting order generated by pairwise.emmc
        icol = jcol = numeric(0)
        k = length(by.out[[i]])
        for (j in 2:k) {
            icol = c(icol, seq_len(j-1))
            jcol = c(jcol, rep(j, j-1))
        }
        labs = paste(icol, jcol, sep="-")
        pb = p.boo[by.rows[[i]]]
        names(pb) = labs
        mcl = .mcletters(pb, Letters = Letters, reversed = reversed)$monospacedLetters
        ltrs[by.out[[i]]] = paste0(" ", mcl[seq_along(by.out[[i]])])
    }
    # any missing estimates get blanks...
    ltrs[excl.rows] = ""
    
    emmtbl[[".group"]] = ltrs
    if(sort && reversed) for (i in seq_len(length(by.out))) {
        r = by.out[[i]]
        emmtbl[r, ] = emmtbl[rev(r), ]
    }
    
    dontusemsg = paste0("NOTE: Compact letter displays can be misleading\n",
                        "      because they show NON-findings rather than findings.\n",
                        "      Consider using 'pairs()', 'pwpp()', or 'pwpm()' instead.")
    
    attr(emmtbl, "mesg") = c(attr(emmtbl,"mesg"), attr(pwtbl, "mesg"), 
                             paste("significance level used: alpha =", alpha), dontusemsg)
    
    if (details)
        list(emmeans = emmtbl, comparisons = pwtbl)
    else
        emmtbl
}

### Registered dynamically in zzz.R
### NOTE: MUST KEEP the rdname of CLD.emmGrid 
###       because it's referenced by augmentedRCBD package
#' @rdname CLD.emmGrid
#' @order 2
#' @method cld emm_list
#' @param which Which element of the \code{emm_list} object to process
#'   (If length exceeds one, only the first one is used)
cld.emm_list = function(object, ..., which = 1) {
    multcomp::cld(object[[which[1]]], ...)
}

