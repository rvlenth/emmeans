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

# confint and test methods

# confint method
#' @rdname summary.emmGrid
#' @param parm (Required argument for \code{confint} methods, but not used)
#' @method confint emmGrid
#' @export
confint.emmGrid = function(object, parm, level = .95, ...) {
    summary(object, infer = c(TRUE, FALSE), level = level, ...)
}

#' @rdname summary.emmGrid
#' @export
test = function(object, null, ...) {
    UseMethod("test")
}

#' @rdname summary.emmGrid
#' @param joint Logical value. If \code{FALSE}, the arguments are passed to 
#'   \code{\link{summary.emmGrid}} with \code{infer=c(FALSE, TRUE)}. If \code{joint = 
#'   TRUE}, a joint test of the hypothesis L beta = null is performed, where L 
#'   is \code{object@linfct} and beta is the vector of fixed effects estimated 
#'   by \code{object@betahat}. This will be either an \emph{F} test or a 
#'   chi-square (Wald) test depending on whether degrees of freedom are 
#'   available.
#' @param verbose Logical value. If \code{TRUE} and \code{joint = TRUE}, a table
#'   of the effects being tested is printed.
#' @param rows Integer values. The rows of L to be tested in the joint test. If
#'   missing, all rows of L are used. If not missing, \code{by} variables are
#'   ignored.
#' @param status logical. If \code{TRUE}, a \code{note} column showing status
#'   flags (for rank deficiencies and estimability issues) is displayed even 
#'   when empty. If \code{FALSE}, the column is included only if there are 
#'   such issues.
#' @method test emmGrid
#' @export
test.emmGrid = function(object, null = 0, 
                    joint = FALSE, verbose = FALSE, rows, by, status = FALSE, ...) {
    # if joint = FALSE, this is a courtesy method for 'contrast'
    # else it computes the F test or Wald test of H0: L*beta = null
    # where L = object@linfct    
    if (!joint) {
        if (missing(by))
            summary(object, infer=c(FALSE,TRUE), null = null, ...)
        else
            summary(object, infer=c(FALSE,TRUE), null = null, by = by, ...)
    }
    else {
        if(verbose) {
            cat("Joint test of the following linear predictions\n")
            print(cbind(object@grid, equals = null))
        } 
        L = object@linfct
        bhat = object@bhat
        estble.idx = which(!is.na(object@bhat))
        bhat = bhat[estble.idx]
        est.flag = !is.na(object@nbasis[1])

        if (!missing(rows))
            by.rows = list(sel.rows = rows)
        else {
            by.rows = list(all = seq_len(nrow(L)))
            if(missing(by)) 
                by = object@misc$by.vars
            if (!is.null(by)) 
                by.rows = .find.by.rows(object@grid, by)
        }
        
        result = lapply(by.rows, function(rows) {
            LL = L[rows, , drop = FALSE]

            # discard any rows that have NAs
            narows = apply(LL, 1, function(x) any(is.na(x)))
            LL = LL[!narows, , drop = FALSE]
            rrflag = 0 + 2 * any(narows)  ## flag for estimability issue
            
            if(est.flag)  { 
                
                if (any(!estimability::is.estble(LL, object@nbasis)))
                    rrflag = bitwOr(rrflag, 2)
                LL = LL[, estble.idx, drop = FALSE]
            }
            # Check rank
            qrLt = qr(t(LL)) # this will work even if LL has 0 rows
            r = qrLt$rank
            if (r == 0)
                return(c(df1 = 0, df2 = NA, F.ratio = NA, p.value = NA, note = 3))
            
            if (r < nrow(LL)) {
                if(!all(null == 0))
                    stop("Rows are linearly dependent - cannot do the test when 'null' != 0")
                rrflag = bitwOr(rrflag, 1)
            }
            tR = t(qr.R(qrLt))[1:r, 1:r, drop = FALSE]
            tQ = t(qr.Q(qrLt))[1:r, , drop = FALSE]
            if(length(null) < r) null = rep(null, r)
            z = tQ %*% bhat - solve(tR, null[1:r])
            zcov = tQ %*% object@V %*% t(tQ)
            F = try(sum(z * solve(zcov, z)) / r)
            if (inherits(F, "try-error"))
                c(df1 = r, df2 = NA,  F.ratio = NA, p.value = NA, note = 1)
            else {
                df2 = object@dffun(tQ, object@dfargs)
                if (is.na(df2))
                    p.value = pchisq(F*r, r, lower.tail = FALSE)
                else
                    p.value = pf(F, r, df2, lower.tail = FALSE)
                c(round(c(df1 = r, df2 = df2), 2), 
                  F.ratio = round(F, 3), p.value = p.value, note = rrflag)
            }
        })
        
        result = as.data.frame(t(as.data.frame(result)))
        if (!missing(by)) {
            fbr = sapply(by.rows, "[", 1)
            result = cbind(object@grid[fbr, by, drop = FALSE], result)
        }
        class(result) = c("summary_emm", "data.frame")
        attr(result, "estName") = "F.ratio"
        if (!status && all(result$note == 0))
            result$note = NULL
        else {
            if (any(result$note %in% c(1,3)))
                attr(result, "mesg") = .dep.msg
            if (any(result$note %in% c(2,3)))
                attr(result, "mesg") = c(attr(result, "mesg"), .est.msg)
            result$note = sapply(result$note, function(x) 
                switch(x + 1, "", " d", "   e", " d e"))
        }
        result
    }
}

# messages (also needed by joint_tests())
.dep.msg = "d: df1 reduced due to linear dependence"
.est.msg = "e: estimability issue may distort result"

# Do all joint tests of contrasts. by, ... passed to emmeans() calls

#' Compute \dQuote{type III} tests of the terms in a model
#'
#' This function produces an analysis-of-variance-like table based on linear
#' functions of predictors in a model or \code{emmGrid} object. Specifically,
#' the function constructs, for each combination of factors, a set of
#' (interaction) contrasts via \code{\link{contrast}}, and then tests them using
#' \code{\link{test}} with \code{joint = TRUE}. Optionally, one or more of the
#' predictors may be used as a \dQuote{by} variable, so that separate tables of
#' tests are produced for each combination of them.
#' 
#' Usually, but not always, these tests correspond to \dQuote{type III} tests
#' a la \pkg{SAS} (see more details later). It is more robust
#' than the model-reduction approach used in \code{\link[car]{Anova}} in the
#' \pkg{car} package. In the latter, the user must fit the model using
#' sum-to-zero contrast specifications such as \code{contr.sum} in order to
#' obtain correct type-III tests. But \code{joint_tests} does them correctly
#' regardless of which contrast functions were used to construct the model
#' matrix. 
#' 
#' While \code{joint_tests} always tests contrasts of EMMs,
#' \pkg{SAS}'s type-III tests are not always constructed that way. 
#' For example, if a factor interacts with a covariate, contrasts of 
#' the factor EMMs give nonzero weights to the interaction effects
#' whereas SAS's type-III estimable functions give zero weight to the 
#' interactions; thus, the joint test of the EMM contrasts is not a 
#' type-III test. In short, type-III tests are tests of terms in a 
#' model, while \code{joint_tests} tests contrasts of EMMs. See the 
#' illustration in the examples below.
#' 
#' 
#' @param object a fitted model or an \code{emmGrid}. If a fitted model, it is
#'    replaced by \code{ref_grid(object, cov.reduce = range, ...)}
#' @param by character names of \code{by} variables. Separate sets of tests are
#'    run for each combination of these.
#' @param ... additional arguments passed to \code{ref_grid} and \code{emmeans}
#'
#' @return a \code{summary_emm} object (same as is produced by 
#'   \code{\link{summary.emmGrid}}). All effects for which there are no
#'   estimable contrasts are omitted from the results.
#'   
#' @seealso \code{\link{test}}
#' @export
#'
#' @examples
#' pigs.lm <- lm(log(conc) ~ source * factor(percent), data = pigs)
#' 
#' joint_tests(pigs.lm)                     ## type III ANOVA
#' 
#' joint_tests(pigs.lm, weights = "outer")  ## differently weighted
#' 
#' joint_tests(pigs.lm, by = "source")      ## separate type III tests of 'percent'
#' 
#' ### Illustration of discrepancies with type III tests
#' toydf = data.frame(
#'     treat = rep(c("A", "B"), c(4, 6)),
#'     female = c(1, 0, 0, 1,   0, 0, 0, 1, 1, 0 ),
#'     resp = c(17, 12, 14, 19, 28, 26, 26, 34, 33, 27))
#' options(contrasts = c("contr.sum", "contr.poly"))
#' toy.lmc = lm(resp ~ treat * female, data = toydf)
#' toy.lmf = lm(resp ~ treat * factor(female), data = toydf)
#' # (These two models have identical fitted values and residuals)
#' 
#' joint_tests(toy.lmc)           # table 1
#' joint_tests(toy.lmf)           # table 2
#' car::Anova(toy.lmc, type = 3)  # table 3
#' car::Anova(toy.lmf, type = 3)  # table 4
#' ## Observe that tables 1 and 2, and 4 are the same, but 
#' ##   table 3 is different in its test for 'treat'
#' ## Exercise: Compare the F ratios with the squares of the t ratios
#' ##           of these models' coefficients.
joint_tests = function(object, by = NULL, ...) {
    if (!inherits(object, "emmGrid"))
        object = ref_grid(object, cov.reduce = range, ...)
    facs = setdiff(names(object@levels), by)
    do.test = function(these, facs, result, ...) {
        if ((k <- length(these)) > 0) {
            emm = emmeans(object, these, by = by, ...)
            tst = test(contrast(emm, interaction = "consec"), 
                       joint = TRUE, status = TRUE)
            tst = cbind(ord = k, `model term` = paste(these, collapse = ":"), tst)
            result = rbind(result, tst)
            last = max(match(these, facs))
        }
        else
            last = 0
        if (last < (n <- length(facs)))
            for (i in last + seq_len(n - last))
                result = do.test(c(these, facs[i]), facs, result, ...)
        result
    }
    result = suppressMessages(do.test(character(0), facs, NULL, ...))
    result = result[order(result[[1]]), -1, drop = FALSE]
    result = result[result$df1 > 0, , drop = FALSE]
    class(result) = c("summary_emm", "data.frame")
    attr(result, "estName") = "F.ratio"
    attr(result, "by.vars") = by
    if (any(result$note != "")) {
        if (any(result$note %in% c("d", "d e")))  
            msg = .dep.msg
        else 
            msg = character(0)
        if (any(result$note %in% c("  e", "d e")))
            msg = c(msg, .est.msg)
        attr(result, "mesg") = msg
    }
    else
        result$note = NULL
    result
}

# provide for displaying in standard 'anova' format (with astars etc.)
# I'm not going there now. Maybe later, probably not

#' #' @export
#' as.anova = function(object, ...)
#'     UseMethod("as.anova")
#' 
#' as.anova.summary_emm = function(object, ...) {
#'     class(object) = c("anova", "data.frame")
#'     row.names(object) = as.character(object[[1]])
#'     names(object) = gsub("p.value", "Pr(>F)", names(object))
#'     object[-1]
#' }


