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
#'   available. See also \code{\link{joint_tests}}.
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
        if(est.flag) {
            est.tol = get_emm_option("estble.tol")
            nbasis = zapsmall(object@nbasis)
        }

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
                if (any(!estimability::is.estble(LL, nbasis, est.tol))) {
                    LL = estimability::estble.subspace(zapsmall(LL), nbasis)
                    rrflag = bitwOr(rrflag, 2)
                }
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
                df2 = min(apply(tQ, 1, function(.) object@dffun(., object@dfargs)))
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
.est.msg = "e: df1 reduced due to non-estimability"

# Do all joint tests of contrasts. by, ... passed to emmeans() calls

#' Compute joint tests of the terms in a model
#'
#' This function produces an analysis-of-variance-like table based on linear
#' functions of predictors in a model or \code{emmGrid} object. Specifically,
#' the function constructs, for each combination of factors (or covariates
#' reduced to two or more levels), a set of (interaction) contrasts via
#' \code{\link{contrast}}, and then tests them using \code{\link{test}} with
#' \code{joint = TRUE}. Optionally, one or more of the predictors may be used as
#' \code{by} variable(s), so that separate tables of tests are produced for
#' each combination of them.
#' 
#' In models with only factors, no covariates, we believe these tests correspond
#' to \dQuote{type III} tests a la \pkg{SAS}, as long as equal-weighted
#' averaging is used and there are no estimability issues. When covariates are
#' present and interact with factors, the results depend on how the covariate is
#' handled in constructing the reference grid. See the example at the end of
#' this documentation. The point that one must always remember is that
#' \code{joint_tests} always tests contrasts among EMMs, in the context of the
#' reference grid, whereas type III tests are tests of model coefficients --
#' which may or may not have anything to do with EMMs or contrasts.
#' 
#' @param object a fitted model or an \code{emmGrid}. If a fitted model, it is
#'    replaced by \code{ref_grid(object, cov.reduce = range, ...)}
#' @param by character names of \code{by} variables. Separate sets of tests are
#'    run for each combination of these.
#' @param show0df logical value; if \code{TRUE}, results with zero numerator
#'    degrees of freedom are displayed, if \code{FALSE} they are skipped
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
#' joint_tests(pigs.lm)                     ## will be same as type III ANOVA
#' 
#' joint_tests(pigs.lm, weights = "outer")  ## differently weighted
#' 
#' joint_tests(pigs.lm, by = "source")      ## separate joint tests of 'percent'
#' 
#' ### Comparisons with type III tests
#' toy = data.frame(
#'     treat = rep(c("A", "B"), c(4, 6)),
#'     female = c(1, 0, 0, 1,   0, 0, 0, 1, 1, 0 ),
#'     resp = c(17, 12, 14, 19, 28, 26, 26, 34, 33, 27))
#' toy.fac = lm(resp ~ treat * factor(female), data = toy)
#' toy.cov = lm(resp ~ treat * female, data = toy)
#' # (These two models have identical fitted values and residuals)
#' 
#' joint_tests(toy.fac)
#' joint_tests(toy.cov)   # female is regarded as a 2-level factor by default
#' 
#' joint_tests(toy.cov, at = list(female = 0.5))
#' joint_tests(toy.cov, cov.keep = 0)   # i.e., female = mean(toy$female)
#' joint_tests(toy.cov, at = list(female = 0))
#' 
#' # -- Compare with SAS output -- female as factor --
#' ## Source          DF    Type III SS    Mean Square   F Value   Pr > F
#' ## treat            1    488.8928571    488.8928571    404.60   <.0001
#' ## female           1     78.8928571     78.8928571     65.29   0.0002
#' ## treat*female     1      1.7500000      1.7500000      1.45   0.2741
#' # 
#' # -- Compare with SAS output -- female as covariate --
#' ## Source          DF    Type III SS    Mean Square   F Value   Pr > F
#' ## treat            1    252.0833333    252.0833333    208.62   <.0001
#' ## female           1     78.8928571     78.8928571     65.29   0.0002
#' ## female*treat     1      1.7500000      1.7500000      1.45   0.2741
joint_tests = function(object, by = NULL, show0df = FALSE, ...) {
    if (!inherits(object, "emmGrid"))
        object = ref_grid(object, ...)
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
    if(!show0df)
        result = result[result$df1 > 0, , drop = FALSE]
    class(result) = c("summary_emm", "data.frame")
    attr(result, "estName") = "F.ratio"
    attr(result, "by.vars") = by
    if (any(result$note != "")) {
        msg = character(0)
        if (any(result$note %in% c(" d", " d e")))  
            msg = .dep.msg
        if (any(result$note %in% c("   e", " d e")))
            msg = c(msg, .est.msg)
        attr(result, "mesg") = msg
    }
    else
        result$note = NULL
    result
}

# provide for displaying in standard 'anova' format (with astars etc.)
# I'm not g

# #' @export
# as.anova = function(object, ...)
#     UseMethod("as.anova")
# 
# as.anova.summary_emm = function(object, ...) {
#     class(object) = c("anova", "data.frame")
#     row.names(object) = as.character(object[[1]])
#     names(object) = gsub("p.value", "Pr(>F)", names(object))
#     object[-1]
# }


