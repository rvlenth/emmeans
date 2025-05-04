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

# summary of a summary - just pass it thru
#' @export
summary.summary_emm = function(object, ...) object

# confint and test methods

# confint method
#' @rdname summary.emmGrid
#' @order 2
#' @param parm (Required argument for \code{confint} methods, but not used)
#' @method confint emmGrid
#' @export
confint.emmGrid = function(object, parm, level = .95, ...) {
    summary(object, infer = c(TRUE, FALSE), level = level, ...)
}

#' @rdname summary.emmGrid
#' @order 3
#' @export
test = function(object, null, ...) {
    UseMethod("test")
}

#' @rdname summary.emmGrid
#' @order 3
#' @param joint Logical value. If \code{FALSE}, the arguments are passed to 
#'   \code{\link{summary.emmGrid}} with \code{infer=c(FALSE, TRUE)}. If \code{joint = 
#'   TRUE}, a joint test of the hypothesis L beta = null is performed, where L 
#'   is \code{linfct(object)} and beta is the vector of fixed effects estimated 
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
        
        # my own zapsmall fcn - sets a hard limit
        zapsm = function(x, tol = 1e-7) {
            x[abs(x) < tol] = 0
            x
        }

        result = lapply(by.rows, function(rows) {
            LL = L[rows, , drop = FALSE]

            # discard any rows that have NAs
            narows = apply(LL, 1, function(x) any(is.na(x)) | all(x == 0))
            LL = LL[!narows, , drop = FALSE]
            rrflag = 0 + 2 * any(narows)  ## flag for estimability issue
            
            if(est.flag)  { 
                if (any(!estimability::is.estble(LL, nbasis, est.tol))) {
                    LL = estimability::estble.subspace(zapsm(LL), nbasis)
                    rrflag = bitwOr(rrflag, 2)
                }
            }
            # Check rank
            qrLt = qr(zapsm(t(LL))) # this will work even if LL has 0 rows
            r = qrLt$rank
            if (r == 0)
                return(c(df1 = 0, df2 = NA, F.ratio = NA, p.value = NA, note = 3))
            
            if (r < nrow(LL)) {
                if(!all(null == 0))
                    stop("Rows are linearly dependent - cannot do the test when 'null' != 0")
                rrflag = bitwOr(rrflag, 1)
            }
            nz = 1:r   # indexes to use
            tR = t(qr.R(qrLt))[nz, nz, drop = FALSE]
            # for some reason I am getting overly optimistic ranks sometimes...
            nz = which(zapsmall(diag(tR)) != 0)
            if (length(nz) < r) {
                r = length(nz)
                tR = tR[nz, nz, drop = FALSE]
            }
            if(length(null) < r) null = rep(null, r)
            tQ = tQ.all = t(qr.Q(qrLt))[nz, , drop = FALSE]
            # tQ.all will have all the columns. tQ may get subsetted
            # NOW get rid of the NA parts...
            if (est.flag)
                tQ = tQ[, estble.idx, drop = FALSE]
            z = try(tQ %*% bhat - solve(tR, null[nz]), silent = TRUE)
            zcov = tQ %*% object@V %*% t(tQ)
            F = try(sum(z * solve(zcov, z)) / r)
            if (inherits(F, "try-error"))
                c(df1 = r, df2 = NA,  F.ratio = NA, p.value = NA, note = 1)
            else {
                if(is.null(df2 <- object@misc$df))
                    df2 = min(apply(tQ, 1, function(.) object@dffun(., object@dfargs)))
                if (is.na(df2)) 
                    df2 = Inf
                p.value = suppressWarnings(pf(F, r, df2, lower.tail = FALSE))
                rtn = c(round(c(df1 = r, df2 = df2), 2), 
                        F.ratio = round(F, 3), p.value = p.value, note = rrflag)
                # Note: Following will screw-up joint_tests if some df2's are finite and others infinite
                if (is.infinite(df2)) {
                    rtn = c(rtn[1:3], Chisq = round(F*r, 3), rtn[4:5])
                }
                attr(rtn, "L") = tQ.all
                rtn
            }
        })
        
        ef = lapply(result, function(r) attr(r, "L"))
        if (length(unique(sapply(result, length))) > 1) { # some results don't have chisq
            result = lapply(result, \(x) {
                if (length(x) == 5) 
                    x = c(x[1:3], Chisq = NA, x[4:5])
                x
            })
        }
        result = as.data.frame(t(as.data.frame(result)))
        if (!missing(by)) {
            fbr = sapply(by.rows, "[", 1)
            result = cbind(object@grid[fbr, by, drop = FALSE], result)
        }
        class(result) = c("summary_emm", "data.frame")
        attr(result, "estName") = "F.ratio"
        attr(result, "est.fcns") = lapply(ef, zapsmall)        
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
#' In models with only factors, no covariates, these tests correspond to
#' \dQuote{type III} tests a la \pkg{SAS}, as long as equal-weighted averaging
#' is used and there are no estimability issues. When covariates are present and
#' they interact with factors, the results depend on how the covariate is
#' handled in constructing the reference grid. See the section on covariates
#' below. The point that one must always remember is that \code{joint_tests}
#' always tests contrasts among EMMs, in the context of the reference grid,
#' whereas SAS's type III tests are tests of model coefficients -- which may or may
#' not have anything to do with EMMs or contrasts.
#' 
#' @param object a fitted model, \code{emmGrid}, or \code{emm_list}. If the
#'   latter, its first element is used.
#' @param cov.reduce a function.
#'    If \code{object} is a fitted model, it is
#'    replaced by \code{ref_grid(object, cov.reduce = cov.reduce, ...)}.
#'    For this purpose, the functions \code{meanint} and \code{symmint} are
#'    available for returning an interval around the mean or around zero,
#'    respectively. Se the section below on covariates.
#' @param by character names of \code{by} variables. Separate sets of tests are
#'    run for each combination of these.
#' @param show0df logical value; if \code{TRUE}, results with zero numerator
#'    degrees of freedom are displayed, if \code{FALSE} they are skipped
#' @param showconf logical value.
#'    When we have models with estimability issues (e.g., missing cells), then with
#'    \code{showconf = TRUE}, we test any remaining effects that are not purely
#'    due to contrasts of a single term. If found, they are labeled
#'    \code{(confounded)}. See
#'    \code{vignette("xplanations")} for more information.
#' @param ... additional arguments passed to \code{ref_grid} and \code{emmeans}
#'
#' @return a \code{summary_emm} object (same as is produced by 
#'   \code{\link{summary.emmGrid}}). All effects for which there are no
#'   estimable contrasts are omitted from the results. 
#'   There may be an additional row named \code{(confounded)} which accounts
#'   for additional degrees of freedom for effects not accounted for in the 
#'   preceding rows.
#'   
#'   The returned object also includes an \code{"est.fcns"} attribute, which is a
#'   named list containing the linear functions associated with each joint test.
#'   Each row of these is standardized to have length 1. 
#'   No estimable functions are included for confounded effects.
#'   
#' @section Dealing with covariates:
#' A covariate (or any other predictor) must have \emph{more than one value in 
#' the reference grid} in order to test its effect and be included in the results.
#' Therefore, when \code{object} is a model, we default to \code{cov.reduce = meanint}
#' which sets each covariate at a symmetric interval about its mean. But
#' when \code{object} is an existing reference grid, it often has only one value 
#' for covariates, in which case they are excluded from the joint tests.
#' 
#' While having two points is sufficient when the covariate term has a linear trend,
#' you need more than two when some kind of curved trend (polynomial, spline, etc.)
#' is present -- else \code{joint_tests()} will not show enough degrees of freedom 
#' for terms involving the covariate. You may specify these points manually using \code{at},
#' or by including an \code{npts} argument in \code{cov.reduce}, via \code{make.meanint}
#' or \code{make.symmint()}. With some kinds of curved trends, the joint tests of
#' covariate terms may become somewhat meaningless.
#' 
#' Covariates present further complications in that their values in the
#' reference grid can affect the joint tests of \emph{other} effects. When
#' covariates are centered around their means (the default), then the tests we
#' obtain can be described as joint tests of covariate-adjusted means; and that
#' is our intended use here. However, some software such as \pkg{SAS} and
#' \code{car::Anova} adopt the convention of centering covariates around zero;
#' and for that purpose, one can use \code{cov.reduce = symmint(0)} when calling
#' with a model object (or in constructing a reference grid). However, adjusted
#' means with covariates set at or around zero do not make much sense in the
#' context of interpreting estimated marginal means, unless the covariate means
#' really are zero.
#' 
#' See the examples below with the \code{toy} dataset.
#' 
#' @note \code{joint_tests} is flaky with models having nested fixed effects. In
#' some cases, terms that could be relevant are not identified, or confounded
#' with unidentifiable terms.
#' 
#' @seealso \code{\link{test}}
#' @order 1
#' @export
#'
#' @examples
#' pigs.lm <- lm(log(conc) ~ source * factor(percent), data = pigs)
#' 
#' (jt <- joint_tests(pigs.lm))             ## will be same as type III ANOVA
#' 
#' ### Estimable functions associated with "percent"
#' attr(jt, "est.fcns") $ "percent"
#' 
#' joint_tests(pigs.lm, weights = "outer")  ## differently weighted
#' 
#' joint_tests(pigs.lm, by = "source")      ## separate joint tests of 'percent'
#' 
#' ### Comparisons with type III tests in SAS
#' toy = data.frame(
#'     treat = rep(c("A", "B"), c(4, 6)),
#'     female = c(1, 0, 0, 1,   0, 0, 0, 1, 1, 0 ),
#'     resp = c(17, 12, 14, 19, 28, 26, 26, 34, 33, 27))
#' toy.fac = lm(resp ~ treat * factor(female), data = toy)
#' toy.cov = lm(resp ~ treat * female, data = toy)
#' # (These two models have identical fitted values and residuals)
#' 
#' # -- SAS output we'd get with toy.fac --
#' ## Source          DF    Type III SS    Mean Square   F Value   Pr > F
#' ## treat            1    488.8928571    488.8928571    404.60   <.0001
#' ## female           1     78.8928571     78.8928571     65.29   0.0002
#' ## treat*female     1      1.7500000      1.7500000      1.45   0.2741
#' # 
#' # -- SAS output we'd get with toy.cov --
#' ## Source          DF    Type III SS    Mean Square   F Value   Pr > F
#' ## treat            1    252.0833333    252.0833333    208.62   <.0001
#' ## female           1     78.8928571     78.8928571     65.29   0.0002
#' ## female*treat     1      1.7500000      1.7500000      1.45   0.2741
#' 
#' joint_tests(toy.fac)
#' joint_tests(toy.cov)   # female is regarded as a 2-level factor by default
#' 
#' ## Treat 'female' as a numeric covariate (via cov.keep = 0)
#' ## ... then tests depend on where we center things
#' 
#' # Center around the mean
#' joint_tests(toy.cov, cov.keep = 0, cov.reduce = make.meanint(delta = 1))
#' # Center around zero (like SAS's results for toy.cov)
#' joint_tests(toy.cov, cov.keep = 0, cov.reduce = make.symmint(ctr = 0, delta = 1))
#' # Center around 0.5 (like SAS's results for toy.fac)
#' joint_tests(toy.cov, cov.keep = 0, cov.reduce = range)
#' 
#' ### Example with empty cells and confounded effects
#' low3 <- unlist(attr(ubds, "cells")[1:3]) 
#' ubds.lm <- lm(y ~ A*B*C, data = ubds, subset = -low3)
#' 
#' # Show overall joint tests by C:
#' ref_grid(ubds.lm, by = "C") |> contrast("consec") |> test(joint = TRUE)
#' 
#' # Break each of the above into smaller components:
#' joint_tests(ubds.lm, by = "C")
#' 
joint_tests = function(object, by = NULL, show0df = FALSE, 
                       showconf = TRUE,
                       cov.reduce = make.meanint(1), ...) {
    
    # hidden defaults for contrast methods and which basis to use for all contrasts
    use.contr = (function(use.contr = c("consec", "consec"), ...) use.contr)(...)

    object = .chk.list(object,...)
    if (!inherits(object, "emmGrid")) {
        args = .zap.args(object = object, cov.reduce = cov.reduce, ..., omit = "submodel")
        object = do.call(ref_grid, args)
    }
    facs = setdiff(names(object@levels), c(by, "1"))

    if(length(facs) == 0)
        stop("There are no factors to test")
    
    # Use "factors" attr if avail to screen-out interactions not in model
    # For any factors not in model (created by emmeans fcns), assume they interact w/ everything
    trmtbl = attr(object@model.info$terms, "factors")
    if (is.null(trmtbl) || (length(trmtbl) == 0))
        trmtbl = matrix(1, nrow = length(facs), dimnames = list(facs, NULL))
    else
        row.names(trmtbl) = sapply(row.names(trmtbl), function(x)
            .all.vars(reformulate(x)))
    xtras = setdiff(facs, row.names(trmtbl))
    if (length(xtras) > 0) {
        xt = matrix(1, nrow = length(xtras), ncol = ncol(trmtbl), dimnames = list(xtras, NULL))
        trmtbl = rbind(trmtbl, xt)
    }
    nesting = object@model.info$nesting
    for (nst in names(nesting)) { # make sure all complete nests are in the table
        trmtbl = cbind(trmtbl, 0)
        n = ncol(trmtbl)
        trmtbl[c(nst, nesting[[nst]]), n] = 1
    }
    est.fcns = list()
    ef.ord = character(0)

    do.test = function(these, facs, result, ...) {
        if ((k <- length(these)) > 0) {
            if(any(apply(trmtbl[these, , drop = FALSE], 2, prod) != 0)) { # term is in model
                nesters = NULL
                if (!is.null(nesting)) {
                    nst = intersect(these, names(nesting))
                    if (length(nst) > 0)
                        nesters = unlist(nesting[nst]) # proceed only if these includes all nesters
                }
                if (is.null(nesting) || length(setdiff(nesters, these)) == 0) {   
                    emm = emmeans(object, these, by = by, ...)
                    tst = test(contrast(emm, interaction = use.contr[1], by = union(by, nesters)), 
                               by = by, joint = TRUE, status = TRUE)
                    ef = attr(tst, "est.fcns") # get this before we subset the results
                    tst = tst[names(tst) != "Chisq"]   # could have some with Chisq and some without
                    mt = paste(these, collapse = ":")
                    if (length(ef) > 1)
                        ef = list(ef)
                    names(ef) = mt
                    est.fcns <<- c(est.fcns, ef)
                    ef.ord <<- c(ef.ord, k)
                    tst = cbind(ord = k, `model term` = mt, tst)
                    result = rbind(result, tst)
                }
            }
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
    if (!is.null(result) && all(is.infinite(result$df2) | is.na(result$df2))) {
        w = which(names(result) == "F.ratio")
        result = cbind(result[, 1:w], Chisq = result$F.ratio * result$df1, 
                       result[, (w+1):ncol(result)])
    }
    
    
    ## look at as-yet-unexplained effects
    if (showconf) {
        tmp = contrast(object, use.contr[2], by = by, name = ".cnt.", ...)
        br = .find.by.rows(tmp@grid, by)
        ef = est.fcns
        if (!is.null(ef)) {
            lf = rows = NULL
            for (i in seq_along(br)) {  # assemble est fcns for each by variable
                r = br[[i]]
                nm = names(br)[i]
                efi = if (!is.null(nm)) lapply(ef, function(e) e[[nm]])
                else ef
                efi = do.call(rbind, efi)
                if(!is.null(efi)) {
                    lf = rbind(lf, efi)     # stack 'em up into lf
                    rows = c(rows, r[seq_len(nrow(efi))])   # rows w/ same by combs
                }
                else {
                    lf = rbind(lf, NA * tmp@linfct[r[1], ])
                    rows = c(rows, r[1])
                }
            }
            tmpe = tmp
            tmpe@linfct = lf
            tmpe@grid = tmp@grid[rows, , drop = FALSE]
            
            ref = conf = test(tmp, joint = TRUE, status = TRUE)
            tst = test(tmpe, joint = TRUE)
            conf$df1 = ref$df1 - tst$df1
            conf$F.ratio = (ref$df1 * ref$F.ratio - tst$df1 * tst$F.ratio) / conf$df1
            conf$p.value = pf(conf$F.ratio, conf$df1, conf$df2, lower.tail = FALSE)
            conf$note = ""
            conf = cbind(ord = 999, `model term` = "(confounded)", conf)
            result = rbind(result, conf)
        }
    }
    
    result = result[order(result[[1]]), -1, drop = FALSE]
    est.fcns = est.fcns[order(ef.ord)]
    if(!show0df) {
        result = result[result$df1 > 0, , drop = FALSE]
        if(!is.null(by))
            est.fcns = lapply(est.fcns, function(x) x[!sapply(x, is.null)])
        est.fcns = est.fcns[!sapply(est.fcns, is.null)]
    }
        
    class(result) = c("summary_emm", "data.frame")
    attr(result, "estName") = "F.ratio"
    attr(result, "by.vars") = by
    nms = colnames(object@linfct)
    if(is.null(by)) 
        est.fcns = lapply(est.fcns, function(x) {
            if(!is.null(x)) colnames(x) = nms; x})
    else 
        est.fcns = lapply(est.fcns, function(L) lapply(L, function(x) {
            if(!is.null(x)) colnames(x) = nms; x}))
    attr(result, "est.fcns") = est.fcns
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


#' @rdname joint_tests
#' @order 6
#'
#' @param delta,ctr arguments for \code{make.meanint} and \code{make.symmint}.
#'   \code{delta} sets the distance each side of the center, so that the
#'   width of the interval is \code{2*delta}.
#' @param npts  number of points to include in the interval
#'
#' @return \code{make.meanint} returns the function 
#' \code{function(x) mean(x) + delta * c(-1, 1)},
#'   and \code{make.symmint(ctr, delta)} returns the function
#' \code{function(x) ctr + delta * c(-1, 1)}
#'         (which does not depend on \code{x}).
#' The cases with \code{delta = 1}, \code{meanint = make.meanint(1)} and
#' \code{symmint(ctr) = make.symmint(ctr, 1)} are retained for
#' back-compatibility reasons. These functions are available primarily for use
#' with \code{cov.reduce}.
#' @export
make.meanint = function(delta = 1, npts = 2) 
    function(x) mean(x) + delta * seq(-1, 1, length.out = npts)

#' @rdname joint_tests
#' @order 7
#' @param x argument for \code{meanint} and \code{symmint}
#' @export
meanint = function(x) mean(x) + c(-1, 1)


# const pm 1
#' @rdname joint_tests
#' @order 8
#' @export
make.symmint = function(ctr, delta = 1, npts = 2) {
    function(x) ctr + delta * seq(-1, 1, length.out = npts) 
}

#' @rdname joint_tests
#' @order 9
#' @export
symmint = function(ctr) 
    make.symmint(ctr, 1)


