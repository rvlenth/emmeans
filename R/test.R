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
                df2 = min(apply(tQ, 1, function(.) object@dffun(., object@dfargs)))
                if (is.na(df2))
                    p.value = suppressWarnings(pchisq(F*r, r, lower.tail = FALSE))
                else
                    p.value = suppressWarnings(pf(F, r, df2, lower.tail = FALSE))
                rtn = c(round(c(df1 = r, df2 = df2), 2), 
                        F.ratio = round(F, 3), p.value = p.value, note = rrflag)
                attr(rtn, "L") = tQ.all
                rtn
            }
        })
        
        ef = lapply(result, function(r) attr(r, "L"))
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
#' @param object,cov.reduce \code{object} is a fitted model or an \code{emmGrid}. 
#'    If a fitted model, it is
#'    replaced by \code{ref_grid(object, cov.reduce = cov.reduce, ...)}
#' @param by character names of \code{by} variables. Separate sets of tests are
#'    run for each combination of these.
#' @param show0df logical value; if \code{TRUE}, results with zero numerator
#'    degrees of freedom are displayed, if \code{FALSE} they are skipped
#' @param showconf logical value; if \code{true}, we look for additional effects 
#'    that are not purely due to contrasts of a single term, and if found, are
#'    labeled \code{(confounded)}. Such effects can occur, e.g., when there are empty
#'    cells in the data. Setting this to \code{FALSE} can save some computation
#'    time, especially when there are a lot of factors. The default is based on
#'    the internal vector \code{facs}, which is the names of factors being
#'    considered.
#' @param ... additional arguments passed to \code{ref_grid} and \code{emmeans}
#'
#' @return a \code{summary_emm} object (same as is produced by 
#'   \code{\link{summary.emmGrid}}). All effects for which there are no
#'   estimable contrasts are omitted from the results. 
#'   There may be an additional row named \code{(confounded)} which accounts
#'   for joint tests of effects that are estimable but are not 
#'   determined by contrasts of any one term.
#'   The returned object also includes an \code{"est.fcns"} attribute, which is a
#'   named list containing the linear functions associated with each joint test. 
#'   The order of this list may not be the same as the order of the summary.
#'   
#' @note
#' When we have models with estimability issues (e.g., missing cells), the results
#' may in rare cases depend on what contrast method is used. 
#' The default is 
#' \code{use.contr = c("consec", "consec")}, meaning that we use \code{"consec"} comparisons
#' for constructing [interaction] contrasts for named terms, and also use \code{"consec"} contrasts
#' for constructing contrasts for confounded effects (we construct the latter overall, then
#' remove any linear dependence on the named contrasts). You may override these defaults
#' via a hidden option: specify \code{use.contr = <character 2-vector>} among the arguments.
#' Examining the \code{"est.fcns"} attribute may help understand what you are testing.
#' functions.
#'   
#' @seealso \code{\link{test}}
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
#' ### Comparisons with type III tests
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
#' # results for toy.cov depend on value we use for 'female'
#' joint_tests(toy.cov, at = list(female = 0.5))
#' joint_tests(toy.cov, cov.keep = 0)   # i.e., female = mean(toy$female)
#' joint_tests(toy.cov, at = list(female = 0))
#' 
#' ### Example with empty cells and confounded effects
#' low3 <- unlist(attr(ubds, "cells")[1:3]) 
#' ubds.lm <- lm(y ~ A*B*C, data = ubds, subset = -low3)
#' joint_tests(ubds.lm, by = "B")
#' 
joint_tests = function(object, by = NULL, show0df = FALSE, 
                       showconf = (!is.na(object@nbasis[1]) && length(facs) < 4),
                       cov.reduce = range, ...) {
    
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
                    mt = paste(these, collapse = ":")
                    ef = attr(tst, "est.fcns")
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
    
    ## look at leftover effects
    if (showconf) {
        tmp = contrast(object, use.contr[2], by = by, name = ".cnt.", ...)
        # rbind all the est.fcns
        ef = est.fcns
        # if (!is.null(by))
        #     ef = lapply(ef, function(x) do.call(rbind, x))
        # ef = do.call(rbind, ef)
        if (!is.null(ef)) {
            lf = tmp@linfct
            br = .find.by.rows(tmp@grid, by)
            for (i in seq_along(br)) {
                r = br[[i]]
                nm = names(br)[i]
                tst = test(tmp[r], by = NULL, joint = TRUE)
                efi = if (!is.null(nm)) lapply(ef, function(e) e[[nm]])
                else ef
                efi = do.call(rbind, efi)
                
                ## which should we regress on efi? I'm not sure. But they're different
                allcon = attr(tst, "est.fcns")[[1]]
                
                lfi = t(qr.resid(qr(t(efi)), t(allcon)))
                slf = .squash(lfi, nbasis = tmp@nbasis)
                k = seq_len(nrow(slf))
                if (nrow(slf) == 0)
                    lf[r, ] = NA
                else {
                    lf[r[k], ] = slf
                    lf[r[-k], ] = NA
                }
            }
            k = !is.na(lf[,1])
            tmp@linfct = lf[k, , drop = FALSE]
            tmp@grid = tmp@grid[k, , drop = FALSE]
        }
        if (nrow(tmp@linfct) > 0) {
            jt = test(tmp, by = by, joint = TRUE, status = TRUE)
            tst = cbind(ord = 999, `model term` = "(confounded)", jt)
            result = rbind(result, tst)
            ef = lapply(attr(jt, "est.fcns"), zapsmall)
            if (length(ef) > 1)
                ef = list(ef)
            names(ef) = "(confounded)"
            est.fcns = c(est.fcns, ef)
            ef.ord = c(ef.ord, 999)
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


### Squash rows of L to best nrows of row space
### If nrows not specified, determine via those with SVs > tol
.squash = function(L, nbasis, tol = 1e-3, nrow) {
    if(!missing(nbasis))
        L = estimability::estble.subspace(L, nbasis)
    svd = try(svd(L, nu = 0), silent = TRUE)
    if(inherits(svd, "try-error")) return(L * 0)
    if (missing(nrow))
        nrow = sum(svd$d > tol)
    r = seq_len(nrow)
    diag(svd$d[r], nrow = nrow) %*% t(svd$v[, r, drop = FALSE])
}

