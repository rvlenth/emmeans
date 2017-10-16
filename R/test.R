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
#' @method test emmGrid
#' @export
test.emmGrid = function(object, null = 0, 
                    joint = FALSE, verbose = FALSE, rows, by, ...) {
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
        
        ### L = L[, estble.idx, drop = FALSE]
        if (!missing(rows))
            by.rows = list(sel.rows = rows)
        else {
            by.rows = list(all = seq_len(nrow(L)))
            if(missing(by)) 
                by = object@misc$by.vars
            if (!is.null(by)) 
                by.rows = .find.by.rows(object@grid, by)
        }
        
        lindep = nonest = FALSE
        
        result = lapply(by.rows, function(rows) {
            LL = L[rows, , drop = FALSE]
            # estract est'ble rows
            if(est.flag) {
                erows = estimability::is.estble(LL, object@nbasis)
                nonest <<- nonest || (sum(erows) < nrow(LL))
                LL = LL[erows, estble.idx, drop = FALSE]
            }
            # Check rank
            qrLt = qr(t(LL)) # this will work even if LL has 0 rows
            r = qrLt$rank
            if (r == 0)
                return(c(df1 = 0, df2 = NA, F.ratio = NA, p.value = NA))
            if (r < nrow(LL)) {
                if(!all(null == 0))
                    stop("Rows are linearly dependent - cannot do the test when 'null' != 0")
                else 
                    lindep <<- TRUE
            }
            tR = t(qr.R(qrLt))[1:r, 1:r, drop = FALSE]
            tQ = t(qr.Q(qrLt))[1:r, , drop = FALSE]
            if(length(null) < r) null = rep(null, r)
            z = tQ %*% bhat - solve(tR, null[1:r])
            zcov = tQ %*% object@V %*% t(tQ)
            F = sum(z * solve(zcov, z)) / r
            df2 = object@dffun(tQ, object@dfargs)
            if (is.na(df2))
                p.value = pchisq(F*r, r, lower.tail = FALSE)
            else
                p.value = pf(F, r, df2, lower.tail = FALSE)
            c(round(c(df1 = r, df2 = df2), 2), F.ratio = round(F, 3), p.value = p.value)
        })
        
        result = as.data.frame(t(as.data.frame(result)))
        if (!missing(by)) {
            fbr = sapply(by.rows, "[", 1)
            result = cbind(object@grid[fbr, by, drop = FALSE], result)
        }
        class(result) = c("summary_emm", "data.frame")
        attr(result, "estName") = "F.ratio"
        if (lindep)
            message("There are linearly dependent rows - df are reduced accordingly")
        if (nonest)
            message("Some rows are non-estimable and were excluded")
        
        result
    }
}
