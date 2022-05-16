##############################################################################
#    Copyright (c) 2012-2021 Russell V. Lenth                                #
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

#' Multivariate contrasts
#' 
#' This function displays tests of multivariate comparisons or contrasts.
#' The contrasts are constructed at each level of the variable in \code{mult.name},
#' and then we do a multivariate test that the vector of estimates is equal to
#' \code{null} (zero by default). The \emph{F} statistic and degrees
#' of freedom are determined via the Hotelling distribution. that is, if there are
#' \eqn{m} error degrees of freedom and multivariate dimensionality \eqn{d}, then
#' the resulting \eqn{F} statistic has degrees of freedom \eqn{(d, m - d + 1)}
#' as shown in Hotelling (1931).
#' 
#' 
#' @param object An object of class \code{emmGrid}
#' @param method A contrast method, per \code{\link{contrast.emmGrid}}
#' @param mult.name Character vector of nNames of the factors whose levels
#'   define the multivariate means to contrast. If the model itself has a
#'   multivariate response, that is what is used. Otherwise, \code{mult.name}
#'   \emph{must} be specified.
#' @param null Scalar or conformable vector of null-hypothesis values to test against
#' @param by Any \code{by} variable(s). These should not include the primary
#'   variables to be contrasted. For convenience, the \code{by} variable is
#'   nulled-out if it would result in no primary factors being contrasted.
#' @param adjust Character value of a multiplicity adjustment method
#'   (\code{"none"} for no adjustment). The available adjustment methods are
#'   more limited that in \code{contrast}, and any default adjustment returned
#'   via \code{method} is ignored.
#' @param show.ests Logical flag determining whether the multivariate means 
#'   are displayed
#' @param ... Additional arguments passed to \code{contrast}
#'
#' @return An object of class \code{summary_emm} containing the multivariate
#'   test results; or a list of the estimates and the tests if \code{show.ests}
#'   is \code{TRUE}. The test results include the Hotelling \eqn{T^2} statistic,
#'   \eqn{F} ratios, degrees of freedom, and \eqn{P} values.
#' @note
#' If some interactions among the primary and \code{mult.name} factors are
#' absent, the covariance of the multivariate means is singular; this situation
#' is accommodated, but the result has reduced degrees of freedom and a message
#' is displayed. If there are other abnormal conditions such as non-estimable
#' results, estimates are shown as \code{NA}.
#'
#' While designed primarily for testing contrasts, multivariate tests of the
#' mean vector itself can be implemented via \code{method = "identity")} (see
#' the examples).
#'
#' @references Hotelling, Harold (1931) "The generalization of Student's ratio", 
#'   \emph{Annals of Mathematical Statistics} 2(3), 360–378. doi:10.1214/aoms/1177732979
#' 
#' @export
#'
#' @examples
#' MOats.lm <- lm(yield ~ Variety + Block, data = MOats)
#' MOats.emm <- emmeans(MOats.lm, ~ Variety | rep.meas)
#' mvcontrast(MOats.emm, "consec", show.ests = TRUE)  # mult.name defaults to rep.meas
#' 
#' # Test each mean against a specified null vector
#' mvcontrast(MOats.emm, "identity", name = "Variety", 
#'            null = c(80, 100, 120, 140), adjust = "none")
#' # (Note 'name' is passed to contrast() and overrides default name "contrast")
#' 
#' # 'mult.name' need not refer to a multivariate response
#' mvcontrast(MOats.emm, "trt.vs.ctrl1", mult.name = "Variety")
#' 
mvcontrast = function(object, method = "eff", mult.name = object@roles$multresp, null = 0,
                      by = object@misc$by.vars, adjust = c("sidak", p.adjust.methods),
                      show.ests = FALSE, ...) {
    object = .chk.list(object)
    if (is.null(mult.name) || length(mult.name) == 0)
        stop("Must specify at least one factor in 'mult.name'")
    if(length(setdiff(names(object@levels), union(by, mult.name))) == 0)
        by = NULL   # avoid the case where we're left with no variables
    con = contrast(object, method = method, by = union(by, mult.name), ...)
    mvnm = paste(mult.name, collapse = " ")
    con = contrast(con, "identity", simple = mult.name, name = mvnm) # just re-orders it
    ese = .est.se.df(con)
    est = ese$est
    df = ese$df
    V = vcov(con)
    rows = .find.by.rows(con@grid, con@misc$by.vars)
    red.rank = FALSE  # flag for red.rank cases
    result = lapply(rows, function(r) {
        QR = try(qr(V[r, r, drop = FALSE]), silent = TRUE)
        if (inherits(QR, "try-error"))
            T2 = df1 = df2 = F = NA
        else {
            df1 = QR$rank
            if(df1 < length(r)) red.rank <<- TRUE
            rawdf = mean(df[r])
            df2 = rawdf - df1 + 1
            qe = qr.coef(QR, est[r] - null)
            qe[is.na(qe)] = 0
            T2 = sum(qe * (est[r] - null))
            F = T2 / df1 * (df2 / rawdf)
        }
        data.frame(T.square = T2, df1 = df1, df2 = df2, F.ratio = F)
    })
    result = cbind(con@grid[sapply(rows, function(r) r[1]), ], do.call(rbind, result))
    result[[mvnm]] = NULL
    class(result) = c("summary_emm", "data.frame")

    by = setdiff(by, mult.name)
    if (length(by) == 0) 
        by = NULL
    rows = .find.by.rows(result, by)
    adjust = match.arg(adjust)
    result$p.value = NA
    for (r in rows) {
        pv = with(result[r, ], pf(F.ratio, df1, df2, lower.tail = FALSE))
        result$p.value[r] = switch(adjust, sidak = 1 - (1 - pv)^length(r),
            p.adjust(pv, adjust))
    }
    attr(result, "estName") = "F.ratio"
    attr(result, "by.vars") = by
    if (adjust == "none") 
        mesg = NULL
    else
        mesg = paste("P value adjustment:", adjust)
    if (red.rank)
        mesg = c(mesg, "NOTE: Some or all d.f. are reduced due to singularities")
    if(any(is.na(result$T.square)))
        mesg = c(mesg, 
            "NAs indicate non-estimabile cases or other errors")
    attr(result, "mesg") = mesg
    if (show.ests)
        list(estimates = con, tests = result)
    else
        result
}