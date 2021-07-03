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
#' @param mult.name Name of the factor that defines the multivariate response
#' @param null Scalar or conformable vector of null-hypothesis values to test against
#' @param by Any \code{by} variable(s). These should not include the primary variables
#'   to be contrasted.
#' @param adjust Character value of a multiplicity adjustment method (\code{"none"} for no adjustment)
#' @param show.ests Logical flag determining whether the multivariate means are displayed
#' @param ... Additional arguments passed to \code{contrast}
#'
#' @return An object of class \code{summary_emm} containing the multivariate
#'   test results; or a list of the estimates and the tests if \code{show.ests}
#'   is \code{TRUE}. The test results include the Mahalanobis distances (\emph{not} squared),
#'   \eqn{F} ratios, degrees of freedom, and \eqn{P} values.
#' @note
#' While designed primarily for testing contrasts, multivariate tests of the mean
#' vector itself can be implemented via \code{method = "identity")} (see the examples).
#'
#' @references Hotelling, Harold (1931) "The generalization of Student's ratio", 
#'   \emph{Annals of Mathematical Statistics} 2(3), 360–378. doi:10.1214/aoms/1177732979
#' 
#' @export
#'
#' @examples
#' MOats.lm <- lm(yield ~ Variety + Block, data = MOats)
#' MOats.emm <- emmeans(MOats.lm, ~ Variety | rep.meas)
#' mvcontrast(MOats.emm, "consec", show.ests = TRUE)
#' 
#' # Test each mean against a specified null vector
#' mvcontrast(MOats.emm, "identity", name = "Variety", 
#'            null = c(80, 100, 120, 140), adjust = "none")
#' # (Note 'name' is passed to contrast() and overrides default name "contrast")
#' 
mvcontrast = function(object, method = "eff", mult.name = "rep.meas", null = 0,
                      by = object@misc$by.vars, adjust = c("sidak", p.adjust.methods),
                      show.ests = FALSE, ...) {
    con = contrast(object, method = method, by = union(by, mult.name), ...)
    con = contrast(con, "identity", simple = mult.name, name = "dimension.") # just re-orders it
    ese = .est.se.df(con)
    est = ese$est
    df = ese$df
    V = vcov(con)
    rows = .find.by.rows(con@grid, con@misc$by.vars)
    result = lapply(rows, function(r) {
        df1 = length(r)
        rawdf = mean(df[r])
        df2 = rawdf - df1 + 1
        D2 = mahalanobis(est[r], null, V[r, r])
        F = D2 / df1 * (df2 / rawdf)
        data.frame(Mahal.dist = sqrt(D2), df1 = df1, df2 = df2, F.ratio = F)
    })
    result = cbind(con@grid[sapply(rows, function(r) r[1]), ], do.call(rbind, result))
    result[["dimension."]] = NULL
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
    if (adjust != "none")
        attr(result, "mesg") = paste("P value adjustment:", adjust)
    if (show.ests)
        list(estimates = summary(con, by = setdiff(con@misc$by.vars, mult.name)), tests = result)
    else
        result
}