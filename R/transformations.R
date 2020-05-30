##############################################################################
#    Copyright (c) 2012-2016 Russell V. Lenth                                #
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

# Code to implement transformations my way

# Implementation of additional transformations, typically ones with parameters
# Returns a list like stats::make.link, but often with an additional "param" member
# types:
#       glog: log(mu + param)


#' Response-transformation extensions
#' 
#' The \code{make.tran} function creates the needed information to perform
#' transformations of the response
#' variable, including inverting the transformation and estimating variances of
#' back-transformed predictions via the delta method. \code{make.tran} is
#' similar to \code{\link{make.link}}, but it covers additional transformations.
#' The result can be used as an environment in which the model is fitted, or as
#' the \code{tran} argument in \code{\link{update.emmGrid}} (when the given
#' transformation was already applied in an existing model).
#' 
#' The functions \code{\link{emmeans}}, \code{\link{ref_grid}}, and related ones
#' automatically detect response transformations that are recognized by
#' examining the model formula. These are \code{log}, \code{log2}, \code{log10},
#' \code{sqrt}, \code{logit}, \code{probit}, \code{cauchit}, \code{cloglog}; as
#' well as (for a response variable \code{y}) \code{asin(sqrt(y))},
#' \code{asinh(sqrt(y))}, and \code{sqrt(y) + sqrt(y+1)}. In addition, any
#' constant multiple of these (e.g., \code{2*sqrt(y)}) is auto-detected and
#' appropriately scaled (see also the \code{tran.mult} argument in
#' \code{\link{update.emmGrid}}).
#' 
#' A few additional character strings may be supplied as the \code{tran}
#' argument in \code{\link{update.emmGrid}}: \code{"identity"},
#' \code{"1/mu^2"}, \code{"inverse"}, \code{"reciprocal"}, \code{"asin.sqrt"},
#' and \code{"asinh.sqrt"}.
#' 
#' More general transformations may be provided as a list of functions and
#' supplied as the \code{tran} argument as documented in
#' \code{\link{update.emmGrid}}. The \code{make.tran} function returns a
#' suitable list of functions for several popular transformations. Besides being
#' usable with \code{update}, the user may use this list as an enclosing
#' environment in fitting the model itself, in which case the transformation is
#' auto-detected when the special name \code{linkfun} (the transformation
#' itself) is used as the response transformation in the call. See the examples
#' below.
#' 
#' Most of the transformations available in "make.tran" require a parameter, 
#' specified in \code{param}; in the following discussion, we use \eqn{p} to
#' denote this parameter, and \eqn{y} to denote the response variable.
#' The \code{type} argument specifies the following transformations:
#' \describe{
#' \item{\code{"genlog"}}{Generalized logarithmic transformation: \eqn{log(y +
#'   p)}, where \eqn{y > -p}}
#' \item{\code{"power"}}{Power transformation: \eqn{y^p}, where \eqn{y > 0}.
#'   When \eqn{p = 0}, \code{"log"} is used instead}
#' \item{\code{"boxcox"}}{The Box-Cox transformation (unscaled by the geometric
#'   mean): \eqn{(y^p - 1) / p}, where \eqn{y > 0}. When \eqn{p = 0}, \eqn{log(y)}
#'   is used.}
#' \item{\code{"sympower"}}{A symmetrized power transformation on the whole real
#'   line:
#'   \eqn{abs(y)^p * sign(y)}. There are no restrictions on \eqn{y}, but we
#'   require \eqn{p > 0} in order for the transformation to be monotone and
#'   continuous.}
#' \item{\code{"asin.sqrt"}}{Arcsin-square-root transformation:
#'   \eqn{sin^(-1)(y/p)^{1/2)}}. Typically, the parameter \eqn{p} is equal to 1 for
#'   a fraction, or 100 for a percentage.}
#' \item{\code{"bcnPower"}}{Box-Cox with negatives allowed, as described for the 
#'   \code{bcnPower} function in the \pkg{car} package. It is defined as the Box-Cox
#'   transformation \eqn{(z^p - 1) / p} of the variable \eqn{z = y + (y^2+g^2)^(1/2)}. 
#'   This requires \code{param} to have two elements:
#'   the power \eqn{p} and the offset \eqn{g > 0}.}
#' \item{\code{"scale"}}{This one is a little different than the others, in that
#'   \code{param} is ignored; instead, \code{param} is determined by calling 
#'   \code{scale(y, ...)}. The user should give as \code{y} the response variable in the
#'   model to be fitted to its scaled version.}
#' }
#' The user may include a second element in \code{param} to specify an
#' alternative origin (other than zero) for the \code{"power"}, \code{"boxcox"},
#' or \code{"sympower"} transformations. For example, \samp{type = "power",
#' param = c(1.5, 4)} specifies the transformation \eqn{(y - 4)^1.5}.
#' In the \code{"genpower"} transformation, a second \code{param} element may be
#' used to specify a base other than the default natural logarithm. For example,
#' \samp{type = "genlog", param = c(.5, 10)} specifies the \eqn{log10(y + .5)}
#' transformation. In the \code{"bcnPower"} transformation, the second element
#' is required and must be positive.
#' 
#' For purposes of back-transformation, the \samp{sqrt(y) + sqrt(y+1)}
#' transformation is treated exactly the same way as \samp{2*sqrt(y)}, because
#' both are regarded as estimates of \eqn{2\sqrt\mu}.
#'
#' @param type The name of the transformation. See Details.
#' @param param Numeric parameter needed for the transformation. Optionally, it 
#'   may be a vector of two numeric values; the second element specifies an
#'   alternative base or origin for certain transformations. See Details.
#' @param y,... Used only with \code{type = "scale"}. These parameters are
#'   passed to \code{\link{scale}} to determine \code{param}.
#'
#' @return A \code{list} having at least the same elements as those returned by
#'   \code{\link{make.link}}. The \code{linkfun} component is the transformation
#'   itself.
#' 
#' @note The \code{genlog} transformation is technically unneeded, because
#'   a response transformation of the form \code{log(y + c)} is now auto-detected 
#'   by \code{\link{ref_grid}}.
#' @note We modify certain \code{\link{make.link}} results in transformations
#'   where there is a restriction on valid prediction values, so that reasonable
#'   inverse predictions are obtained, no matter what. For example, if a
#'   \code{sqrt} transformation was used but a predicted value is negative, the
#'   inverse transformation is zero rather than the square of the prediction. A
#'   side effect of this is that it is possible for one or both confidence
#'   limits, or even a standard error, to be zero.
#' @export
#'
#' @examples
#' # Fit a model using an oddball transformation:
#' bctran <- make.tran("boxcox", 0.368)
#' warp.bc <- with(bctran, 
#'     lm(linkfun(breaks) ~ wool * tension, data = warpbreaks))
#' # Obtain back-transformed LS means:    
#' emmeans(warp.bc, ~ tension | wool, type = "response")
#' 
#' ### Using a scaled response...
#' # Case where it is auto-detected:
#' fib.lm <- lm(scale(strength) ~ diameter + machine, data = fiber)
#' ref_grid(fib.lm) 
#' 
#' # Case where scaling is not auto-detected -- and what to do about it:
#' fib.aov <- aov(scale(strength) ~ diameter + Error(machine), data = fiber)
#' fib.rg <- suppressWarnings(ref_grid(fib.aov, at = list(diameter = c(20, 30))))
#' 
#' # Scaling was not retrieved, so we can do:
#' fib.rg = update(fib.rg, tran = make.tran("scale", y = fiber$strength))
#' emmeans(fib.rg, "diameter")

#' 
#' \dontrun{
#' ### An existing model 'mod' was fitted with a y^(2/3) transformation...
#'   ptran = make.tran("power", 2/3)
#'   emmeans(mod, "treatment", tran = ptran)
#' }
make.tran = function(type = c("genlog", "power", "boxcox", "sympower", 
                              "asin.sqrt", "bcnPower", "scale"), param = 1, y, ...) {
    type = match.arg(type)
    origin = 0
    mu.lbl = "mu"
    if (length(param) > 1) {
        origin = param[2]
        param = param[1]
        mu.lbl = paste0("(mu - ", round(origin, 3), ")")
    }
    if(type == "scale") {
        sy = scale(y, ...)
        if(is.null(origin <- attr(sy, "scaled:center")))
            origin = 0
        if(is.null(param <- attr(sy, "scaled:scale")))
            param = 1
        remove(list = c("y", "sy")) # remove baggage from env
    }
    switch(type,
           genlog = {
               if((origin < 0) || (origin == 1))
                   stop('"genlog" transformation must have a positive base != 1')
               logbase = ifelse(origin == 0, 1, log(origin))
               xlab = ifelse(origin == 0, "", paste0(" (base ", round(origin, 3), ")"))
               list(linkfun = function(mu) log(pmax(mu + param, 0)) / logbase,
                    linkinv = function(eta) pmax(exp(logbase * eta), .Machine$double.eps) - param,
                    mu.eta = function(eta) logbase * pmax(exp(logbase * eta), .Machine$double.eps),
                    valideta = function(eta) TRUE,
                    param = c(param, origin),
                    name = paste0("log(mu + ", round(param,3), ")", xlab)
               )
           },
           power = {
               if (param == 0) {
                   if(origin == 0) make.link("log")
                   else make.tran("genlog", -origin)
               }
               else list(
                   linkfun = function(mu) pmax(mu - origin, 0)^param,
                   linkinv = function(eta) origin + pmax(eta, 0)^(1/param),
                   mu.eta = function(eta) pmax(eta, 0)^(1/param - 1) / param,
                   valideta = function(eta) all(eta > 0),
                   param = c(param, origin),
                   name = ifelse(param > 0, 
                                 paste0(mu.lbl, "^", round(param,3)),
                                 paste0(mu.lbl, "^(", round(param,3), ")"))
               )
           },
           boxcox = {
               if (param == 0) {
                   result = if(origin == 0) make.link("log")
                   else make.tran("genlog", -origin)
                   return (result)
               }
               min.eta = ifelse(param > 0, -1 / param, -Inf)
               xlab = ifelse(origin == 0, "", paste0(" with origin at ", round(origin, 3)))
               list(
                   linkfun = function(mu) ((mu - origin)^param - 1) / param,
                   linkinv = function(eta) origin + (1 + param * pmax(eta, min.eta))^(1/param),
                   mu.eta = function(eta) (1 + param * pmax(eta, min.eta))^(1/param - 1),
                   valideta = function(eta) all(eta > min.eta),
                   param = c(param, origin),
                   name = paste0("Box-Cox (lambda = ", round(param, 3), ")", xlab)
               )
           },
           sympower = {
               if (param <= 0) 
                   stop('"sympower" transformation requires positive param')
               if (origin == 0) 
                   mu.lbl = paste0("(", mu.lbl, ")")
               absmu.lbl = gsub("\\(|\\)", "|", mu.lbl)
               list(linkfun = function(mu) sign(mu - origin) * abs(mu - origin)^param,
                    linkinv = function(eta) origin + sign(eta) * abs(eta)^(1/param),
                    mu.eta = function(eta) (abs(eta))^(1/param - 1),
                    valideta = function(eta) all(eta > min.eta),
                    param = c(param, origin),
                    name = paste0(absmu.lbl, "^", round(param,3), " * sign", mu.lbl)
               )
           },
           asin.sqrt = {
               mu.lbl = ifelse(param == 1, "mu", paste0("mu/", round(param,3)))
               list(linkfun = function(mu) asin(sqrt(mu/param)),
                    linkinv = function(eta) param * sin(pmax(pmin(eta, pi/2), 0))^2,
                    mu.eta = function(eta) param * sin(2*pmax(pmin(eta, pi/2), 0)),
                    valideta = function(eta) all(eta <= pi/2) && all(eta >= 0),
                    name = paste0("asin(sqrt(", mu.lbl, "))")
               )
           },
           bcnPower = {
               if(origin <= 0)
                   stop ("The second parameter for 'bcnPower' must be strictly positive.")
               list(
                   linkfun = function(mu) {
                       s = sqrt(mu^2 + origin^2)
                       if (abs(param) < 1e-10) log(.5*(mu + s))
                       else ((0.5 * (mu + s))^param - 1) / param  },
                   linkinv = function(eta) {
                       q = if (abs(param) < 1e-10) 2 * exp(eta)
                           else 2 * (param * eta + 1) ^ (1/param)
                       (q^2 - origin^2) / (2 * q) },
                   mu.eta = function(eta) {
                       if (abs(param) < 1e-10) { q = 2 * exp(eta); dq = q }
                       else { q = 2 * (param * eta + 1) ^ (1/param)
                           dq = 2 * (param * eta + 1)^(1/param - 1) }
                       0.5 * (1 + (origin/q)^2) * dq },
                   valideta = function(eta) all(eta > 0),
                   param = c(param, origin),
                   name = paste0("bcnPower(", signif(param,3), ", ", signif(origin,3), ")")
               )
           },
           scale = list(
               linkfun = function(mu) (mu - origin) / param,
               linkinv = function(eta) param * eta + origin,
               mu.eta = function(eta) rep(param, length(eta)),
               valideta = function(eta) TRUE,
               name = paste0("scale(", signif(origin, 3), ", ", signif(param, 3), ")"),
               param = c(param, origin)
           )
    )
}




### My modification/expansion of stats:make.link()
### Also, if not found, returns make.link("identity") modified with
##     unknown = TRUE, name = link
## In addition, I make all links truly monotone on (-Inf, Inf) in
##     lieu of valideta
##
## Extensions to make.link results:
##     unknown: set to TRUE if link is unknown
##     mult: scalar multiple of transformation
##
.make.link = function(link) {
    if (link %in% c("logit", "probit", "cauchit", "cloglog", "identity", "log"))
        result = stats::make.link(link)
    else result = switch(link,
         sqrt = { 
             tmp = make.link("sqrt") 
             tmp$linkinv = function(eta) pmax(0, eta)^2
             tmp$mu.eta = function(eta) 2*pmax(0, eta)
             tmp },
         `1/mu^2` = { 
             tmp = make.link("1/mu^2") 
             tmp$linkinv = function(eta) 1/sqrt(pmax(0, eta))
             tmp$mu.eta = function(eta) -1/(2*pmax(0, eta)^1.5)
             tmp },
         inverse = { 
             tmp = make.link("inverse") 
             tmp$linkinv = function(eta) 1/pmax(0, eta)
             tmp$mu.eta = function(eta) -1/pmax(0, eta)^2
             tmp },
         `/` = .make.link("inverse"),
         reciprocal = .make.link("inverse"),
         log10 = list(
             linkinv = function(eta) 10^eta,
             mu.eta = function(eta) 10^eta * log(10),
             name = "log10"
         ),
         log2 = list(
             linkinv = function(eta) 2^eta,
             mu.eta = function(eta) 2^eta * log(2),
             name = "log2"
         ),
         asin.sqrt = make.tran("asin.sqrt"),
         `asin.sqrt./` = make.tran("asin.sqrt", 100),
         asinh.sqrt = list(
             linkinv = function(eta) sinh(eta)^2,
             mu.eta = function(eta) sinh(2 * eta),
             name = "asinh(sqrt(mu))"
         ),
         exp = list(
           linkinv = function(eta) log(eta),
           mu.eta = function(eta) 1/eta,
           name = "exp"
         ),
         `+.sqrt` = {
             tmp = .make.link("sqrt")
             tmp$mult = 2
             tmp
         },
         log.o.r. = {
             tmp = make.link("log")
             tmp$name = "log odds ratio"
             tmp
         },
         
         { # default if not included, flags it as unknown
             tmp = stats::make.link("identity")
             tmp$unknown = TRUE
             tmp$name = link
             tmp
         }
    )
    result
}

