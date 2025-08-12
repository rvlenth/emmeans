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
# Returns a list like stats::make.link


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
#'
#' @param type The name of a standard transformation supported by \code{stat::make.link},
#' or of a special transformation described under Details.
#' @param alpha,beta Numeric parameters needed for special transformations.
#' @param param If non-missing, this specifies either
#'   \code{alpha} or \code{c(alpha, beta)} (provided for backward compatibility).
#'   Also, for the same reason, if \code{alpha} is of length more than 1,
#'   it is taken as \code{param}.
#' @param y A numeric response variable used (\emph{and required}) with \code{type = "scale"}, 
#'   where \code{scale(y)} determines \code{alpha} and \code{beta}.
#' @param inner another transformation. See the section on compound transformations
#' @param ... Additional arguments passed to other functions/methods
#'
#'
#' @section Details:
#' The \code{make.tran} function returns a
#' suitable list of functions for several popular transformations. Besides being
#' usable with \code{update}, the user may use this list as an enclosing
#' environment in fitting the model itself, in which case the transformation is
#' auto-detected when the special name \code{linkfun} (the transformation
#' itself) is used as the response transformation in the call. See the examples
#' below.
#' 
#' The primary purpose of \code{make.tran} is to support transformations that
#' require additional parameters, specified as \code{alpha} and \code{beta};
#' these are the onse shown in the argument-matching list. However, standard
#' transformations supported by \code{stats::make.link} are also supported.
#' In the following discussion of ones requiring parameters, 
#' we use \eqn{\alpha} and \eqn{\beta} to
#' denote \code{alpha} and \code{beta}, and \eqn{y} to denote the response variable.
#' The \code{type} argument specifies the following transformations:
#' \describe{
#' \item{\code{"genlog"}}{Generalized logarithmic transformation: \eqn{\log_\beta(y +
#'   \alpha)}, where \eqn{y > -\alpha}.
#'   When \eqn{\beta = 0} (the default), we use \eqn{\log_e(y + \alpha)}}
#' \item{\code{"power"}}{Power transformation: \eqn{(y-\beta)^\alpha}, where \eqn{y > \beta}.
#'   When \eqn{\alpha = 0}, \eqn{\log(y-\beta)} is used instead.}
#' \item{\code{"boxcox"}}{The Box-Cox transformation (unscaled by the geometric
#'   mean): \eqn{((y - \beta)^\alpha - 1) / \alpha}, where \eqn{y > \beta}. 
#'   When \eqn{\alpha = 0}, \eqn{\log(y - \beta)}
#'   is used.}
#' \item{\code{"sympower"}}{A symmetrized power transformation on the whole real
#'   line:
#'   \eqn{|y - \beta|^\alpha\cdot sign(y - \beta)}. There are no restrictions on \eqn{y}, but we
#'   require \eqn{\alpha > 0} in order for the transformation to be monotone and
#'   continuous.}
#' \item{\code{"asin.sqrt"}}{Arcsin-square-root transformation:
#'   \eqn{\sin^{-1}(y/\alpha)^{1/2}}. Typically, \code{alpha} will be either 1 (default) or 100.}
#' \item{\code{"atanh"}}{Arctanh transformation:
#'   \eqn{\tanh^{-1}(y/\alpha)}. Typically, \code{alpha} will be either 1 (default) or 100.}
#' \item{\code{"bcnPower"}}{Box-Cox with negatives allowed, as described for the 
#'   \code{bcnPower} function in the \pkg{car} package. It is defined as the Box-Cox
#'   transformation \eqn{(z^\alpha - 1) / \alpha} of the variable \eqn{z = y + (y^2+\beta^2)^{1/2}}. 
#'   Note that this requires both parameters and that \code{beta > 0}.}
#' \item{\code{"scale"}}{This one is a little different than the others, in that
#'   \code{alpha} and \code{beta} are ignored; instead, they are determined by calling 
#'   \code{scale(y, ...)}. The user should give as \code{y} the response variable in the
#'   model to be fitted to its scaled version.}
#' 
#' }
#' 
#' Note that with the \code{"power"}, \code{"boxcox"}, or \code{"sympower"} transformations, 
#' the argument \code{beta} specifies a location shift. 
#' In the \code{"genpower"} transformation, \code{beta} specifies
#' the base of the logarithm -- however, quirkily, the default of \code{beta = 0}
#' is taken to be the natural logarithm. For example,
#' \code{make.tran(0.5, 10)} sets up the \eqn{\log_{10}(y + \frac12)}
#' transformation. In the \code{"bcnPower"} transformation, \code{beta}
#' must be specified as a positive value.
#' 
#' For purposes of back-transformation, the \samp{sqrt(y) + sqrt(y+1)}
#' transformation is treated exactly the same way as \samp{2*sqrt(y)}, because
#' both are regarded as estimates of \eqn{2\sqrt\mu}.
#' 
#' @section Cases where \code{make.tran} may not be needed:
#' For standard transformations with no parameters, we usually don't need to use
#' \code{make.tran}; just the name of the transformation is all that is needed.
#' The functions \code{\link{emmeans}}, \code{\link{ref_grid}}, and related ones
#' automatically detect response transformations that are recognized by
#' examining the model formula. These are \code{log}, \code{log2}, \code{log10},
#' \code{log1p},
#' \code{sqrt}, \code{logit}, \code{probit}, \code{cauchit}, \code{cloglog}; as
#' well as (for a response variable \code{y}) \code{asin(sqrt(y))},
#' \code{asinh(sqrt(y))}, \code{atanh(y)}, and \code{sqrt(y) + sqrt(y+1)}. 
#' In addition, any
#' constant multiple of these (e.g., \code{2*sqrt(y)}) is auto-detected and
#' appropriately scaled (see also the \code{tran.mult} argument in
#' \code{\link{update.emmGrid}}).
#' 
#' A few additional transformations may be specified as character strings and
#' are auto-detected: \code{"identity"}, \code{"1/mu^2"},
#' \code{"inverse"}, \code{"reciprocal"}, \code{"log10"}, \code{"log2"},
#' \code{"asin.sqrt"}, \code{"asinh.sqrt"}, and \code{"atanh"}.
#' 
#'
#' @section Compound transformations:
#' A transformation that is a function of another function can be created by
#' specifying \code{inner} for the other function. For example, the
#' transformation \eqn{1/\sqrt{y}} can be created either by
#' \code{make.tran("inverse", inner = "sqrt")} or by \code{make.tran("power",
#' -0.5)}. In principle, transformations can be compounded to any depth.
#' Also, if \code{type} is \code{"scale"}, \code{y} is replaced by 
#' \code{inner$linkfun(y)}, because that will be the variable that is scaled.
#' 
#' @return A \code{list} having at least the same elements as those returned by
#'   \code{\link{make.link}}. The \code{linkfun} component is the transformation
#'   itself. Each of the functions is associated with an environment where any 
#'   parameter values are defined.
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
#' mod <- lm(scale(yield[, 1]) ~ Variety, data = MOats)
#' emmeans(mod, "Variety", type = "response")
#' 
#' # Case where scaling is not auto-detected -- and what to do about it:
#' copt <- options(contrasts = c("contr.sum", "contr.poly"))
#' mod.aov <- aov(scale(yield[, 1]) ~ Variety + Error(Block), data = MOats)
#' emm.aov <- suppressWarnings(emmeans(mod.aov, "Variety", type = "response"))
#' 
#' # Scaling was not retrieved, but we can do:
#' emm.aov <- update(emm.aov, tran = make.tran("scale", y = MOats$yield[, 1]))
#' emmeans(emm.aov, "Variety", type = "response")
#' 
#' ### Compound transformations
#' # The following amount to the same thing:
#' t1 <- make.tran("inverse", inner = "sqrt")
#' t2 <- make.tran("power", -0.5)
#' 
#' options(copt)
#'
#' 
#' \dontrun{
#' ### An existing model 'mod' was fitted with a y^(2/3) transformation...
#'   ptran = make.tran("power", 2/3)
#'   emmeans(mod, "treatment", tran = ptran)
#' }
make.tran = function(type = c("genlog", "power", "boxcox", "sympower", 
                              "asin.sqrt", "atanh", "bcnPower", "scale"), 
                              alpha = 1, beta = 0, param, y, inner, ...) {
    if(!missing(inner)){
        if(type == "scale") {
            if(is.character(inner))
                inner = make.tran(inner)
            y = inner$linkfun(y)
        }
        return(.make.compound.link(make.tran(type, alpha, beta, param, y), inner))
    }
    
    txt = type
    type = try(match.arg(type), silent = TRUE)
    if(inherits(type, "try-error"))
        return(.make.link(txt))
    remove(list = "txt")
    
    # backward compat
    if(length(alpha) > 1) param = alpha  # unnamed parem
    if (!missing(param)) {
        if (length(param) > 1) beta = param[2]
        alpha = param[1]
     }
    mu.lbl = "mu"
    if (beta != 0)
        mu.lbl = paste0("(mu - ", round(beta, 3), ")")
    if(type == "scale") {
        sy = scale(y, ...)
        if(is.null(alpha <- attr(sy, "scaled:center")))
            alpha = 0
        if(is.null(beta <- attr(sy, "scaled:scale")))
            beta = 1
        remove(list = c("y", "sy")) # remove baggage from env
    }
    switch(type,
           genlog = { # beta serves in the role of the base of the log
               if((beta < 0) || (beta == 1))
                   stop('"genlog" transformation must have a positive base != 1')
               logbase = ifelse(beta == 0, 1, log(beta))
               xlab = ifelse(beta == 0, "", paste0(" (base ", round(beta, 3), ")"))
               list(linkfun = function(mu) log(pmax(mu + alpha, 0)) / logbase,
                    linkinv = function(eta) pmax(exp(logbase * eta), .Machine$double.eps) - alpha,
                    mu.eta = function(eta) logbase * pmax(exp(logbase * eta), .Machine$double.eps),
                    valideta = function(eta) TRUE,
                    alpha = alpha, logbase = logbase,
                    name = paste0("log(mu + ", round(alpha,3), ")", xlab)
               )
           },
           power = {
               if (alpha == 0) {
                   if(beta == 0) make.link("log")
                   else make.tran("genlog", -beta)
               }
               else list(
                   linkfun = function(mu) pmax(mu - beta, 0)^alpha,
                   linkinv = function(eta) beta + pmax(eta, 0)^(1/alpha),
                   mu.eta = function(eta) pmax(eta, 0)^(1/alpha - 1) / alpha,
                   valideta = function(eta) all(eta > 0),
                   alpha = alpha, beta = beta, 
                   name = ifelse(alpha > 0, 
                                 paste0(mu.lbl, "^", round(alpha,3)),
                                 paste0(mu.lbl, "^(", round(alpha,3), ")"))
               )
           },
           boxcox = {
               if (alpha == 0) {
                   result = if(beta == 0) make.link("log")
                   else make.tran("genlog", -beta)
                   return (result)
               }
               min.eta = ifelse(alpha > 0, -1 / alpha, -Inf)
               xlab = ifelse(beta == 0, "", paste0(" of (y - ", round(beta, 3), ")"))
               list(
                   linkfun = function(mu) ((mu - beta)^alpha - 1) / alpha,
                   linkinv = function(eta) beta + (1 + alpha * pmax(eta, min.eta))^(1/alpha),
                   mu.eta = function(eta) (1 + alpha * pmax(eta, min.eta))^(1/alpha - 1),
                   valideta = function(eta) all(eta > min.eta),
                   alpha = alpha, beta = beta, 
                   name = paste0("Box-Cox (lambda = ", round(alpha, 3), ")", xlab)
               )
           },
           sympower = {
               if (alpha <= 0) 
                   stop('"sympower" transformation requires positive alpha')
               absmu.lbl = gsub("\\(|\\)", "|", mu.lbl)
               if (beta == 0) mu.lbl = paste0("(", mu.lbl, ")")
               list(linkfun = function(mu) sign(mu - beta) * abs(mu - beta)^alpha,
                    linkinv = function(eta) beta + sign(eta) * abs(eta)^(1/alpha),
                    mu.eta = function(eta) (abs(eta))^(1/alpha - 1),
                    valideta = function(eta) all(eta > min.eta),
                    alpha = alpha, beta = beta, 
                    name = paste0(absmu.lbl, "^", round(alpha,3), " * sign", mu.lbl)
               )
           },
           asin.sqrt = {
               mu.lbl = ifelse(alpha == 1, "mu", paste0("mu/", round(alpha,3)))
               list(linkfun = function(mu) asin(sqrt(mu/alpha)),
                    linkinv = function(eta) alpha * sin(pmax(pmin(eta, pi/2), 0))^2,
                    mu.eta = function(eta) alpha * sin(2*pmax(pmin(eta, pi/2), 0)),
                    valideta = function(eta) all(eta <= pi/2) && all(eta >= 0),
                    alpha = alpha,
                    name = paste0("asin(sqrt(", mu.lbl, "))")
               )
           },
           atanh = {
               mu.lbl = ifelse(alpha == 1, "mu", paste0("mu/", round(alpha,3)))
               list(linkfun = function(mu) atanh(mu/alpha),
                    linkinv = function(eta) alpha * tanh(eta),
                    mu.eta = function(eta) alpha * (1 - tanh^2(eta)),
                    valideta = function (eta) all(is.finite(eta)) && all(eta > -1) && all(eta < 1),
                    alpha = alpha,
                    name = paste0("atanh(", mu.lbl, ")")
               )
           },
           bcnPower = {
               if(beta <= 0)
                   stop ("The value of 'beta' must be strictly positive.")
               list(
                   linkfun = function(mu) {
                       s = sqrt(mu^2 + beta^2)
                       if (abs(alpha) < 1e-10) log(.5*(mu + s))
                       else ((0.5 * (mu + s))^alpha - 1) / alpha  },
                   linkinv = function(eta) {
                       q = if (abs(alpha) < 1e-10) 2 * exp(eta)
                           else 2 * (alpha * eta + 1) ^ (1/alpha)
                       (q^2 - beta^2) / (2 * q) },
                   mu.eta = function(eta) {
                       if (abs(alpha) < 1e-10) { q = 2 * exp(eta); dq = q }
                       else { q = 2 * (alpha * eta + 1) ^ (1/alpha)
                           dq = 2 * (alpha * eta + 1)^(1/alpha - 1) }
                       0.5 * (1 + (beta/q)^2) * dq },
                   valideta = function(eta) all(eta > 0),
                   alpha = alpha, beta = beta, 
                   name = paste0("bcnPower(", signif(alpha,3), ", ", signif(beta,3), ")")
               )
           },
           scale = list(
               linkfun = function(mu) (mu - alpha) / beta,
               linkinv = function(eta) beta * eta + alpha,
               mu.eta = function(eta) rep(beta, length(eta)),
               valideta = function(eta) TRUE,
               alpha = alpha, beta = beta, 
               name = paste0("scale(", signif(alpha, 3), ", ", signif(beta, 3), ")")
           )
    )
}

# Create a compound link function f(g(y)) where outer and inner are make.link results
.make.compound.link = function(outer, inner, ...) {
    if(is.character(inner))
        inner = make.tran(inner)
    name = outer$name
    if(length(grep("mu", name)) == 0) name = paste0(name, "(mu)")
    name = gsub("mu", inner$name, name)
    list(
        linkfun = function(mu) outer$linkfun(inner$linkfun(mu)),
        linkinv = function(eta) inner$linkinv(outer$linkinv(eta)),
        mu.eta = function(eta) inner$mu.eta(outer$linkinv(eta))*outer$mu.eta(eta),
        valideta = function(eta) inner$valideta(outer$linkinv(eta)),
        name = name
    )
}


#' @rdname make.tran
#' @export
#' @return \code{inverse} returns the reciprocal of its argument. It allows
#'   the \code{"inverse"} link to be auto-detected as a response transformation.
#' @examples
#' 
#' pigs.lm <- lm(inverse(conc) ~ source + factor(percent), data = pigs)
#' emmeans(pigs.lm, "source", type = "response")
inverse = function(y) 1/y


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
             linkfun = log10,
             linkinv = function(eta) 10^eta,
             mu.eta = function(eta) 10^eta * log(10),
             name = "log10"
         ),
         log2 = list(
             linkfun = log2,
             linkinv = function(eta) 2^eta,
             mu.eta = function(eta) 2^eta * log(2),
             name = "log2"
         ),
         log1p = list(
             linkfun = log1p,
             linkinv = expm1,
             mu.eta = exp,
             name = "log1p"
         ),
         asin.sqrt = make.tran("asin.sqrt"),
         `asin.sqrt./` = make.tran("asin.sqrt", 100),
         asinh.sqrt = list(
             linkinv = function(eta) sinh(eta)^2,
             mu.eta = function(eta) sinh(2 * eta),
             name = "asinh(sqrt(mu))"
         ),
         atanh = make.tran("atanh"),
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


### Internal routine to make a scales::trans object
.make.scale = function(misc) {
    if (!requireNamespace("scales", quietly = TRUE)) 
        stop("type = \"scale\" requires the 'scales' package to be installed")
    tran = misc$tran
    if (is.null(tran))
        return(scales::identity_trans())
    if (is.character(tran)) {
        # is it a canned scale?
        if ((length(intersect(names(misc), c("tran.mult", "tran.offset"))) == 0) && 
            tran %in% c("log", "log1p", "log2", "log10", "sqrt", "logit", "probit", 
                        "exp", "identity"))
            return(get(paste(tran, "trans", sep = "_"), envir = asNamespace("scales"))())
        # not built-in, so let's get a list
        tran = .make.link(tran)
    }
    # tran is a list. we'll incorporate any scaling
    tran$mult = ifelse(is.null(misc$tran.mult), 1, misc$tran.mult)
    tran$offset = ifelse(is.null(misc$tran.offset), 0, misc$tran.offset)
    with(tran, 
         scales::trans_new(name, 
                           function(x) mult * linkfun(x + offset), 
                           function(z) linkinv(z / mult) - offset))
}

