##############################################################################
#    Copyright (c) 2018 Russell V. Lenth                                     #
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

### Quick-and-dirty support for otherwise unsupported models

#' Quick and dirty reference grid
#' 
#' This function may make it possible to compute a reference grid for a model 
#' object that is otherwise not supported.
#' 
#' If \code{object} is specified, it is used to try to obtain certain 
#' other arguments, as detailed below. The user should ensure that these defaults
#' will work. The default values for the arguments are as follows:
#' \itemize{
#'   \item{\code{formula}: }{Required unless obtainable via \code{formula(object)}}
#'   \item{\code{data}: }{Required if variables are not in \code{parent.frame()} or 
#'       obtainable via \code{object$data}}
#'   \item{\code{coef}: }{\code{coef(object)}}
#'   \item{\code{mcmc}: }{\code{object$sample}}
#'   \item{\code{vcov}: }{\code{vcov(object)}}
#'   \item{\code{df}: }{Set to \code{Inf} if not available in \code{object$df.residual}}
#'   \item{\code{subset}: }{\code{NULL} (so that all observations in \code{data} are used)}
#'   \item{\code{contrasts}: }{\code{NULL} (so that \code{getOption("contrasts")} is used)}
#' }
#' 
#' The functions \code{\link{qdrg}} and \code{emmobj} are close cousins, in that
#' they both produce \code{emmGrid} objects. When starting with summary
#' statistics for an existing grid, \code{emmobj} is more useful, while
#' \code{qdrg} is more useful when starting from a fitted model.
#'
#' @param formula Formula for the fixed effects
#' @param data Dataset containing the variables in the model
#' @param coef Fixed-effect regression coefficients (must conform to formula)
#' @param mcmc Posterior sample of fixed-effect coefficients
#' @param vcov Variance-covariance matrix of the fixed effects
#' @param object Optional model object. If provided, it is used to set 
#'     certain other arguments, if not specified. See Details.
#' @param df Error degrees of freedom
#' @param subset Subset of \code{data} used in fitting the model
#' @param weights Weights used in fitting the model
#' @param contrasts List of contrasts specified in fitting the model
#' @param link Link function (character or list) used, if a generalized linear model.
#'     (Note: response transformations are auto-detected from \code{formula})
#' @param qr QR decomposition of the model matrix; needed only if there are \code{NA}s
#'     in \code{coef}.
#' @param ... Optional arguments passed to \code{\link{ref_grid}}
#'
#' @return An \code{emmGrid} object constructed from the arguments
#' 
#' @seealso \code{\link{emmobj}} for an alternative way to construct an \code{emmGrid}.
#' 
#' @export
#' @examples
#' if (require(biglm)) {
#'   # Post hoc analysis of a "biglm" object -- not supported by emmeans
#'   bigmod <- biglm(log(conc) ~ source + factor(percent), data = pigs)
#'    
#'   rg2 <- qdrg(object = bigmod, data = pigs)
#'   summary(emmeans(rg2, "source"), type = "response")
#' }
#' if(require(coda) && require(lme4)) {
#'   # Use a stored example having a posterior sample
#'   # Model is based on the data in lme4::cbpp
#'   
#'   post <- readRDS(system.file("extdata", "cbpplist", package = "emmeans"))$post.beta
#'   rg1 <- qdrg(~ size + period, data = lme4::cbpp, mcmc = post, link = "logit")
#'   summary(rg1, type = "response")
#' }
#'
qdrg = function(formula, data, coef, mcmc, vcov, object,
                df, subset, weights, contrasts, link, qr, ...) {
    result = match.call(expand.dots = FALSE)
    if (!missing(object)) {
        if (missing(formula)) 
            result$formula = stats::formula(object)
        if (missing(data)) {
            data = object$data
            if (is.null(data)) data = parent.frame()
            result$data = data
        }
        if (missing(coef)) 
            result$coef = stats::coef(object)
        if (missing(mcmc)) 
            result$mcmc = object$sample
        if (missing(vcov)) 
            result$vcov = stats::vcov(object)
        if(missing(df))
            result$df = object$df.residual
        if(any(is.na(result$coef)))
            result$qr = object$qr
    }
    else {
        if(missing(formula))
            stop("When 'object' is missing, must at least provide 'formula'")
        result$formula = formula
        if(missing(data))
            result$data = parent.frame()
        else
            result$data = data
        if (!missing(coef)) result$coef = coef
        if (!missing(vcov)) result$vcov = vcov
        if(!missing(df)) result$df = df
    }
    if(is.null(result$df))
        result$df = Inf
    if(!missing(mcmc)) result$mcmc = mcmc
    if(!missing(subset)) result$subset = subset
    if(!missing(weights)) result$weights = weights
    if(!missing(contrasts)) result$contrasts = contrasts
    if(!missing(link)) result$link = link
    if(!missing(qr) && any(is.na(result$coef))) result$qr = qr
    
    # make sure "formula" is 2nd element so that response transformation can be found
    fpos = grep("formula", names(result))[1]
    result = result[c(1, fpos, seq_along(result)[-c(1, fpos)])]

    class(result) = c("qdrg", "call")
    ref_grid(result, ...)
}

recover_data.qdrg = function(object, ...) {
    recover_data.call(object, delete.response(terms(object$formula)), object$na.action, ...)
}

vcov.qdrg = function(object, ...)
    object$vcov

emm_basis.qdrg = function(object, trms, xlev, grid, ...) {
    bhat = object$coef
    m = suppressWarnings(model.frame(trms, grid, na.action = na.pass, xlev = xlev))
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    bhat = as.numeric(bhat) 
    V = .my.vcov(object, ...)
    if (!is.null(object$mcmc)) {
        if (is.null(object$coef)) bhat = apply(object$mcmc, 2, mean)
        if (is.null(object$vcov)) V = cov(object$mcmc)
    }
    
    nbasis = estimability::all.estble
    if (sum(is.na(bhat)) > 0) {
        if(!is.na(object$qr))
            nbasis = estimability::nonest.basis(object$qr)
        else
            warning("Non-estimable cases can't be determined.\n",
                "To rectify, provide appropriate 'qr' in call to qdrg()")
    }
    
    misc = list()
    if (!is.null(object$link)) {
        misc = .std.link.labels(eval(list(link = object$link)), misc)
        dffun = function(k, dfargs) Inf
        dfargs = list()
    }
    else {
        dfargs = list(df = object$df)
        dffun = function(k, dfargs) dfargs$df
    }
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, 
         misc=misc, post.beta=object$mcmc)
    }
    