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

# Support for zeroinfl and hurdle models (pscl)

# We'll support two optional arguments:
# mode     -- type of result required
# lin.pred -- TRUE:  keep linear predictor and link
#             FALSE: back-transform (default)
#
# With lin.pred = FALSE and mode %in% c("response", "count", "zero"), we
# will return comparable results to predict(..., type = mode)
# with mode = "prob0", same results as predict(..., type = "prob")[, 1]
#
# lin.pred only affects results for mode %in% c("count", "zero"). 
# When lin.pred = TRUE, we get the actual linear predictor and link function
# for that part of the model.



# ----- zeroinfl objects -----

#' @exportS3Method recover_data zeroinfl
recover_data.zeroinfl = function(object, mode = c("response", "count", "zero", "prob0"), ...) {
    fcall = object$call
    mode = match.arg(mode)
    if (mode %in% c("count", "zero"))
        trms = delete.response(terms(object, model = mode))
    else ### mode = %in% c("response", "prob0")
        trms = delete.response(object$terms$full)
    # Make sure there's an offset function available
    env = new.env(parent = attr(trms, ".Environment"))
    env$offset = function(x) x
    attr(trms, ".Environment") = env
    recover_data(fcall, trms, object$na.action, ...)
}


#' @exportS3Method emm_basis zeroinfl     
emm_basis.zeroinfl = function(object, trms, xlev, grid, 
        mode = c("response", "count", "zero", "prob0"), lin.pred = FALSE, ...) 
{
    mode = match.arg(mode)
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    if (mode %in% c("count", "zero")) {
        contr = object$contrasts[[mode]]
        X = model.matrix(trms, m, contrasts.arg = contr)
        bhat = coef(object, model = mode)
        V = .pscl.vcov(object, model = mode, ...)
        if (mode == "count")
            misc = list(tran = "log", inv.lbl = "count")
        else
            misc = list(tran = object$link, inv.lbl = "prob")
        if (!lin.pred) { # back-transform the results
            lp = as.numeric(X %*% bhat + .get.offset(trms, grid))
            lnk = make.link(misc$tran)
            bhat = lnk$linkinv(lp)
#            delta = .diag(lnk$mu.eta(lp)) %*% X
            delta = sweep(X, 1, lnk$mu.eta(lp), "*")
            V = delta %*% tcrossprod(V, delta)
            X = diag(1, length(bhat))
            misc = list(offset.mult = 0)
        }
    }
    else { ## "response", "prob0"
        trms1 = delete.response(terms(object, model = "count"))
        off1 = .get.offset(trms1, grid)
        contr1 = object$contrasts[["count"]]
        X1 = model.matrix(trms1, m, contrasts.arg = contr1)
        b1 = coef(object, model = "count")
        lp1 = as.numeric(X1 %*% b1 + off1)
        mu1 = exp(lp1)
        
        trms2 = delete.response(terms(object, model = "zero"))
        off2 = .get.offset(trms2, grid)
        contr2 = object$contrasts[["zero"]]
        X2 = model.matrix(trms2, m, contrasts.arg = contr2)
        b2 = coef(object, model = "zero")
        lp2 = as.numeric(X2 %*% b2) + off2
        mu2 = object$linkinv(lp2)
        mu2prime = stats::make.link(object$link)$mu.eta(lp2)
        
        if(mode == "response") {
            delta = cbind(sweep(X1, 1, mu1 * (1 - mu2), "*"), 
                          sweep(X2, 1, -mu1 * mu2prime, "*"))
            bhat = (1 - mu2) * mu1
        }
        else { # mode = "prob0"
            tmp = .zi.support(mu1, object$theta, .make.p0(object$dist))
            p0 = tmp[1, ]; dp0 = tmp[2, ]
            bhat = (1 - mu2) * p0 + mu2
            delta = cbind(sweep(X1, 1, (1 - mu2) * dp0, "*"),
                          sweep(X2, 1, mu2prime * (1 - p0), "*"))
        }
        V = delta %*% tcrossprod(.pscl.vcov(object, model = "full", ...), delta)
        X = diag(1, length(bhat))
        misc = list(offset.mult = 0)
    }
    nbasis = estimability::all.estble
    dffun = function(k, dfargs) Inf
    dfargs = list()
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, 
         dffun = dffun, dfargs = dfargs, misc = misc)
}



#### Support for hurdle models

#' @exportS3Method recover_data hurdle
recover_data.hurdle = function(object, mode = c("response", "count", "zero", "prob0"), ...) {
    fcall = object$call
    mode = match.arg(mode)
    if (mode %in% c("count", "zero"))
        trms = delete.response(terms(object, model = mode))
    else ### mode = "mean" or "prob.ratio"
        trms = delete.response(object$terms$full)
    # Make sure there's an offset function available
    env = new.env(parent = attr(trms, ".Environment"))
    env$offset = function(x) x
    attr(trms, ".Environment") = env
    recover_data(fcall, trms, object$na.action, ...)
}

# see expl notes afterward for notations in some of this
#' @exportS3Method emm_basis hurdle       
emm_basis.hurdle = function(object, trms, xlev, grid, 
                            mode = c("response", "count", "zero", "prob0"), 
                            lin.pred = FALSE, ...) 
{
    mode = match.arg(mode)
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    if ((lin.pred && mode %in% c("count", "zero")) || (!lin.pred && mode %in% c("count", "prob0"))) {
        model = ifelse(mode == "count", "count", "zero")
        contr = object$contrasts[[model]]
        X = model.matrix(trms, m, contrasts.arg = contr)
        bhat = coef(object, model = model)
        V = .pscl.vcov(object, model = model, ...)
        misc = switch(object$dist[[model]],
                        binomial = list(tran = object$link, inv.lbl = "prob"),
                        list(tran = "log", inv.lbl = "count"))
        if (!lin.pred) { # back-transform
            lp = as.numeric(X %*% bhat + .get.offset(trms, grid))
            lnk = make.link(misc$tran)
            bhat = lnk$linkinv(lp)
            mult = lnk$mu.eta(lp)
            if(mode == "prob0") {
                shape = 
                tmp = .zi.support(bhat, object$theta["zero"], .make.p0(object$dist$zero))
                bhat = tmp[1, ]
                mult = mult * tmp[2, ]
            }
            delta = sweep(X, 1, mult, "*")
            V = delta %*% tcrossprod(V, delta)
            X = diag(1, length(bhat))
            misc = list(offset.mult = 0)
        }
    }
    else {   ### "zero" or "response" with implied lin.pred = FALSE
        trms1 = delete.response(terms(object, model = "count"))
        off1 = .get.offset(trms1, grid)
        contr1 = object$contrasts[["count"]]
        X1 = model.matrix(trms1, m, contrasts.arg = contr1)
        b1 = coef(object, model = "count")
        mu1 = mu1prime = as.numeric(exp(X1 %*% b1 + off1))
        theta1 = object$theta["count"]

        trms2 = delete.response(terms(object, model = "zero"))
        off2 = .get.offset(trms2, grid)
        contr2 = object$contrasts[["zero"]]
        X2 = model.matrix(trms2, m, contrasts.arg = contr2)
        b2 = coef(object, model = "zero")
        lp2 = as.numeric(X2 %*% b2 + off2)
        mu2 = switch(object$dist$zero, 
                     binomial = object$linkinv(lp2),
                     exp(lp2)  )
        theta2 = object$theta["zero"]
        mu2prime = switch(object$dist$zero, 
                          binomial = stats::make.link("logit")$mu.eta(lp2),
                          exp(lp2) )
        if (mode == "response") {
            tmp = .hurdle.support(mu1, theta1, .make.p0(object$dist$count), \(mu, shape) mu, 
                                  mu2, theta2, .make.p0(object$dist$zero))
            bhat = tmp[1, ]
            delta = cbind(sweep(X1, 1, tmp[2, ] * mu1prime, "*"),
                          sweep(X2, 1, tmp[3, ] * mu2prime, "*"))
        }
        else {  ## mode == "zero"
            tmp1 = .zi.support(mu1, theta1, .make.p0(object$dist$count))
            tmp2 = .zi.support(mu2, theta2, .make.p0(object$dist$zero))
            p1 = 1 - tmp1[1, ]; dp1 = -tmp1[2, ]
            p2 = 1 - tmp2[1, ]; dp2 = -tmp2[2, ]
            bhat = p2 / p1
            # delta = cbind(.diag(-p2 * dp1 / p1^2) %*% X1,
            #               .diag(dp2 / p1) %*% X2)
            delta = cbind(sweep(X1, 1, -p2 * dp1, "*"),
                          sweep(X2, 1, dp2 / p1, "*"))
        }
        V = delta %*% tcrossprod(.pscl.vcov(object, model = "full", ...), delta)
        X = .diag(1, length(bhat))
        
        misc = list(estName = mode, offset.mult = 0)
    }
    nbasis = estimability::all.estble
    dffun = function(k, dfargs) dfargs$df
    dfargs = list(df = object$df.residual)
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, 
         dffun = dffun, dfargs = dfargs, misc = misc)
}


#' @rdname extending-emmeans
#' @order 93
#'
#' @param cmu,zmu In \code{.hurdle.support} and \code{.zi.support}, 
#'   these specify a vector of back-transformed 
#'   estimates for the count and zero model, respectively
#' @param cshape,zshape Shape parameter for the count and zero model, respectively 
#' @param cp0,zp0 Function of \code{(mu, shape)} for computing Prob(Y = 0)
#'        for the count and zero model, respectively
#' @param cmean Function of \code{(mu, shape)} for computing the mean of the
#'        count model. Typically, this just returns \code{mu}
#'
#' @return \code{.hurdle.support} returns a matrix with 3 rows containing the
#'       estimated mean responses and the differentials wrt \code{cmu} and \code{zmu},
#'       resp. 
#' @export
.hurdle.support = function(cmu, cshape, cp0, cmean, 
                           zmu, zshape, zp0) {
    if (is.null(cshape) || is.na(cshape)) cshape = 1
    if (is.null(zshape) || is.na(zshape)) zshape = 1
    mfcn = function(x) (1 - zp0(x[2], zshape)) * cmean(x[1], cshape) / 
        (1 - cp0(x[1], cshape))
    result = sapply(seq_along(cmu), function(i) {
        x = c(cmu[i], zmu[i])
        d = numDeriv::jacobian(mfcn, x)
        c(mfcn(x), d)
    })
    rownames(result) = c("mean", "dcmu", "dzmu")
    result
}

#' @rdname extending-emmeans
#' @order 94
#' @return \code{.zi.support} returns a matrix with 2 rows containing the
#'       estimated probabilities of 0 and the differentials wrt \code{mu}.
#'       See the section on hurdle and zero-inflated models.
#' @section Support for Hurdle and Zero-inflated models:
#'   The functions \code{.hurdle.support} and \code{.zi.support} help facilitate
#'   calculations needed to estimate the mean response (count model and zero model
#'   combined) of these models. \code{.hurdle.support} returns a matrix of three rows.
#'   The first is the estimated mean for a hurdle model, and the 2nd and 3rd rows are
#'   differentials for the count and zero models, which needed for delta-method
#'   calculations. To use these, regard the \code{@linfct} slot as comprising
#'   two sets of columns, for the count and zero models respectively. To do
#'   the delta method calculations, multiply the rows of the count part by its 
#'   differentials times \code{link$mu.eta} evcaluated at that part of the linear predictor.
#'   Do the same for the zero part, using its differentials and \code{mu.eta}.
#'   If the resulting matrix is \bold{A}, then the covariance of the mean response
#'   is \bold{AVA'} where \bold{V}is the \code{@V} slot of the object.
#'   
#'   The function \code{zi.support} works the same way, only it is much simpler,
#'   and is used to estimate the probability of 0 and its differential for either 
#'   part of a zero-inflated model or hurdle model.
#'   
#'   See the code for \code{emm_basis.zeroinfl} and \code{emm_basis.hurdle}
#'   for how these are used with models fitted by the \pkg{pscl} package.
#' @export
#'
.zi.support = function(zmu, zshape, zp0) {
    if (is.null(zshape) || is.na(zshape)) zshape = 1
    pfcn = function(x) zp0(x, zshape)
    result = sapply(seq_along(zmu), function(i) {
        c(pfcn(zmu[i]), numDeriv::jacobian(pfcn, zmu[i]))
    })
    rownames(result) = c("p0", "dmu")
    result
}


.make.p0 = function(dist) {
    switch(dist,
           bernoulli = \(mu, shape) 1 - mu,
           binomial = \(mu, shape) (1 - mu)^shape,
           poisson = \(mu, shape) exp(-mu),
           negbin = \(mu, shape) (shape / (mu + shape))^shape,
           geometric = \(mu, shape) 1 / (1 + mu)
    )
}



# special version of .my.vcov that accepts (and requires!) model argument
.pscl.vcov = function(object, model, vcov. = stats::vcov, ...) {
    if (is.function(vcov.))
        vcov. = vcov.(object, model = model)
    else if (!is.matrix(vcov.))
        stop("vcov. must be a function or a square matrix")
    vcov.
}

