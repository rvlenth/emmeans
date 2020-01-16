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
            delta = .diag(lnk$mu.eta(lp)) %*% X
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
            delta = .diag(mu1) %*% cbind(.diag(1 - mu2) %*% X1, .diag(-mu2prime) %*% X2)
            bhat = (1 - mu2) * mu1
        }
        else { # mode = "prob0"
            p0 = 1 - .prob.gt.0(object$dist, mu1, object$theta)
            dp0 = - .dprob.gt.0(object$dist, mu1, object$theta, "log", lp1)
            bhat = (1 - mu2) * p0 + mu2
            delta = cbind(.diag((1 - mu2) * dp0) %*% X1, .diag(mu2prime * (1 - p0)) %*% X2)
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
            if (mode != "prob0") {
               delta = .diag(lnk$mu.eta(lp)) %*% X
            }
            else {
                bhat = 1 - .prob.gt.0(object$dist$zero, bhat, object$theta["zero"])
                db = - .dprob.gt.0(object$dist$zero, bhat, object$theta["zero"], misc$tran, lp)
                delta = .diag(db) %*% X
            }
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
        mu1 = as.numeric(exp(X1 %*% b1 + off1))
        theta1 = object$theta["count"]
        p1 = .prob.gt.0(object$dist$count, mu1, theta1)
        dp1 = .dprob.gt.0(object$dist$count, mu1, theta1, "", 0) # binomial won't happen

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
        p2 = .prob.gt.0(object$dist$zero, mu2, theta2)
        dp2 = .dprob.gt.0(object$dist$zero, mu2, theta2, object$link, lp2)

        if (mode == "response") {
            bhat = p2 * mu1 / p1
            delta = cbind(.diag(bhat*(1 - mu1 * dp1 / p1)) %*% X1,
                          .diag(mu1 * dp2 / p1) %*% X2)
        }
        else {  ## mode == "zero"
            bhat = p2 / p1
            delta = cbind(.diag(-p2 * dp1 / p1^2) %*% X1,
                          .diag(dp2 / p1) %*% X2)
        }
        V = delta %*% tcrossprod(.pscl.vcov(object, model = "full", ...), delta)
        X = .diag(1, length(bhat))
        
        misc = list(estName = mode, offset.mult = 0)
    }
    nbasis = estimability::all.estble
    dffun = function(k, dfargs) object$df.residual
    dfargs = list()
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, 
         dffun = dffun, dfargs = dfargs, misc = misc)
}

# utility for prob (Y > 0 | dist, mu, theta)
.prob.gt.0 = function(dist, mu, theta) {
    switch(dist,
       binomial = mu,
       poisson = 1 - exp(-mu),
       negbin = 1 - (theta / (mu + theta))^theta,
       geometric = 1 - 1 / (1 + mu)
    )
}

# utility for d/d(eta) prob (Y > 0 | dist, mu, theta)
.dprob.gt.0 = function(dist, mu, theta, link, lp) {
    switch(dist,
       binomial = stats::make.link(link)$mu.eta(lp),
       poisson = mu * exp(-mu),
       negbin = mu * (theta /(mu + theta))^(1 + theta),
       geometric = mu / (1 + mu)^2  
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

# Explanatory notes for hurdle models
# -----------------------------------
#     We have a linear predictor eta = X%*%beta + offset
#     mu = h(eta) where h is inverse link (usually exp but not always)
#     Define p = P(Y > 0 | mu). This comes out to...
#         binomial: mu
#         poisson: 1 - exp(-mu)
#         negbin: 1 - (theta/(mu+theta))^theta
#         geometric: 1 - 1/(mu+1)
#     Define dp = dp/d(eta). Note - when h(mu)=exp(mu) we have dp = mu*dp/d(mu)
#         binomial: h'(eta)
#         poisson: mu*exp(-mu)
#         negbin: mu*(theta/(mu+theta))^(theta+1)
#         geometric: mu/(mu+1)^2
#     
#     This gives us what we need to find the estimates and apply the delta method
#     In the code we index these notations with 1 (count model) and 2 (zero model)
#     And we treat theta1 and theta2 as constants
#
#!!! In theory, above seems correct, and estimates match those from predict.hurdle.
#!!! But SEs don't seem right. 
#!!! They do seem right though if I omit the factor of mu in dp
#!!! when link is log


### Simulation-based approach for ZI and hurdle models
### This returns a list of ther same form as an emm_basis() method
##   X is model matrix for both models (X.count | X.zi)
##   ncoef.c # cols in X.count
##   bhat is all regression coefs -- OR a matrix w/ posterior samples
##   V is combined vcov matrix (ignored if bhat is a matrix)
##   links is vector (or list) of links for the 2 parts of the model
##   fams is list of family names
##   hurdle is a named list(famc, thetac, famz, thetaz) or empty if ZI
##   N.sim is # simulations desired
##   keep.sim set to TRUE to return sims in post.beta
##   df is d.f. to return
## NOTE if bhat is a matrix, keep.sim is set to TRUE and N.sim is ignored
#' @export
.zi.simulate = function(X, ncoef.c, bhat, V, links, hurdle = list(), 
                        N.sim = 1000, keep.sim = FALSE,
                        df = Inf, misc = list()) 
{
    if(is.matrix(bhat)) {
        keep.sim = TRUE
        B = bhat
        bhat = NA
    }
    else {
        if (length(bhat) != ncol(V))
            stop("Non-estimable cases not yet supported for zero-inflation calculations")
        B = mvtnorm::rmvnorm(N.sim, bhat, V)
    }
    back.tran = function(idx, link) {
        W = B[, idx] %*% t(X[, idx])
        if (is.character(link))
            link = make.link(link)
        link$linkinv(W)
    }
    if(!is.na(bhat[1]))
        B = rbind(bhat, B) # put our pt est in 1st row
    ic = seq_len(ncoef.c)
    iz = setdiff(seq_len(ncol(B)), ic)
    C = back.tran(ic, links[[1]])
    Z = back.tran(iz, links[[2]])
    if (length(hurdle) >= 4)
        R = C * .prob.gt.0(hurdle$famz, Z, hurdle$thetaz) / 
                .prob.gt.0(hurdle$famc, C, hurdle$thetac)
    else
        R = (1 - Z) * C
    if(!is.na(bhat[1])){
        bhat = R[1, ]
        R = R[-1, ]
    }
    else
        bhat = apply(R, 2, mean)
    V = cov(R)
    dffun = function(k, dfargs) dfargs$df
    dfargs = list(df = df)
    if(is.null(misc$est.name)) 
        misc$estName = "emmean"
    if(!keep.sim)
        post.beta = matrix(NA)
    list(X = diag(length(bhat)), bhat = bhat, nbasis = estimability::all.estble, V = V, 
         dffun = dffun, dfargs = dfargs, misc = misc, post.beta = R)
}
    
