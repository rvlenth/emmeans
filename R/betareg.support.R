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

# Support for 'betareg' class

# mode is same as 'type' in predict.betareg, PLUS 
# mode = "phi.link" refers to link function before back-transforming to "precision"

recover_data.betareg = function(object, mode = c("response", "link", "precision", "phi.link", "variance", "quantile"), ...) {
    fcall = object$call
    mode = match.arg(mode)
    if (mode  %in% c("response", "link"))
        mode = "mean"
    if (mode == "phi.link") 
        mode = "precision"
    if(mode %in% c("mean", "precision"))
        trms = delete.response(terms(object, model = mode))
    else
        trms = delete.response(object$terms$full)
    # Make sure there's an offset function available
    env = new.env(parent = attr(trms, ".Environment"))
    env$offset = function(x) x
    attr(trms, ".Environment") = env
    recover_data(fcall, trms, object$na.action, ...)
}

# PRELIMINARY...
# Currently works correctly only for "resp", "link", "precision", "phi" modes
emm_basis.betareg = function(object, trms, xlev, grid, 
        mode = c("response", "link", "precision", "phi.link", "variance", "quantile"), 
        quantile = .5, ...) {
    mode = match.arg(mode)
#     if (mode %in% c("variance", "quantile"))
#         stop(paste0('"', mode, '" mode is not yet supported.'))
    
    # figure out which parameters we need
    model = if (mode %in% c("response", "link")) "mean"
        else if (mode %in% c("precision", "phi.link")) "precision"
        else "full"
    V = .pscl.vcov(object, model = model) # borrowed from pscl methods
    bhat = coef(object, model = model)
    
    nbasis = estimability::all.estble
    dffun = function(k, dfargs) Inf
    dfargs = list()
    
    
    if (mode %in% c("response", "link", "precision", "phi.link")) {
        m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
        X = model.matrix(trms, m, contrasts.arg = object$contrasts[[model]])
        misc = list(tran = object$link[[model]]$name)
        if (mode %in% c("response", "precision")) {
            misc$postGridHook = ".betareg.pg"
        }
    }
    else { ### (mode %in% c("variance", "quantile"))
        m.trms = delete.response(terms(object, "mean"))
        m.m = model.frame(m.trms, grid, na.action = na.pass, xlev = xlev)
        X = model.matrix(m.trms, m.m, contrasts.arg = object$contrasts$mean)
        m.idx = seq_len(ncol(X))
        m.lp = as.numeric(X %*% bhat[m.idx] + .get.offset(m.trms, grid))
        mu = object$link$mean$linkinv(m.lp)
            
        p.trms = delete.response(terms(object, "precision"))
        p.m = model.frame(m.trms, grid, na.action = na.pass, xlev = xlev)
        Z = model.matrix(p.trms, p.m, contrasts.arg = object$contrasts$precision)
        p.lp = as.numeric(Z %*% bhat[-m.idx] + .get.offset(p.trms, grid))
        phi = object$link$precision$linkinv(p.lp)
        
        if (mode == "variance") {
            bhat = mu * (1 - mu) / (1 + phi)
            dbhat.dm = (1 - 2 * mu) / (1 + phi)
            dbhat.dp = -bhat / (1 + phi)
            delta = cbind(diag(dbhat.dm) %*% X, diag(dbhat.dp) %*% Z)
            V = delta %*% tcrossprod(V, delta)
            misc = list()
        }
        else {  ### (mode = "quantile")
            bhat = as.numeric(sapply(quantile, function(q)
                stats::qbeta(q, phi * mu, phi * (1 - mu))))
            V = matrix(NA, nrow = length(bhat), ncol = length(bhat))
            misc = list(ylevs = list(quantile = quantile))
        }
        X = diag(1, length(bhat))
    }
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}

# Post-grid hook for simple back-transforming
.betareg.pg = function(object, ...) {
    object@misc$postGridHook = NULL
    regrid(object, transform = TRUE)
}


### predict methods
# link: X%*%beta + off_m
# response: mu = h_m(link)
# 
# phi.link: Z%*%gamma + off_p
# precision: phi = h_p(phi.link)
# 
# variance: mu*(1 - mu) / (1 + phi)
# quantile: qbeta(p, mu*phi, (1 - mu)*phi)
#
# Defns:
#   phi = a + b
#   mu = a / (a + b)
# so that phi*mu = a and phi*(1 - mu) = b,
#   Variance = ab / [(a + b)^2 * (a + b + 1)]
