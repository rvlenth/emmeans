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

# Support for objects in the *rms* package

recover_data.rms = function(object, ...) {
    fcall = object$call
    recover_data(fcall, delete.response(terms(object)), object$na.action$omit, ...)
}

# TODO: 
# 1. If multivariate - like mlm method?
# 2. orm cases?

emm_basis.rms = function(object, trms, xlev, grid, 
        mode = c("middle", "latent", "linear.predictor", "cum.prob", "exc.prob", "prob", "mean.class"), 
        vcov., ...) {
    mode = match.arg(mode)
    bhat = coef(object) 
    if (missing(vcov.))
        V = vcov(object, intercepts = "all")
    else
        V = .my.vcov(object, vcov.)
    misc = list()
    
    X = predict(object, newdata = grid, type = "x")
    #xnames = dimnames(X)[[2]]
    #intcpts = setdiff(names(bhat), xnames)
    nint = length(bhat) - ncol(X)
    intcpts = names(bhat)[seq_len(nint)]
    xnames = setdiff(names(bhat), intcpts)
        
    if (length(intcpts) == 1) 
        mode = "single" # stealth mode for ordinary single-intercept case
    if (mode %in% c("single", "middle", "latent")) {
        X = cbind(1, X)
        mididx = ifelse(mode != "middle", 1, as.integer((1 + length(intcpts)) / 2))
        dimnames(X)[[2]][1] = switch(mode,
            single = intcpts,
            middle = intcpts[mididx],
            latent = "avg.intercept")
        if (mode == "middle") {
            nms = c(intcpts[mididx], xnames)
            bhat = bhat[nms]
            V = V[nms, nms, drop = FALSE]
        }
        else if (mode == "latent") {
            bhat = c(mean(bhat[intcpts]), bhat[xnames])
            nx = length(xnames)
            J1 = rbind(rep(1/nint, nint), 
                       matrix(0, nrow = nx, ncol = nint))
            J2 = rbind(0, diag(1, nx))
            J = cbind(J1, J2)
            V = J %*% V %*% t(J)
        }
        ### else mode == "single" and all is OK as it is
    }
    else { # mode %in% c("linear.predictor", "cum.prob", "exc.prob", "prob", "mean.class")
        misc$ylevs = list(cut = intcpts)
        I = diag(1, nint)
        J = matrix(1, nrow = nrow(X))
        JJ = matrix(1, nrow=nint)
        X = cbind(kronecker(I, J), kronecker(JJ, X))
        # Note V is correct as-is
        dimnames(X)[[2]] = c(intcpts, xnames)
        if (mode != "linear.predictor") {
            misc$mode = mode
            misc$postGridHook = .clm.postGrid
            misc$respName = as.character.default(object$terms)[2]
        }
    }
    
    # I think rms does not allow rank deficiency...
    nbasis = estimability::all.estble
    if (!is.null(object$family)) {
        if (!is.character(object$family))
            misc = .std.link.labels(object$family, misc)
        else {
            misc$tran = object$family
            if (misc$tran == "logistic") misc$tran = "logit"
            misc$inv.lbl = switch(class(object)[1], 
                orm = "exc.prob",
                lrm = ifelse(nint == 1, "prob", "exc.prob"),
                "response")
        }
        dffun = function(k, dfargs) Inf
        dfargs = list()
    }
    else {
        dfargs = list(df = object$df.residual)
        if (is.null(dfargs$df)) 
            dfargs$df = Inf
        dffun = function(k, dfargs) dfargs$df
    }
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, 
         dffun=dffun, dfargs=dfargs, misc=misc)
}


