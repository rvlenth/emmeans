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

### Support for 'gam' objects
# Note: This is a mess, because both packages 'gam' and 'mgcv' produce these,
# and they are different. Both inherit from glm and lm, though, so recover_data.lm
# still serves for these (I hope)


# gam::Gam objects...

# We have two args:
#   nboot   # of bootstrap reps to get variances of smooths
emm_basis.Gam = function(object, trms, xlev, grid, nboot = 800, ...) {
    result = emm_basis.lm(object, trms, xlev, grid, ...)
    old.smooth = object$smooth
    if (is.null(old.smooth))  # "just an ordinary glm" (My Fair Lady)
        return(result)
    # else we need to add-in some smoothers
    smooth.frame = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    data = object$smooth.frame
    labs = names(data)
    w = object$weights
    resid = object$residuals
    for (i in seq_along(labs)) {
        lab = labs[i]
        sig = apply(smooth.frame[, i, drop=FALSE], 1, paste, collapse = ":")
        usig = unique(sig)
        rows = lapply(usig, function(s) which(sig == s))
        xeval = smooth.frame[sapply(rows, "[", 1), lab]
        bsel = matrix(0, nrow = length(sig), ncol = length(usig))
        for (j in seq_along(rows))
            bsel[rows[[j]], j] = 1
        
        cl = attr(data[[i]], "call")
        cl$xeval = substitute(xeval)
        z = resid + old.smooth[, lab]
        bh = as.numeric(eval(cl))
        m = length(bh)
        n = length(result$bhat)
        result$bhat = c(result$bhat, bh)
        result$X = cbind(result$X, bsel)
        boot = replicate(nboot, {
                z = sample(resid, replace = TRUE) + old.smooth[, lab]
                as.numeric(eval(cl))
            })
        covar = if(m == 1) var(boot) 
                else       cov(t(boot))
        result$V = rbind(cbind(result$V, matrix(0, nrow = n, ncol = m)),
                         cbind(matrix(0, nrow = m, ncol = n), covar))
    }
    result
}


### emm_basis method for mgcv::gam objects
### extra arg `unconditional` and `freq` as in `vcov.gam`
emm_basis.gam = function(object, trms, xlev, grid,
                         freq = FALSE, unconditional = FALSE,
                         what = 1, ...) {
    # coef() works right for lm but coef.aov tosses out NAs
    bhat = object$coefficients
#    m = suppressWarnings(model.frame(trms, grid, na.action = na.pass, xlev = xlev))
#    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    X = predict(object, newdata = grid, type = "lpmatrix", 
                newdata.guaranteed = TRUE)
    bhat = as.numeric(bhat) 
    # stretches it out if multivariate - see mlm method
    V = .my.vcov(object, freq = freq, unconditional = unconditional, ...)

    sel = attr(X, "lpi")[[what]]

    if (!is.null(sel)) {
        bhat = bhat[sel]
        X = X[, sel]
        V = V[sel, sel]
    }

    # if (sum(is.na(bhat)) > 0)
    #     nbasis = estimability::nonest.basis(object$qr)
    # else
        nbasis = estimability::all.estble

        if (!is.null(sel)) {
            object$family$link = object$family$link[what]

            if (object$family$link == "logb") {
                object$family$link = "log"
            }
        }

        misc = .std.link.labels(object$family, list())
        # dffun = function(k, dfargs) Inf
        # dfargs = list()
#    else {
        dfargs = list(df = object$df.residual)
        dffun = function(k, dfargs) dfargs$df
#    }
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}


### mgcv::gamm objects...
recover_data.gamm = function(object, data = NULL, call = object$gam$call, ...) {
    gam = object$gam
    class(gam) = c("gam", "glm", "lm")
    if (!is.null(data)) {
        gam$call = quote(gamm())
        return(recover_data(gam, data = data, ...))
    }
    else {
        if (is.null(call))
            return("Must supply either 'data' or 'call' with gamm objects")
        gam$call = call
        recover_data(gam, ...)
    }
}

emm_basis.gamm = function(object, ...)
    emm_basis(object$gam, ...)


###===================================================================
# Support for gamlss objects
# 'what' parameter mimics predict.gamlss

recover_data.gamlss = function(object, what = c("mu", "sigma", "nu", "tau"), ...) {
    fcall = object$call
    what = match.arg(what)
    trms = terms(formula(object, what = what))
    recover_data(fcall, delete.response(trms), object$na.action, ...)
}

emm_basis.gamlss = function(object, trms, xlev, grid, 
                            what = c("mu", "sigma", "nu", "tau"), vcov., ...) {
    what = match.arg(what)
    smo.mat = object[[paste0(what, ".s")]]
    if (!is.null(smo.mat))
        stop("gamlss models with smoothing are not yet supported in 'emmeans'",
             call. = NULL)
    
    object$coefficients = object[[paste0(what, ".coefficients")]]
    if (missing(vcov.)) {
        # tedious code to pull needed vcov elements
        # Gotta do this before messing up the object
        V = suppressWarnings(vcov(object))
        len = sapply(object$parameters, function(p) length(object[[paste0(p, ".coefficients")]]))
        before = which(object$parameters == what) - 1
        if (before > 0)  before = sum(len[seq_len(before)])
        idx = before + seq_along(object$coefficients)
        vcov. = V[idx, idx, drop = FALSE]
    }
    if (!is.null(link <- object[[paste0(what, ".link")]]))
        object$family = list(link = link)
    object$qr = object[[paste0(what, ".qr")]]
    NextMethod("emm_basis", vcov. = vcov., ...)
}
