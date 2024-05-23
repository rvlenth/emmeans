##############################################################################
#    Copyright (c) 2012-2022 Russell V. Lenth                                #
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

# We have one args:
#   nboot   # of bootstrap reps to get variances of smooths
#' @exportS3Method emm_basis Gam
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


### This addition contributed by Hannes Riebl (#303)
.emm_basis.gam_multinom = function(object, trms, xlev, grid,
                                   freq = FALSE, unconditional = FALSE,
                                   mode = c("prob", "latent"), ...) {
    mode = match.arg(mode)

    X = mgcv::predict.gam(object, newdata = grid, type = "lpmatrix",
                          newdata.guaranteed = TRUE)

    k = length(attr(X, "lpi"))
    nbhat = vapply(attr(X, "lpi"), length, FUN.VALUE = integer(1))
    pat = (rbind(0, diag(k + 1, k)) - 1) / (k + 1)

    X = apply(pat, 1, function(row) {
        y = rep.int(row, times = nbhat)
        out = apply(X, 1, "*", y = y, simplify = FALSE)
        do.call(rbind, out)
    }, simplify = FALSE)

    X = do.call(rbind, X)

    bhat = as.numeric(coef(object))
    V = .my.vcov(object, freq = freq, unconditional = unconditional, ...)
    nbasis = kronecker(rep.int(1, times = k), estimability::all.estble)

    dfargs = list(df = sum(object$edf))
    dffun = function(k, dfargs) dfargs$df

    misc = list(tran = "log", inv.lbl = "e^y")

    ylevs = list(class = seq.int(0, k))
    names(ylevs) = as.character(object$formula[[1]][[2]])
    misc$ylevs = ylevs

    if (mode == "prob") {
        misc$postGridHook = .multinom.postGrid
    }

    list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun,
         dfargs = dfargs, misc = misc)
}

# Local utility for identifying random smooths
.smooth.is.random = function(s) {
    rcls = c("random.effect", "fs.interaction") ### what to look for
    cls = c(class(s), unname(unlist(lapply(s$margin, class))))
    any(cls %in% rcls)
}

### for mgcv::gam objects
### Many thanks to Maarten Jung for help on sorting-out fixed and random effects
#' @exportS3Method recover_data gam
recover_data.gam = function(object, ...) {
    if (length(object$smooth) > 0) { # get rid of random terms
        fixnm = unlist(lapply(object$smooth, function(s) {
                              if(.smooth.is.random(s)) ""
                              else c(s$term, s$by)
            }))
        fixnm = union(.all.vars(delete.response(object$pterms)), fixnm)
        fixnm = setdiff(fixnm, c("1", "", "NA"))
        object$terms = terms(.reformulate(fixnm, env = environment(terms(object))))
    }
    recover_data.lm(object, ...)
}


### emm_basis method for mgcv::gam objects
### extra arg `unconditional` and `freq` as in `vcov.gam`
#' @exportS3Method emm_basis gam          
emm_basis.gam = function(object, trms, xlev, grid,
                         freq = FALSE, unconditional = FALSE,
                         what = c("location", "scale", "shape", "rate", "prob.gt.0"), 
                         ...) {
    if (length(object$smooth) > 0) { # get rid of random terms 
        rand = sapply(object$smooth, function(s) {ifelse(.smooth.is.random(s), s$label, NA)})
        rand = if (all(is.na(rand))) NULL else rand[!is.na(rand)]
    }
    else
        rand = NULL
    X = mgcv::predict.gam(object, newdata = grid, type = "lpmatrix",
                          exclude = rand, newdata.guaranteed = TRUE)
    keep = if (is.null(rand)) rep(TRUE, ncol(X)) else apply(X, 2, function(x) !all(x == 0))
    X = X[, keep, drop = FALSE]
    bhat = as.numeric(object$coefficients[keep])
    V = .my.vcov(object, freq = freq, unconditional = unconditional, ...)[keep, keep]

    fam_name = object$family$family
    what_num = what
    
    if (fam_name == "multinom") {
        return(.emm_basis.gam_multinom(object, trms, xlev, grid, freq,
                                       unconditional, ...))
    }
    else if (fam_name == "mvn") {
        if (!is.numeric(what)) {
            stop("Family 'mvn' requires a numeric argument 'what'")
        }
    } 
    else if (is.character(what)) {
        what = match.arg(what)
        if (fam_name == "ziplss") {
            what_num = switch(what, location = 1, rate = 1, prob.gt.0 = 2)
        } 
        else {
            what_num = switch(what, location = 1, scale = 2, shape = 3)
        }
    }
    
    select = attr(X, "lpi")
    if (is.null(select)) select = list(seq_along(bhat))
    select = try(select[[what_num]], silent = TRUE)
    
    if (inherits(select, "try-error")) {
        stop("Model does not have a linear predictor 'what = ", what, "'")
    }
    
    bhat = bhat[select]
    X = X[, select, drop = FALSE]
    V = V[select, select, drop = FALSE]
    
    nbasis = estimability::all.estble
    link = object$family$link[what_num]
    if(link == "identity") # they may be lying
        link = switch(fam_name,
                      ocat = "logit",
                      ziP = "log",
                      cox.ph = "log", ## ???
                      ziplss = c("log", "cloglog")[what_num],
                      gevlss = c("identity", "log", "logit")[what_num],
                      "identity")
    misc = .std.link.labels(list(link = link, family = fam_name), list())
    
    if (!is.null(misc$tran) && misc$tran == "logb")  # the way this is documented is truly bizarre but I think this is right
        misc$tran = make.tran("genlog", - environment(object$family$linfo[[what_num]]$linkfun)$b)

    dfargs = list(df = object$df.residual)
    dffun = function(k, dfargs) dfargs$df
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}


### mgcv::gamm objects...
#' @exportS3Method recover_data gamm
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

#' @exportS3Method emm_basis gamm         
emm_basis.gamm = function(object, ...)
    emm_basis(object$gam, ...)


###===================================================================
# Support for gamlss objects
# 'what' parameter mimics predict.gamlss

#' @exportS3Method recover_data gamlss
recover_data.gamlss = function(object, what = c("mu", "sigma", "nu", "tau"), ...) {
    fcall = object$call
    what = match.arg(what)
    trms = terms(formula(object, what = what))
    recover_data(fcall, delete.response(trms), object$na.action, ...)
}

#' @exportS3Method emm_basis gamlss       
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
    if (!is.null(link <- object[[paste0(what, ".link")]])) {
        # Decide whether to use d.f. or not
        use.df = c("BCCG", "BCPE", "BCT", "GA", "GT", "NO", "NOF", "TF")
        fam = ifelse((what == "mu") && (object$family[1] %in% use.df), 
                     "gaussian", "other")
        object$family = list(family = fam, link = link)
    }
    ###object$qr = object[[paste0(what, ".qr")]]
    ###NextMethod("emm_basis", vcov. = vcov., ...)
    emm_basis.lm(object, trms, xlev, grid, vcov. = vcov., ...)
}
