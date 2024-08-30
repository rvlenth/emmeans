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

### Helper functions for emmeans
### Here we have 'recover_data' and 'emm_basis' methods
### For models that this package supports.

#--------------------------------------------------------------
### lm objects (and also aov, rlm, others that inherit) -- but NOT aovList
### Recent additional arhument 'frame' should point to where the model frame 
###   might be available, or NULL otherwise
#' @exportS3Method recover_data lm
#' @export
recover_data.lm = function(object, frame = object$model, ...) {
        fcall = object$call
    recover_data(fcall, delete.response(terms(object)), object$na.action, 
                 frame = frame, pwts = weights(object), ...)
}

#' @export
emm_basis.lm = function(object, trms, xlev, grid, ...) {
    # coef() works right for lm but coef.aov tosses out NAs
    bhat = object$coefficients
    nm = if(is.null(names(bhat))) row.names(bhat) else names(bhat)
    m = suppressWarnings(model.frame(trms, grid, na.action = na.pass, xlev = xlev))
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    assign = attr(X, "assign")
    X = X[, nm, drop = FALSE]
    bhat = as.numeric(bhat) 
    # stretches it out if multivariate - see mlm method
    V = .my.vcov(object, ...)
    
    if (sum(is.na(bhat)) > 0)
        nbasis = estimability::nonest.basis(object$qr)
    else
        nbasis = estimability::all.estble
    misc = list()
    if (inherits(object, "glm")) {
        misc = .std.link.labels(object$family, misc)
        dffun = function(k, dfargs) dfargs$df
        dfargs = list(df = ifelse(object$family$family %in% c("gaussian", "Gamma"), 
                                  object$df.residual, Inf))
    }
    else {
        dfargs = list(df = object$df.residual)
        dffun = function(k, dfargs) dfargs$df
    }
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc,
         model.matrix = .cmpMM(object$qr, assign = assign))
}



#--------------------------------------------------------------
### mlm objects
# (recover_data.lm works just fine)

#' @export
emm_basis.mlm = function(object, trms, xlev, grid, ...) {
    class(object) = c("mlm", "lm") # avoids error in vcov for "maov" objects
    bas = emm_basis.lm(object, trms, xlev, grid, ...)
    bhat = coef(object)
    k = ncol(bhat)
    bas$X = kronecker(diag(rep(1,k)), bas$X)
    bas$nbasis = kronecker(rep(1,k), bas$nbasis)
    ylevs = dimnames(bhat)[[2]]
    if (is.null(ylevs)) ylevs = seq_len(k)
    bas$misc$ylevs = list(rep.meas = ylevs)
    bas
}

#----------------------------------------------------------
# manova objects
#' @exportS3Method recover_data manova
recover_data.manova = function(object, ...) {
    fcall = match.call(aov, object$call)   # need to borrow arg matching from aov()
    recover_data(fcall, delete.response(terms(object)), object$na.action, 
                 frame = object$model, ...)
}




#--------------------------------------------------------------
### merMod objects (lme4 package)
#' @export
recover_data.merMod = function(object, ...) {
    if(!lme4::isLMM(object) && !lme4::isGLMM(object)) 
        return("Can't handle a nonlinear mixed model")
    fcall = object@call
    recover_data(fcall, delete.response(terms(object)), 
                 attr(object@frame, "na.action"), frame = object@frame, pwts = weights(object), ...)
}

#' @export
emm_basis.merMod = function(object, trms, xlev, grid, 
                            mode = get_emm_option("lmer.df"), lmer.df, 
                            disable.pbkrtest = get_emm_option("disable.pbkrtest"), 
                            pbkrtest.limit = get_emm_option("pbkrtest.limit"), 
                            disable.lmerTest = get_emm_option("disable.lmerTest"), 
                            lmerTest.limit = get_emm_option("lmerTest.limit"), 
                            options, ...) {
    V = .my.vcov(object, ...)
    dfargs = misc = list()
    
    if (lme4::isLMM(object)) {
        # Allow lmer.df in lieu of mode
        if (!missing(lmer.df))
            mode = lmer.df

        mode = match.arg(tolower(mode), c("satterthwaite", "kenward-roger", "asymptotic"))
        # if we're gonna override the df anyway, keep it simple 
        # OTOH, if K-R, documentation promises we'll adjust V
        if (!is.null(options$df) && (mode != "kenward-roger")) 
            mode = "asymptotic"
        
        
        # set flags
        objN = lme4::getME(object, "N")
        tooBig.k = (objN > pbkrtest.limit)
        tooBig.s = (objN > lmerTest.limit)
        
        tooBigMsg = function(pkg, limit) {  
            message("Note: D.f. calculations have been",
                    " disabled because the number of observations exceeds ", limit, ".\n",
                    "To enable adjustments, add the argument '", pkg, ".limit = ", objN, "' (or larger)\n",
                    "[or, globally, 'set emm_options(", pkg, ".limit = ", objN, ")' or larger];\n",
                    "but be warned that this may result in large computation time and memory use.")
        }

        # pick the lowest-hanging apples first
        if (mode == "kenward-roger") {
            if (disable.pbkrtest || tooBig.k || !.requireNS("pbkrtest", 
                    "Cannot use mode = \"kenward-roger\" because *pbkrtest* package is not installed", 
                    fail = message))
                mode = "satterthwaite"
            if (!disable.pbkrtest && tooBig.k)
                tooBigMsg("pbkrtest", pbkrtest.limit)
        }
        if (mode == "satterthwaite") {
            if (disable.lmerTest || tooBig.s || !.requireNS("lmerTest", 
                    "Cannot use mode = \"satterthwaite\" because *lmerTest* package is not installed", 
                    fail = message))
                mode = ifelse(!disable.pbkrtest && !tooBig.k && 
                                  .requireNS("pbkrtest", fail = .nothing), 
                              "kenward-roger", "asymptotic")
            if (!disable.lmerTest && tooBig.s)
                tooBigMsg("lmerTest", lmerTest.limit)
        }
        # if my logic isn't flawed, we are guaranteed that mode is both desired and possible
        
        if (mode == "kenward-roger") {
            if ((\(vcov., ...) missing(vcov.))(...)) {   # if (vcov. not in ...)
                dfargs = list(unadjV = V, 
                              adjV = pbkrtest::vcovAdj.lmerMod(object, 0))
                V = as.matrix(dfargs$adjV)
                tst = try(pbkrtest::Lb_ddf)
                if(!inherits(tst, "try-error"))
                    dffun = function(k, dfargs) pbkrtest::Lb_ddf (k, dfargs$unadjV, dfargs$adjV)
                else {
                    mode = "asymptotic"
                    warning("Failure in loading pbkrtest routines",
                            " - reverted to \"asymptotic\"")
                }
            }
            else {
                message("Kenward-Roger method can't be used with user-supplied covariances")
                mode = "satterthwaite"
            }
        }
        
        if (mode == "satterthwaite") {
            dfargs = list(object = object)
            dffun = function(k, dfargs) 
                suppressMessages(lmerTest::calcSatterth(dfargs$object, k)$denom)
        }
        
        if (mode == "asymptotic") {
            dffun = function(k, dfargs) Inf
        }
        
        attr(dffun, "mesg") = mode
    }
    else if (lme4::isGLMM(object)) {
        dffun = function(k, dfargs) Inf
        misc = .std.link.labels(family(object), misc)
    }
    else 
        stop("Can't handle a nonlinear mixed model")
    
    contrasts = attr(object@pp$X, "contrasts")
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    bhat = lme4::fixef(object)
    
    if (length(bhat) < ncol(X)) {
        # Newer versions of lmer can handle rank deficiency, but we need to do a couple of
        # backflips to put the pieces together right,
        # First, figure out which columns were retained
        kept = match(names(bhat), dimnames(X)[[2]])
        # Now re-do bhat with NAs in the right places
        bhat = NA * X[1, ]
        bhat[kept] = lme4::fixef(object)
        # we have to reconstruct the model matrix
        modmat = model.matrix(trms, object@frame, contrasts.arg=contrasts)
        nbasis = estimability::nonest.basis(modmat)
    }
    else
        nbasis=estimability::all.estble
    
    mm = .cmpMM(object@pp$X, object@pp$Xwts^2,
                attr(object@pp$X, "assign"))
    
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc,
         model.matrix = mm)
}




#--------------------------------------------------------------
### lme objects (nlme package)
#' @export
recover_data.lme = function(object, data, ...) {
    fcall = object$call
    if (!is.null(eval(fcall$weights, environment(terms(object)), parent.frame()))) {  # painful -- we only get weights for complete cases
        if (!is.null(object$na.action)) {
            w = nlme::varWeights(object$modelStruct)
            wts = rep(0, length(w) + length(object$na.action))
            wts[-object$na.action] = w
            fcall$weights = wts
        }
        else
            fcall$weights = nlme::varWeights(object$modelStruct)
    }
    dat = recover_data(fcall, delete.response(object$terms), object$na.action, data = data, ...)
    attr(dat, "pass.it.on") = TRUE
    dat
}

#' @export
emm_basis.lme = function(object, trms, xlev, grid, 
        mode = c("containment", "satterthwaite", "appx-satterthwaite", "auto", "boot-satterthwaite", "asymptotic"), 
        sigmaAdjust = TRUE, options, extra.iter = 0, ...) {
    mode = match.arg(mode)
    if (mode == "boot-satterthwaite") mode = "appx-satterthwaite"  # backward compatibility
    if (mode == "asymptotic")
        options$df = Inf
    if (!is.null(options$df)) # if we're gonna override the df anyway, keep it simple!
        mode = "fixed"
    if (mode == "auto")
        mode = ifelse(is.null(object$apVar), "containment", "appx-satterthwaite")
    if (is.null(object$apVar))
        mode = "containment"
    contrasts = object$contrasts
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    bhat = nlme::fixef(object)
    V = .my.vcov(object, ...)
    if (sigmaAdjust && object$method == "ML") 
        V = V * object$dims$N / (object$dims$N - nrow(V))
    misc = list()
    if (!is.null(object$family)) {
        misc = .std.link.labels(object$family, misc)
    }
    nbasis = estimability::all.estble
    
    if (mode == "fixed") { # hack to just put in df from options
        dfargs = list(df = options$df)
        dffun = function(k, dfargs) dfargs$df
    }
    else if (mode %in% c("satterthwaite", "appx-satterthwaite")) {
        mode = "appx-satterthwaite"
        G = try(gradV.kludge(object, extra.iter = extra.iter), silent = TRUE)
        ###! not yet, doesn't work G = try(lme_grad(object, object$call, object$data, V))
        if (inherits(G, "try-error"))
            stop("Unable to estimate Satterthwaite parameters")
        dfargs = list(V = V, A = object$apVar, G = G)
        dffun = function(k, dfargs) {
            est = tcrossprod(crossprod(k, dfargs$V), k)
            g = sapply(dfargs$G, function(M) tcrossprod(crossprod(k, M), k))
            varest = tcrossprod(crossprod(g, dfargs$A), g)
            2 * est^2 / varest
        }
    }
    else { # containment df
        dfx = object$fixDF$X
        if (names(bhat[1]) == "(Intercept)")
            dfx[1] = length(levels(object$groups[[1]])) - 1
        ### Correct apparent error in lme containment algorithm
        dffun = function(k, dfargs) {
            idx = which(abs(k) > 1e-4)
            ifelse(length(idx) > 0, min(dfargs$dfx[idx]), NA)
        }
        dfargs = list(dfx = dfx)
    }
    attr(dffun, "mesg") = mode
    
    # submodel support (not great -- omits any weights)
    m = model.frame(trms, attr(object, "data"), na.action = na.pass, xlev = xlev)
    mm = model.matrix(trms, m, contrasts.arg = contrasts)
    mm = .cmpMM(mm, assign = attr(mm, "assign"))
    
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, 
         dffun = dffun, dfargs = dfargs, misc = misc,
         model.matrix = mm)
}

# Here is a total hack, but it works pretty well
# We estimate the gradient of the V matrix by fitting the
# model with a few random perturbations of y, then
# regressing the changes in V against the changes in the 
# covariance parameters
gradV.kludge = function(object, Vname = "varFix", call = formula(object$terms), 
                        data = object$data, extra.iter = 0) {
    # check consistency of contrasts
    #### This code doesn't work with coerced factors. Hardly seems messing with, so I commented it out
    # cnm = names(object$contrasts)
    # cdiff = sapply(cnm, function(.) max(abs(contrasts(data[[.]]) - object$contrasts[[.]])))
    # if (max(cdiff) > 1e-6) {
    #     message("Contrasts don't match those used when the model was fitted. Fix this and re-run")
    #     stop()
    # }
    
    if(is.null(data)) {
        vars = all.vars(eval(object$call[[2]]))
        lst = lapply(vars, get)
        names(lst) = vars
        data = data.frame(lst)
    }
    A = object$apVar
    theta = attr(A, "Pars")
    V = object[[Vname]]
    sig = .01 * object$sigma
    #data = object$data
    yname = all.vars(eval(call))[1]
    y = data[[yname]]
    n = length(y)
    object$call[[2]] = call
    # we're gonna carefully check simulation results to make sure we didn't
    # lose (or gain?) any params so that everything conforms
    niter = 0
    nsim = 2 + extra.iter + length(theta)
    dat = t(replicate(nsim, {
        simt = simv = numeric(0)
        niter = niter + 1
        while((niter < 3 * nsim) && 
              ((length(simv) != length(V)) || (length(simt) != length(theta)))) {
            data[[yname]] = y + sig * rnorm(n)
            mod = update(object, data = data)
            simt = attr(mod$apVar, "Pars")
            simv = mod[[Vname]]
        }
        niter = niter + 1
        if (niter > 3 * nsim) 
            stop("Too many simulations due to inconsistency")
        c(simt - theta, as.numeric(simv - V))
    }))
    dimnames(dat) = c(NULL, NULL)
    xcols = seq_along(theta)
    B = lm.fit(dat[, xcols], dat[,-xcols])$coefficients
    grad = lapply(seq_len(nrow(B)), function(i) matrix(B[i, ], nrow=nrow(V)))
    grad
}

# ### new way to get gradients for lme models
# # (not ready for primetime...)
# lme_grad = function(object, call, data, V) {
#     obj = object$modelStruct
#     conLin = object
#     class(conLin) = class(obj)
#     X = model.matrix(eval(call$fixed), data = data)
#     y = data[[all.vars(call)[1]]]
#     conLin$Xy = cbind(X, y)
#     conLin$fixedSigma = FALSE
#     grps = object$groups # May have to re-order these fancily?
#     MEest = get("MEestimate", getNamespace("nlme")) ## workaround its not being exported
#     func = function(x) {
#         coef(obj) = x
#         tmp = MEest(obj, grps, conLin)
#         crossprod(tmp$sigma * tmp$varFix)
#     }
#     res = numDeriv::jacobian(func, coef(obj))
#     G = lapply(seq_len(ncol(res)), function(j) matrix(res[, j], ncol = ncol(V)))
#     G[[1 + length(G)]] = 2 * V  # gradient wrt log sigma
#     G
# }



#--------------------------------------------------------------

### new way to get jacobians for gls models
gls_grad = function(object, call, data, V) {
    if (is.null(data)) 
        data = environment(eval(call$model))
    obj = object$modelStruct
    conLin = object
    class(conLin) = class(obj)
    X = model.matrix(eval(call$model), data = data)
    y = eval(eval(call$model)[[2]], envir = data)
    conLin$Xy = cbind(X, y)
    conLin$fixedSigma = FALSE
    func = function(x) {
        obj = nlme::`coef<-`(obj, value = x)
        tmp = nlme::glsEstimate(obj, conLin)
        .get.lt(crossprod(tmp$sigma * tmp$varBeta))  # lower triangular form
    }
    res = numDeriv::jacobian(func, coef(obj))
    G = lapply(seq_len(ncol(res)), function(j) .lt2mat(res[, j]))
    G[[1 + length(G)]] = 2 * V  # gradient wrt log sigma
    G
}

### gls objects (nlme package)
#' @exportS3Method recover_data gls
recover_data.gls = function(object, data, ...) {
    fcall = object$call
    if (!is.null(wts <- fcall$weights)) {
        wts = nlme::varWeights(object$modelStruct)
        fcall$weights = NULL
    }
    trms = delete.response(terms(nlme::getCovariateFormula(object)))
    attr(trms, "predvars") = attr(delete.response(terms(object)), "predvars")
    # above copies scaling info for scale() and poly() terms
    
    result = recover_data.call(fcall, trms, object$na.action, data = data, ...)
    if (!is.null(wts))
        result[["(weights)"]] = wts
    if (!missing(data))
        attr(result, "misc") = list(data = data)
    attr(result, "pass.it.on") = TRUE
    result
}

#' @exportS3Method emm_basis gls          
emm_basis.gls = function(object, trms, xlev, grid, 
                         mode = c("auto", "df.error", "satterthwaite", "appx-satterthwaite", "boot-satterthwaite", "asymptotic"), 
                         extra.iter = 0, options, misc, ...) {
    contrasts = object$contrasts
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)

    # submodel support (not great because I don't know how to retrieve the weights)
    mf = model.frame(trms, attr(object, "data"), na.action = na.pass, xlev = xlev)
    mm = model.matrix(trms, mf, contrasts.arg = contrasts)
    mm = .cmpMM(mm, assign = attr(mm, "assign"))
    
    tmp = coef(object)
    bhat = rep(NA, ncol(X))
    bhat[match(names(tmp), colnames(X), nomatch = 0)] = tmp
    V = .my.vcov(object, ...)
    if (any(is.na(bhat)))
        nbasis = estimability::nonest.basis(mm)
    else
        nbasis = estimability::all.estble
    mode = match.arg(mode)
    if (mode == "boot-satterthwaite") mode = "appx-satterthwaite"  # backward compatibility
    if (!is.null(options$df)) # if we're gonna override the df anyway, keep it simple!
        mode = "df.error"
    if (mode == "auto")
        mode = ifelse(is.null(object$apVar), "df.error", "satterthwaite")
    if (!is.matrix(object$apVar))
        mode = "df.error"
    if (mode %in% c("satterthwaite", "appx-satterthwaite")) {
        data = if(is.null(misc$data))
            eval(object$call$data, parent.frame(2))
        else
            misc$data
        misc = list()
        chk = attr(object$apVar, "Pars")
        if(max(abs(coef(object$modelStruct) - chk[-length(chk)])) > .001) {
            message("Analytical Satterthwaite method not available; using appx-satterthwaite")
            mode = "appx-satterthwaite"
        }
        if (mode == "appx-satterthwaite") {
            G = try(gradV.kludge(object, "varBeta", call = object$call$model,
                                 data = data, extra.iter = extra.iter),
                    silent = TRUE)
        }
        else
            G = try(gls_grad(object, object$call, data, V))
        if (inherits(G, "try-error")) {
            sugg = ifelse(mode == "satterthwaite", "appx-satterthwaite", "df.error")
            stop("Can't estimate Satterthwaite parameters.\n",
                 "  Try adding the argument 'mode = \"", sugg, "\"'", call. = FALSE)
        }
        dfargs = list(V = V, A = object$apVar, G = G)
        dffun = function(k, dfargs) {
            est = tcrossprod(crossprod(k, dfargs$V), k)
            g = sapply(dfargs$G, function(M) tcrossprod(crossprod(k, M), k))
            varest = tcrossprod(crossprod(g, dfargs$A), g)
            2 * est^2 / varest
        }
    }
    else if (mode %in%  c("df.error", "asymptotic")) {
        df = ifelse(mode == "asymptotic",
                    Inf, 
                    object$dims$N - object$dims$p - length(unlist(object$modelStruct)))
        dfargs = list(df = df)
        dffun = function(k, dfargs) dfargs$df
    }
    attr(dffun, "mesg") = mode
    
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc,
         model.matrix = mm)
}



#--------------------------------------------------------------
### polr objects (MASS package)
#' @exportS3Method recover_data polr
recover_data.polr = function(object, ...)
    recover_data.clm(object, ...)

#' @exportS3Method emm_basis polr         
emm_basis.polr = function(object, trms, xlev, grid, 
                          mode = c("latent", "linear.predictor", "cum.prob", "exc.prob", "prob", "mean.class"), 
                          rescale = c(0,1), ...) {
    mode = match.arg(mode)
    contrasts = object$contrasts
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    # Strip out the intercept (borrowed code from predict.polr)
    xint = match("(Intercept)", colnames(X), nomatch = 0L)
    if (xint > 0L) 
        X = X[, -xint, drop = FALSE]
    bhat = c(coef(object), object$zeta)
    V = .my.vcov(object, ...)
    k = length(object$zeta)
    if (mode == "latent") {
        X = rescale[2] * cbind(X, matrix(- 1/k, nrow = nrow(X), ncol = k))
        bhat = c(coef(object), object$zeta - rescale[1] / rescale[2])
        misc = list(offset.mult = rescale[2])
    }
    else {
        j = matrix(1, nrow=k, ncol=1)
        J = matrix(1, nrow=nrow(X), ncol=1)
        X = cbind(kronecker(-j, X), kronecker(diag(1,k), J))
        link = object$method
        if (link == "logistic") link = "logit"
        misc = list(ylevs = list(cut = names(object$zeta)), 
                    tran = link, inv.lbl = "cumprob", offset.mult = -1)
        if (mode != "linear.predictor") {
            # just use the machinery we already have for the 'ordinal' package
            misc$mode = mode
            misc$postGridHook = ".clm.postGrid"
        }
    }
    misc$respName = as.character.default(terms(object))[2]
    nbasis = estimability::all.estble
    dffun = function(...) Inf
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=list(), misc=misc)
}




#--------------------------------------------------------------
### survreg objects (survival package)
#' @exportS3Method recover_data survreg
recover_data.survreg = function(object, ...) {
    fcall = object$call
    trms = delete.response(terms(object))
    # I'm gonna delete any terms involving cluster(), frailty(); keep strata()
    mod.elts = dimnames(attr(trms, "factor"))[[2]]
    tmp = grep("cluster\\(|frailty", mod.elts)
    if (length(tmp))
        trms = trms[-tmp]
    recover_data(fcall, trms, object$na.action, pwts = weights(object), ...)
}

# Seems to work right in a little testing.
# However, it fails sometimes if I update the model 
# with a subset argument. Workaround: just fitting a new model
#' @exportS3Method emm_basis survreg      
emm_basis.survreg = function(object, trms, xlev, grid, ...) {
    # Much of this code is adapted from predict.survreg
    bhat = object$coefficients
    k = length(bhat) - sum(is.na(bhat))
    V = .my.vcov(object, ...)[seq_len(k), seq_len(k), drop=FALSE]
    # ??? not used... is.fixeds = (k == ncol(object$var))
    ### zap-out factors in xlev not needed by model.frame
    xlev[setdiff(names(xlev), rownames(attr(trms, "factors")))] = NULL
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)    
    # X = model.matrix(object, m) # This is what predict.survreg does
    # But I have manipulated trms, so need to make sure things are consistent
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    usecols = intersect(colnames(X), names(bhat))
    bhat = bhat[usecols]  # in case ref_grid code excluded some levels...
    X = X[, usecols, drop=FALSE]
    nbasis = estimability::nonest.basis(model.matrix(object)[, usecols, drop = FALSE])
    dfargs = list(df = object$df.residual)
    dffun = function(k, dfargs) dfargs$df
    if (object$dist %in% c("exponential","weibull","loglogistic","loggaussian","lognormal")) 
        misc = list(tran = "log", inv.lbl = "response")
    else 
        misc = list()
    misc$postGridHook = .notran2   # removes "Surv()" as response transformation
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}



#--------------------------------------------------------------
###  coxph objects (survival package)
#' @exportS3Method recover_data coxph
recover_data.coxph = function(object, ...) 
    recover_data.survreg(object, ...)

#' @exportS3Method emm_basis coxph        
emm_basis.coxph = function (object, trms, xlev, grid, ...) 
{
    object$dist = "doesn't matter"
    result = emm_basis.survreg(object, trms, xlev, grid, ...)
    result$dfargs$df = Inf
    nms = colnames(result$X)
    # delete columns for intercept and main effects of strata
    zaps = which(nms %in% setdiff(nms, names(result$bhat)))
    if(length(zaps) > 0)
        result$X = result$X[, -zaps, drop = FALSE]
    ### result$X = result$X - rep(object$means, each = nrow(result$X))
    result$misc$tran = "log"
    result$misc$inv.lbl = "hazard"
    result
}

.notran2 = function(object, ...) {
    for (nm in c("tran", "tran2"))
        if(!is.null(object@misc[[nm]]) && object@misc[[nm]] == "Surv") object@misc[[nm]] = NULL
    object
}

# Note: Very brief experimentation suggests coxph.penal also works.
# This is an extension of coxph


#--------------------------------------------------------------
###  coxme objects ####
### Greatly revised 6-15-15 (after version 2.18)
#' @exportS3Method recover_data coxme
recover_data.coxme = function(object, ...) 
    recover_data.survreg(object, ...)

#' @exportS3Method emm_basis coxme        
emm_basis.coxme = function(object, trms, xlev, grid, ...) {
    bhat = coxme::fixef(object)
    k = length(bhat)
    V = .my.vcov(object, ...)[seq_len(k), seq_len(k), drop = FALSE]
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m)
    X = X[, -1, drop = FALSE] # remove the intercept
    # scale the linear predictor
    for (j in seq_along(X[1, ]))
        X[, j] = (X[, j] - object$means[j]) ### / object$scale[j]
    nbasis = estimability::all.estble
    dffun = function(k, dfargs) Inf
    misc = list(tran = "log", inv.lbl = "hazard")
    misc$postGridHook = .notran2   # removes "Surv()" as response transformation
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
         dfargs = list(), misc = misc)
}


###  special vcov prototype for cases where there are several vcov options
###  e.g., gee, geeglm, geese
.named.vcov = function(object, method, ...)
    UseMethod(".named.vcov")

# default has optional idx of same length as valid and if so, idx indicating 
#   which elt of valid to use if matched
# Ex: valid = c("mammal", "fish", "rat", "dog", "trout", "perch")
#     idx   = c(   1,        2,     1,     1,       2,       2)
#     -- so ultimately results can only be "mammal" or "fish"
# nonmatches revert to 1st elt.
.named.vcov.default = function(object, method, valid, idx = seq_along(valid), ...) {
    if (!is.character(method)) { # in case vcov. arg was matched by vcov.method {
        V = .my.vcov(object, method)
        method = "user-supplied"
    }
    else {
        i = pmatch(method, valid, 1)
        method = valid[idx[i]]
        V = object[[method]]
    }
    attr(V, "methMesg") = paste("Covariance estimate used:", method)
    V
}

# general-purpose emm_basis function for GEEs
.emmb.geeGP = function(object, trms, xlev, grid, vcov.method, valid, idx = seq_along(valid), ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    bhat = coef(object)
    V = .named.vcov(object, vcov.method, valid, idx, ...)
    
    if (sum(is.na(bhat)) > 0)
        nbasis = estimability::nonest.basis(object$qr)
    else
        nbasis = estimability::all.estble
    
    misc = .std.link.labels(object$family, list())
    misc$initMesg = attr(V, "methMesg")
    dffun = function(k, dfargs) Inf
    dfargs = list()
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}

#---------------------------------------------------------------
###  gee objects  ####


#' @exportS3Method recover_data gee
recover_data.gee = function(object, ...)
    recover_data.lm(object, frame = NULL, ...)

#' @exportS3Method emm_basis gee          
emm_basis.gee = function(object, trms, xlev, grid, vcov.method = "robust.variance", ...)
    .emmb.geeGP(object, trms, xlev, grid, vcov.method, 
                valid = c("robust.variance", "naive.variance"))

###  geepack objects  ####
#' @exportS3Method recover_data geeglm
recover_data.geeglm = function(object, ...)
    recover_data.lm(object, ...)

#' @exportS3Method emm_basis geeglm       
emm_basis.geeglm = function(object, trms, xlev, grid, vcov.method = "vbeta", ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    bhat = coef(object)
    V = .named.vcov(object$geese, vcov.method, 
                    valid = c("vbeta", "vbeta.naiv","vbeta.j1s","vbeta.fij","robust","naive"), 
                    idx = c(1,2,3,4,1,2))
    
    if (sum(is.na(bhat)) > 0)
        nbasis = estimability::nonest.basis(object$qr)
    else
        nbasis = estimability::all.estble
    
    misc = .std.link.labels(object$family, list())
    misc$initMesg = attr(V, "methMesg")
    dffun = function(k, dfargs) dfargs$df
    dfargs = list(df = df.residual(object))
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}


#' @exportS3Method recover_data geese
recover_data.geese = function(object, ...) {
    fcall = object$call
    # what a pain - we need to reconstruct the terms component
    args = as.list(fcall[-1])
    na.action = object$na.action
    #trms = terms.formula(fcall$formula)
    if (!is.null(args$data)) {
        data = eval(args$data, parent.frame())
        trms = terms(model.frame(fcall$formula, data = data))
    } else {
        trms = terms(model.frame(fcall$formula))
    }
    recover_data(fcall, delete.response(trms), na.action, ...)
}

#' @exportS3Method emm_basis geese        
emm_basis.geese = function(object, trms, xlev, grid, vcov.method = "vbeta", ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    bhat = object$beta
    V = .named.vcov(object, vcov.method, 
                    valid = c("vbeta", "vbeta.naiv","vbeta.j1s","vbeta.fij","robust","naive"), 
                    idx = c(1,2,3,4,1,2))

    # We don't have the qr component - I'm gonna punt for now
     if (sum(is.na(bhat)) > 0)
         warning("There are non-estimable functions, but estimability is NOT being checked")
#         nbasis = estimability::nonest.basis(object$qr)
#     else
        nbasis = estimability::all.estble
    
    misc = list()
    if (!is.null(fam <- object$call$family))
        misc = .std.link.labels(eval(fam)(), misc)
    misc$initMesg = attr(V, "methMesg")
    dffun = function(k, dfargs) Inf
    dfargs = list()
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, misc=misc)
}


### glmtoolbox package -- we just pretend it's lm and it seems to work

#' @exportS3Method recover_data glmgee
recover_data.glmgee = function(object, ...) {
    recover_data.lm(object, ...)
}

#' @exportS3Method emm_basis glmgee
emm_basis.glmgee = function(object, trms, xlev, grid, vcov.method = "robust", ...) {
    vcov. = if(is.character(vcov.method))
        vcov(object, type = vcov.method, ...)
    else
        vcov.method
    class(object) = c("glm", "lm") # pretend it's glm
    rtn = emm_basis(object, trms = trms, xlev = xlev, grid = grid, vcov. = vcov., ...)
    rtn$misc$initMesg = paste("Covariance estimate used:",
        ifelse(is.character(vcov.method), vcov.method, "user-supplied"))
    rtn
}



### survey package
# svyglm class
#' @exportS3Method recover_data svyglm
recover_data.svyglm = function(object, data = NULL, ...) {
    if (is.null(data)) {
        env = environment(terms(object))
        des = eval(object$call$design, envir = env)
        data = eval(des$call$data, envir = env)
    }
    recover_data.lm(object, data = data, frame = object$model, ...)
}

# inherited emm_basis.lm method works fine




### ----- Auxiliary routines -------------------------
# Provide for vcov. argument in ref_grid call, which could be a function or a matrix

.statsvcov = function(object, ...)
    stats::vcov(object, complete = FALSE, ...)

#' @rdname extending-emmeans
#' @order 36
#' @param vcov. Function or matrix that returns a suitable covariance matrix.
#' The default is \code{.statsvcov} which is \code{stats::vcov}. The \code{.my.vcov} 
#' function should be called in place of \code{\link{vcov}}, and it supports the user 
#' being able to specify a different matrix or function via the
#' optional \code{vcov.} argument.
#' @export
.my.vcov = function(object, vcov. = .statsvcov, ...) {
    if (is.function(vcov.))
        vcov. = vcov.(object, ...)
    vcov. = try(as.matrix(vcov.), silent = TRUE)
    if (inherits(vcov., "try-error") || (nrow(vcov.) != ncol(vcov.)))
        stop("vcov. must be a function or a square matrix")
    vcov.
}

#' @rdname extending-emmeans
#' @order 30
#' @param fam Result of call to \code{family(object)}
#' @param misc A \code{list} intended for the \code{@misc} slot of an \code{emmGrid} object
#' @return \code{.std.link.llabels} returns a modified version of \code{misc}
#'   with the appropriate information included corresponding to the information in \code{fam}
#' @export
.std.link.labels = function(fam, misc) {
    if (is.null(fam) || !is.list(fam))
        return(misc)
    if (fam$link == "identity")
        return(misc)
    misc$tran = fam$link
    misc$inv.lbl = "response"
    if (length(grep("binomial", fam$family)) == 1)
        misc$inv.lbl = "prob"
    else if (length(grep("poisson", fam$family)) == 1)
        misc$inv.lbl = "rate"
    if(length(grep("gaussian", fam$family)) == 0)
        misc$sigma = NA
    misc
}


