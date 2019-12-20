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

### support for the ordinal package

recover_data.clm = function(object, mode = "latent", ...) {
    if (!is.na(pmatch(mode, "scale"))) {
        if (is.null(trms <- object$S.terms))
            return("Specified mode=\"scale\", but no scale model is present") # ref_grid's error handler takes it from here
        recover_data(object$call, trms, object$na.action, ...)
    }
    else if (is.null(object$S.terms) && is.null(object$nom.terms))
        recover_data.lm(object, ...)
    else { # bring-in predictors from loc, scale, and nom models
        trms = delete.response(object$terms)
        x.preds = union(.all.vars(object$S.terms), .all.vars(object$nom.terms))
        x.trms = terms(update(trms, .reformulate(c(".", x.preds))))
        recover_data(object$call, x.trms, object$na.action, ...)
    }
}

# For now at least, clmm doesn't cover scale, nominal options
recover_data.clmm = recover_data.lm

# Note: For ALL thresholds, object$Theta has all the threshold values
# for the different cuts (same as object$alpha when threshold=="flexible")
# and object$tJac is s.t. tJac %*% alpha = Theta
# Note also that some functions of cut are constrained to be zero when
# threshold != "flexible". Can get basis using nonest.basis(t(tJac))
#
# opt arg 'mode' - determines what goes into ref_grid
#         'rescale' - (loc, scale) for linear transformation of latent result

emm_basis.clm = function (object, trms, xlev, grid, 
                          mode = c("latent", "linear.predictor", "cum.prob", "exc.prob", "prob", "mean.class", "scale"), 
                          rescale = c(0,1), ...) {
    # general stuff
    mode = match.arg(mode)
    if (mode == "scale")
        return (.emm_basis.clm.scale(object, trms, xlev, grid, ...))
    
    # if (is.null(object$contrasts))
    #     warning("Contrasts used to fit the model are unknown.\n",
    #             "Defaulting to system option, but results may be wrong.")
    
    bhat = coef(object)
    V = .my.vcov(object, ...)
    tJac = object$tJac
    dffun = function(...) Inf
    link = as.character(object$info$link)
    cnm = dimnames(object$tJac)[[1]]
    if (is.null(cnm))
        cnm = paste(seq_len(nrow(tJac)), "|", 1 + seq_len(nrow(tJac)), sep = "")
    misc = list()
    
    # My strategy is to piece together the needed matrices for each threshold parameter
    # Then assemble the results
    
    ### ----- Location part ----- ###
    contrasts = object$contrasts
    # Remember trms was trumped-up to include scale and nominal predictors.
    # Recover the actual terms for the principal model
    trms = delete.response(object$terms)
    m = model.frame(trms, grid, na.action = na.pass, xlev = object$xlevels)
    X = model.matrix(trms, m, contrasts.arg = contrasts)
    xint = match("(Intercept)", colnames(X), nomatch = 0L)
    if (xint > 0L) {
        X = X[, -xint, drop = FALSE]
    }
    
    ### ----- Nominal part ----- ###
    if (is.null(object$nom.terms))
        NOM = matrix(1, nrow = nrow(X))
    else {
        mn = model.frame(object$nom.terms, grid, na.action = na.pass, xlev = object$nom.xlevels)
        NOM = model.matrix(object$nom.terms, mn, contrasts.arg = object$nom.contrasts)
    }
    bigNom = kronecker(tJac, NOM)
    # cols are in wrong order... I'll get the indexes by transposing a matrix of subscripts
    if (ncol(NOM) > 1)
        bigNom = bigNom[, as.numeric(t(matrix(seq_len(ncol(bigNom)), nrow=ncol(NOM))))]
    
    ### ----- Scale part ----- ###
    if (!is.null(object$S.terms)) {
        ms = model.frame(object$S.terms, grid, na.action = na.pass, xlev = object$S.xlevels)
        S = model.matrix(object$S.terms, ms, contrasts.arg = object$S.contrasts)
        S = S[, names(object$zeta), drop = FALSE]
        if (!is.null(attr(object$S.terms, "offset"))) {
            soff = .get.offset(object$S.terms, grid)
            # we'll add a column to S and adjust bhat and V accordingly
            S = cbind(S, offset = soff)
            bhat = c(bhat, offset = 1)
            V = rbind(cbind(V, offset = 0), offset = 0)
        }
        si = misc$scale.idx = length(object$alpha) + length(object$beta) + seq_len(ncol(S))
        # Make sure there are no name clashes
        names(bhat)[si] = paste(".S", names(object$zeta), sep=".")
        misc$estHook = ".clm.estHook"
        misc$vcovHook = ".clm.vcovHook"
    }
    else
        S = NULL
    
    ### ----- Get non-estimability basis ----- ###
    nbasis = snbasis = estimability::all.estble
    if (any(is.na(bhat))) {
        mm = model.matrix(object)
        # note: mm has components X, NOM, and S
        if (any(is.na(c(object$alpha, object$beta)))) {
            NOMX = if (is.null(mm$NOM)) mm$X
                   else                 cbind(mm$NOM, mm$X[, -1])
            nbasis = estimability::nonest.basis(NOMX)
            # replicate and reverse the sign of the NOM parts
            nomcols = seq_len(ncol(NOM))
            nbasis = apply(nbasis, 2, function(x)
                c(rep(-x[nomcols], each = nrow(NOM)), x[-nomcols]))
        }
        if (!is.null(mm$S)) {
            if (any(is.na(object$zeta))) {
                snbasis = estimability::nonest.basis(mm$S)
                # put intercept part at end
                snbasis = rbind(snbasis[-1, , drop=FALSE], snbasis[1, ])
                if (!is.null(attr(object$S.terms, "offset")))
                    snbasis = rbind(snbasis, 0)
                snbasis = rbind(matrix(0, ncol=ncol(snbasis), nrow=min(si)-1), snbasis)
                 # Note scale intercept is included, so tack it on to the end of everything
                S = cbind(S, .S.intcpt = 1)
                bhat = c(bhat, .S.intcpt = 0)
                V = rbind(cbind(V, .S.intcpt = 0), .S.intcpt = 0)
                si = misc$scale.idx = c(si, 1 + max(si))
            }
        }
        if (is.na(nbasis[1])) # then only nonest part is scale
            nbasis = snbasis
        else { 
            if (!is.null(S)) # pad nbasis with zeros when there's a scale model
                nbasis = rbind(nbasis, matrix(0, nrow=length(si), ncol=ncol(nbasis)))
            if (!is.na(snbasis[1]))
                nbasis = cbind(nbasis, snbasis)
        }
    }
    
    if (mode == "latent") {
        # Create constant columns for means of scale and nominal parts
        J = matrix(1, nrow = nrow(X))
        nomm = rescale[2] * apply(bigNom, 2, mean)
        X = rescale[2] * X
        if (!is.null(S)) {
            sm = apply(S, 2, mean)
            X = cbind(X, kronecker(-J, matrix(sm, nrow = 1)))
        }
        bigX = cbind(kronecker(-J, matrix(nomm, nrow = 1)), X)
        misc$offset.mult = misc$offset.mult * rescale[2]
        intcpt = seq_len(ncol(tJac))
        bhat[intcpt] = bhat[intcpt] - rescale[1] / rescale[2]
    }
    else { ### ----- Piece together big matrix for each threshold ----- ###
        misc$ylevs = list(cut = cnm)
        misc$tran = link
        misc$inv.lbl = "cumprob"
        misc$offset.mult = -1
        if (!is.null(S))
            X = cbind(X, S)
        J = matrix(1, nrow=nrow(tJac))
        bigX = cbind(bigNom, kronecker(-J, X))
        if (mode != "linear.predictor") {
            misc$mode = mode
            misc$respName = as.character.default(object$terms)[2]
            misc$postGridHook = ".clm.postGrid"
        }
    }
    
    dimnames(bigX)[[2]] = names(bhat)
    
    list(X = bigX, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
         dfargs = list(), misc = misc)
}

# fuction called at end of ref_grid
# I use this for polr as well
# Also used for stanreg result of stan_polr & potentially other MCMC ordinal models
.clm.postGrid = function(object, ...) {
    mode = object@misc$mode
    object@misc$postGridHook = object@misc$mode = NULL
    object = regrid(object, transform = "response", ...)
    if(object@misc$estName == "exc.prob") { # back-transforming yields exceedance probs
        object@bhat = 1 - object@bhat
        if(!is.null(object@post.beta[1]))
            object@post.beta = 1 - object@post.beta
        object@misc$estName = "cum.prob"
    }
    if (mode == "prob") {
        object = .clm.prob.grid(object, ...)
    }
    else if (mode == "mean.class") {
        object = .clm.mean.class(object, ...)
    }
    else if (mode == "exc.prob") {
        object@bhat = 1 - object@bhat
        if(!is.null(object@post.beta[1]))
            object@post.beta = 1 - object@post.beta
        object@misc$estName = "exc.prob"        
    }
    # (else mode == "cum.prob" and it's all OK)
    object@misc$respName = NULL # cleanup
    object
}


# Make the linear-predictor ref_grid into one for class probabilities
# This assumes that object has already been re-gridded and back-transformed
.clm.prob.grid = function(object, thresh = "cut", newname = object@misc$respName, ...) {
    byv = setdiff(names(object@levels), thresh)
    newrg = contrast(object, ".diff_cum", by = byv, ...)
    if (!is.null(wgt <- object@grid[[".wgt."]])) {
        km1 = length(object@levels[[thresh]])
        wgt = wgt[seq_len(length(wgt) / km1)] # unique weights for byv combs
        newrg@grid[[".wgt."]] = rep(wgt, each = km1 + 1)
    }
    # proceed to disavow that this was ever exposed to 'emmeans' or 'contrast'
    ## class(newrg) = "ref.grid"
    misc = newrg@misc
    misc$is.new.rg = TRUE
    misc$infer = c(FALSE,FALSE)
    misc$estName = "prob"
    misc$pri.vars = misc$by.vars = misc$con.coef = misc$orig.grid = NULL
    newrg@misc = misc
    names(newrg@levels)[1] = names(newrg@grid)[1] = newname
    newrg@roles = object@roles
    newrg@roles$multresp = newname
    newrg
}

.clm.mean.class = function(object, ...) {
    prg = .clm.prob.grid(object, newname = "class", ...)
    byv = setdiff(names(prg@levels), "class")
    lf = as.numeric(prg@levels$class)
    newrg = contrast(prg, list(mean = lf), by = byv, ...)
    newrg = update(newrg, infer = c(FALSE, FALSE), 
        pri.vars = NULL, by.vars = NULL, estName = "mean.class")
    newrg@levels$contrast = newrg@grid$contrast = NULL
    prg@roles$multresp = NULL
    newrg@roles = prg@roles
    ## class(newrg) = "ref.grid"
    update(newrg, is.new.rg = TRUE)
}

# Contrast fcn for turning estimates of cumulative probabilities
# into cell probabilities
.diff_cum.emmc = function(levs, sep = "|", ...) {
    plevs = unique(setdiff(unlist(strsplit(levs, sep, TRUE)), sep))
    k = 1 + length(levs)
    if (length(plevs) != k)
        plevs = seq_len(k)
    M = matrix(0, nrow = length(levs), ncol = k)
    for (i in seq_along(levs))
        M[i, c(i,i+1)] = c(1,-1)
    dimnames(M) = list(levs, plevs)
    M = as.data.frame(M)
    attr(M, "desc") = "Differences of cumulative probabilities"
    attr(M, "adjust") = "none"
    attr(M, "offset") = c(rep(0, k-1), 1)
    M
}

#### replacement estimation routines for cases with a scale param

## workhorse for estHook and vcovHook functions
.clm.hook = function(object, tol = 1e-8, ...) {
    scols = object@misc$scale.idx
    bhat = object@bhat
    active = !is.na(bhat)
    bhat[!active] = 0
    linfct = object@linfct
    estble = estimability::is.estble(linfct, object@nbasis, tol) ###apply(linfct, 1, .is.estble, object@nbasis, tol)
    estble[!estble] = NA
    rsigma = estble * as.numeric(linfct[, scols, drop = FALSE] %*% object@bhat[scols])
    rsigma = exp(rsigma) * estble
    # I'll do the scaling later
    eta = as.numeric(linfct[, -scols, drop = FALSE] %*% bhat[-scols])
    if (!is.null(object@grid$.offset.))
        eta = eta + object@grid$.offset.
    for (j in scols) linfct[, j] = eta * linfct[, j]
    linfct = (.diag(rsigma) %*% linfct) [, active, drop = FALSE]
    list(est = eta * rsigma, V = linfct %*% tcrossprod(object@V, linfct))
}

.clm.estHook = function(object, do.se = TRUE, tol = 1e-8, ...) {
    raw.matl = .clm.hook(object, tol, ...)
    SE = if (do.se) sqrt(diag(raw.matl$V))  else NA
    cbind(est = raw.matl$est, SE = SE, df = Inf)
}

.clm.vcovHook = function(object, tol = 1e-8, ...) {
    .clm.hook(object, tol, ...)$V
}

### Special emm_basis fcn for the scale model
.emm_basis.clm.scale = function(object, trms, xlev, grid, ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$S.contrasts)
    bhat = c(`(intercept)` = 0, object$zeta)
    nbasis = estimability::all.estble
    if (any(is.na(bhat)))
        nbasis = estimability::nonest.basis(model.matrix(object)$S)
    k = sum(!is.na(bhat)) - 1
    V = .my.vcov(object, ...)
    pick = nrow(V) - k + seq_len(k)
    V = V[pick, pick, drop = FALSE]
    V = cbind(0, rbind(0,V))
    misc = list(tran = "log")
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, 
         dffun = function(...) Inf, dfargs = list(), misc = misc)
}

emm_basis.clmm = function (object, trms, xlev, grid, ...) {
    if(is.null(object$Hessian)) {
        message("Updating the model to obtain the Hessian...")
        object = update(object, Hess = TRUE)
    }
    # borrowed from Maxime's code -- need to understand this better, e.g. when it happens
    H = object$Hessian
    if (any(apply(object$Hessian, 1, function(x) all(x == 0)))) {
        H = H[names(coef(object)), names(coef(object))]
        object$Hessian = H
    }
    result = emm_basis.clm(object, trms, xlev, grid, ...)
    # strip off covariances of random effects
    keep = seq_along(result$bhat)
    result$V = result$V[keep,keep]
    result
}
