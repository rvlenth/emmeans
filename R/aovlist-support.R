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

# emmeans support for aovlist objects

#' @export
recover_data.aovlist = function(object, ...) {
    fcall = match.call(aov, attr(object, "call")) # matches even if called via `manova()`
    trms = terms(object)
    # Find the Error terms
    lbls = attr(trms, "term.labels")
    err.idx = grep("^Error\\(", lbls)
    newf = as.formula(paste(c(".~.", lbls[err.idx]), collapse = "-"))
    trms = terms(update(trms, newf))
    dat = recover_data(fcall, delete.response(trms), na.action = attr(object, "na.action"), ...)
    attr(dat, "pass.it.on") = TRUE
    dat
}

# This works great for balanced experiments, and goes horribly wrong
# even for slightly unbalanced ones. So I abort on these kinds of cases
#' @export
emm_basis.aovlist = function (object, trms, xlev, grid, vcov., ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    contr = attr(object, "contrasts")
    non.orth = sapply(contr, function(x) {
        !is.character(x) || !(x %in% c("contr.helmert", "contr.poly", "contr.sum"))
    })
    if (any(non.orth)) { # refit with sum-to-zero contrasts
        contr[non.orth] = "contr.sum"
        message("Note: re-fitting model with sum-to-zero contrasts")
        cl = attr(object, "call")
        cl$contrasts = contr
        object = eval(cl)
    }
    X = model.matrix(trms, m, contrasts.arg = contr)
    xnms = dimnames(X)[[2]]
    
    # Check for situations we can't handle...
    colsums = apply(X[, setdiff(xnms, "(Intercept)"), drop=FALSE], 2, sum)
    if (any(round(colsums,3) != 0))
        warning("Some predictors are correlated with the intercept - results may be very biased")
    if (length(unlist(lapply(object, function(x) names(coef(x))))) > length(xnms))
        message("NOTE: Results are based on intra-block estimates and are biased.")
    
    # initialize arrays
    nonint = setdiff(names(object), "(Intercept)")
    
    k = npar = length(xnms)
    bhat1 = rep(NA, k) # I'll use NAs in 1st dim of bhat to track which slots I've filled
    # check for multivariate response
    m = ifelse (is.matrix(coefm <- object[[1]]$coefficients), ncol(coefm), 1)
    bhat = rep(NA, k*m)
    zmm1 = seq_len(m) - 1   # seq 0 : (m-1)
    
    
    #### utility functions...
    # return indices of as.numeric(x[i, ]) where x has nr rows
    indx = function(i, nr = npar) 
        as.numeric(sapply(zmm1, function(j) nr * j + i))
    
    # get names or rownames
    mynames = function(x)
        if (m == 1) names(x) else rownames(x)
    
    V = matrix(0, nrow = k*m, ncol = k*m)
    names(bhat1) = xnms
    allxnms = xnms
    if (m > 1) {
        ylevs = colnames(coefm)
        allxnms = as.character(sapply(ylevs, function(.) paste0(xnms, .)))
    }
    names(bhat) = allxnms
    dimnames(V) = list(allxnms, allxnms)
    empty.list = as.list(nonint)
    names(empty.list) = nonint
    Vmats = Vidx = Vdf = empty.list
    wts = matrix(0, nrow = length(nonint), ncol = k*m)
    dimnames(wts) = list(nonint, allxnms)
    # NOTE: At present, I just do intra-block analysis: wts are all 0 and 1
    btemp = bhat1 #++ temp for tracking indexes
    #++Work thru strata in reverse order
    for (nm in rev(nonint)) {
        x = object[[nm]]
        if (m > 1) class(x) = c("mlm", "lm")   # because vcov.aov is NOT suitable
        bi = coef(x)
        rn = mynames(bi)
        nr = length(rn)
        idx = which(!is.na(bi[seq_along(rn)]))
        bi = bi[indx(idx, nr)]
        ii = match(rn[idx], xnms)
        use = setdiff(ii, which(!is.na(bhat1))) #++ omit elts already filled
        if(length(use) > 0) {
            ii.left = seq_along(ii)[!is.na(match(ii,use))]
            wts[nm, indx(use)] = 1
            bhat1[use] = bi[ii.left]
            allii.left = indx(ii.left, nr)
            alluse = Vidx[[nm]] = indx(use)
            bhat[alluse] = bi[allii.left]
            # following is OK now that we have class(x) = "mlm"
            Vi = vcov(x, complete = FALSE)[allii.left, allii.left, drop = FALSE]
            Vmats[[nm]] = Vi
            V[alluse, alluse] = Vi
        }
        else {
            Vmats[[nm]] = matrix(0, nrow=0, ncol=0)
            Vidx[[nm]] = integer(0)
        }
        # Any cases with 0 df will have NaN for covariances. I make df = -1 
        # in those cases so I don't divide by 0 later in Satterthwaite calcs
        Vdf[[nm]] = ifelse(x$df > 0, x$df, -1)
    }
    
    x <- object[["(Intercept)"]]
    if (!is.null(x)) {
        # The intercept belongs in the 1st error stratum
        # So have to add a row and column to its covariance matrix
        idx1 = indx(1)
        bhat[idx1] = as.numeric(coef(x))
        wts[1, idx1] = 1
        x = object[[nonint[1]]] # 1st non-intercept stratum
        Vidx[[1]] = ii = sort(c(idx1, Vidx[[1]]))
        k = length(ii)
        vv = matrix(0, nrow = k, ncol = k)
        i2k = indx(2:(k / m), k / m)
        if (k > m) vv[i2k, i2k] = Vmats[[1]]
        # Variance of intercept is EMS of this stratum divided by N
        # Here I'm assuming there are no weights
        N = sum(sapply(object, function(x) length(x$residuals))) / m
        i1 = indx(1, k / m)
        if (m > 1)
            V[idx1, idx1] = vv[i1,i1] = estVar(x) / N
        else
            V[1,1] = vv[1,1] = sum(resid(x)^2) / x$df / N
        #dimnames(vv) = list(c(xnms[ii], xnms[ii]))
        Vmats[[1]] = vv
    }
    # override V if vcov. is supplied
    if(!missing(vcov.)) {
        V = .my.vcov(object, vcov.)
        dfargs = list()
        dffun = function(k, dfargs) Inf
    }
    else {
        dfargs = list(Vmats=Vmats, Vidx=Vidx, Vdf=unlist(Vdf), wts = wts)
        dffun = function(k, dfargs) {
            emmeans::.aovlist.dffun(k, dfargs)
        }
    }
    nbasis = estimability::all.estble  # Consider this further?
    misc = list(initMesg = "Warning: EMMs are biased unless design is perfectly balanced")
    if (m > 1) {
        misc$ylevs = list(rep.meas = ylevs)
        X = kronecker(diag(1, m), X)
    }
    
    # submodel support
    mm = NULL
    if(!is.null(dat <- attr(object, "data"))) {
        m = model.frame(trms, dat, na.action = na.pass, xlev = xlev)
        mm = model.matrix(trms, m, contrasts.arg = contr)
        mm = .cmpMM(mm, assign = attr(mm, "assign"))
    }
    
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
         dfargs = dfargs, misc = misc, model.matrix = mm)
}

#' @rdname extending-emmeans
#' @order 32
#' @param k,dfargs Arguments to \code{.aovlist.dffun}, which is made available as a 
#'   convenience to developers providing support similar to that provided for 
#'   \code{aovlist} objects
#' @export
.aovlist.dffun = function(k, dfargs) {
    if(is.matrix(k) && (nrow(k) > 1)) {
        dfs = apply(k, 1, .aovlist.dffun, dfargs)
        min(dfs)
    }
    else {
        v = sapply(seq_along(dfargs$Vdf), function(j) {
            ii = dfargs$Vidx[[j]]
            kk = (k * dfargs$wts[j, ])[ii]            
            #sum(kk * .mat.times.vec(dfargs$Vmats[[j]], kk))
            .qf.non0(dfargs$Vmats[[j]], kk)
        })
        sum(v)^2 / sum(v^2 / dfargs$Vdf) # Good ole Satterthwaite
    }
}