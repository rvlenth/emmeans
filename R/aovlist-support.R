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
    fcall = attr(object, "call")
    trms = terms(object)
    # Find the Error terms
    lbls = attr(trms, "term.labels")
    err.idx = grep("^Error\\(", lbls)
    newf = as.formula(paste(c(".~.", lbls[err.idx]), collapse = "-"))
    trms = terms(update(trms, newf))
    recover_data(fcall, delete.response(trms), na.action = attr(object, "na.action"), ...)
}

# This works great for balanced experiments, and goes horribly wrong
# even for slightly unbalanced ones. So I abort on these kinds of cases
#' @export
emm_basis.aovlist = function (object, trms, xlev, grid, vcov., ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    contr = attr(object, "contrasts")
    X = model.matrix(trms, m, contrasts.arg = contr)
    xnms = dimnames(X)[[2]]
    
    # Check for situations we can't handle...
    colsums = apply(X[, setdiff(xnms, "(Intercept)"), drop=FALSE], 2, sum)
    if (any(round(colsums,3) != 0))
        warning("Some predictors are correlated with the intercept - results are biased.\n",
                "May help to re-fit with different contrasts, e.g. 'contr.sum'")
    if (length(unlist(lapply(object, function(x) names(coef(x))))) > length(xnms))
        message("NOTE: Results are based on intra-block estimates.")
    
    # initialize arrays
    nonint = setdiff(names(object), "(Intercept)")
    
    k = length(xnms)
    bhat = rep(NA, k) # I'll use NAs to track which slots I've filled
    V = matrix(0, nrow=k, ncol=k)
    names(bhat) = xnms
    dimnames(V) = list(xnms, xnms)
    empty.list = as.list(nonint)
    names(empty.list) = nonint
    Vmats = Vidx = Vdf = empty.list
    wts = matrix(0, nrow = length(nonint), ncol = k)
    dimnames(wts) = list(nonint, xnms)
    # NOTE: At present, I just do intra-block analysis: wts are all 0 and 1
    btemp = bhat #++ temp for tracking indexes
    #++Work thru strata in reverse order
    for (nm in rev(nonint)) {
        x = object[[nm]]
        bi = coef(x)
        bi = bi[!is.na(bi)]
        ii = match(names(bi), xnms)
        Vidx[[nm]] = use = setdiff(ii, which(!is.na(bhat))) #++ omit elts already filled
        if(length(use) > 0) {
            ii.left = seq_along(ii)[!is.na(match(ii,use))]
            wts[nm, use] = 1
            bhat[use] = bi[ii.left]
            Vi = vcov(x, complete = FALSE)[ii.left, ii.left, drop=FALSE]
            Vmats[[nm]] = Vi
            V[use,use] = Vi
        }
        else {
            Vmats[[nm]] = matrix(0, nrow=0, ncol=0)
        }
        # Any cases with 0 df will have NaN for covariances. I make df = -1 
        # in those cases so I don't divide by 0 later in Satterthwaite calcs
        Vdf[[nm]] = ifelse(x$df > 0, x$df, -1)
    }
    
    x <- object[["(Intercept)"]]
    if (!is.null(x)) {
        # The intercept belongs in the 1st error stratum
        # So have to add a row and column to its covariance matrix
        bhat[1] = x$coefficients[1]
        wts[1,1] = 1
        Vidx[[1]] = ii = c(1, Vidx[[1]])
        k = length(ii)
        vv = matrix(0, nrow=k, ncol=k)
        if (k > 1) vv[2:k,2:k] = Vmats[[1]]
        # Variance of intercept is EMS of this stratum divided by N
        # Here I'm assuming there are no weights
        N = sum(sapply(object, function(x) length(x$residuals)))
        V[1,1] = vv[1,1] = sum(object[[2]]$residuals^2) / object[[2]]$df / N
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
    misc = list()
    
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, dffun = dffun, 
         dfargs = dfargs, misc = misc)
}

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