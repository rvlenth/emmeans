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

### Multinomial modeling

### Example for testing
### From: http://www.ats.ucla.edu/stat/r/dae/mlogit.htm
# library(foreign)
# ml <- read.dta("http://www.ats.ucla.edu/stat/data/hsbdemo.dta")
# library(nnet)
# ml$prog2 <- relevel(ml$prog, ref = "academic")
# test <- multinom(prog2 ~ ses + write, data = ml)
# 

# same as recover_data.lm
recover_data.multinom = function(object, ...) {
    fcall = object$call
    recover_data(fcall, delete.response(terms(object)), object$na.action, ...)
}

emm_basis.multinom = function(object, trms, xlev, grid, 
                              mode = c("prob", "latent"), ...) {
    mode = match.arg(mode)
    bhat = t(coef(object))
    V = .my.vcov(object, ...)
    k = ifelse(is.matrix(coef(object)), ncol(bhat), 1)
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    # recenter for latent predictions
    pat = (rbind(0, diag(k + 1, k)) - 1) / (k + 1)
    X = kronecker(pat, X)
    nbasis = estimability::all.estble
    nbasis = kronecker(rep(1,k), nbasis)
    misc = list(tran = "log", inv.lbl = "e^y")
    dfargs = list(df = object$edf)
    dffun = function(k, dfargs) dfargs$df
    ylevs = list(class = object$lev)
    if (is.null(ylevs)) ylevs = list(class = seq_len(k))
    names(ylevs) = as.character.default(object$call$formula[[2]])
    misc$ylevs = ylevs
    if (mode == "prob")
        misc$postGridHook = .multinom.postGrid
    list(X = X, bhat = as.numeric(bhat), nbasis = nbasis, V = V, 
         dffun = dffun, dfargs = dfargs, misc = misc)
}

# post-processing of ref_grid for "prob" mode
.multinom.postGrid = function(object, ...) {
    linfct = object@linfct
    misc = object@misc
    # grid will have multresp as slowest-varying factor...
    idx = matrix(seq_along(linfct[, 1]), 
                 ncol = length(object@levels[[object@roles$multresp]]))
    bhat = as.numeric(idx) # right length, contents will be replaced
    for (i in 1:nrow(idx)) {
        rows = idx[i, ]
        exp.psi = exp(linfct[rows, , drop = FALSE] %*% object@bhat)
        p = as.numeric(exp.psi / sum(exp.psi))
        bhat[rows] = p
        A = .diag(p) - outer(p, p)    # partial derivs
        linfct[rows, ] = A %*% linfct[rows, ]
    }
    misc$postGridHook = misc$tran = misc$inv.lbl = NULL
    misc$estName = "prob"
    
    object@bhat = bhat
    object@V = linfct %*% tcrossprod(object@V, linfct)
    object@linfct = diag(1, length(bhat))
    object@misc = misc
    object
}