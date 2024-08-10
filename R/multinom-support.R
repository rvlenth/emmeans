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
#' @exportS3Method recover_data multinom
recover_data.multinom = function(object, ...) {
    fcall = object$call
    recover_data(fcall, delete.response(terms(object)), object$na.action, ...)
}

#' @exportS3Method emm_basis multinom     
emm_basis.multinom = function(object, trms, xlev, grid, 
                              mode = c("prob", "latent"), ...) {
    mode = match.arg(mode)
    bhat = t(coef(object))
    V = .my.vcov(object, ...)
    # NOTE: entries in vcov(object) come out in same order as
    # in as.numeric(bhat), even though latter has been transposed
    k = ifelse(is.matrix(coef(object)), ncol(bhat), 1)
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    # recenter for latent predictions
    pat = (rbind(0, diag(k + 1, k)) - 1) / (k + 1)
    X = kronecker(pat, X)
    nbasis = estimability::all.estble
    nbasis = kronecker(rep(1,k), nbasis)
    misc = list(tran = "clr")   ### misc = list(tran = "log", inv.lbl = "e^y")
    dfargs = list(df = object$edf)
    dffun = function(k, dfargs) dfargs$df
    if(is.null(ylevs <- object$lev))
        ylevs = seq_len(k + 1)
    ylevs = list(class = ylevs)
    if (is.null(ylevs)) ylevs = list(class = seq_len(k))
    names(ylevs) = as.character.default(eval(object$call$formula, environment(trms))[[2]])[1]
    misc$ylevs = ylevs
    if (mode == "prob")
        misc$postGridHook = .multinom.postGrid
    list(X = X, bhat = as.numeric(bhat), nbasis = nbasis, V = V, 
         dffun = dffun, dfargs = dfargs, misc = misc)
}

# post-processing of ref_grid for "prob" mode
# also allows simulated outcomes
## Note - now that we have mvregrid, we could just pass the latent grid
## through with the clrInv transform. But this works so I'm just leaving it here
.multinom.postGrid = function(object, N.sim, ...) {
    linfct = object@linfct
    misc = object@misc
    # grid will have multresp as slowest-varying factor...
    idx = matrix(seq_along(linfct[, 1]), 
                 ncol = length(object@levels[[object@roles$multresp]]))
    bhat = as.numeric(idx) # right length, contents will be replaced
    if(sim <- !missing(N.sim)) {
        message("Simulating a sample of size ", N.sim)
        bsamp = mvtnorm::rmvnorm(N.sim, object@bhat, object@V)
        postb = matrix(0, nrow = N.sim, ncol = length(bhat))
    }
    for (i in 1:nrow(idx)) {
        rows = idx[i, ]
        exp.psi = exp(linfct[rows, , drop = FALSE] %*% object@bhat)
        p = as.numeric(exp.psi / sum(exp.psi))
        bhat[rows] = p
        if (sim) {
            ex = exp(linfct[rows, , drop = FALSE] %*% t(bsamp))  # p x N
            px = t(apply(ex, 2, function(x) x / sum(x)))
            postb[, rows] = px
        }
        A = .diag(p) - outer(p, p)    # partial derivs
        linfct[rows, ] = A %*% linfct[rows, ]
    }
    misc$postGridHook = misc$tran = misc$inv.lbl = NULL
    misc$estName = "prob"
    
    object@bhat = bhat
    object@V = linfct %*% tcrossprod(object@V, linfct)
    object@linfct = diag(1, length(bhat))
    object@misc = misc
    if (sim)
        object@post.beta = postb
    object
}


### Support for mclogit::mblogit models???
#' @exportS3Method recover_data mblogit
recover_data.mblogit = function (object, ...) 
{
    recover_data.multinom(object, ...)
}

#' @exportS3Method emm_basis mblogit      
emm_basis.mblogit = function(object, ..., vcov.) {
    object$coefficients = object$coefmat
    object$lev = levels(object$model[[1]])
    object$edf = Inf
    # we have to arrange the vcov elements in row-major order
    if(missing(vcov.))
        vcov. = vcov(object)
    perm = matrix(seq_along(as.numeric(object$coefmat)), 
                  ncol = ncol(object$coefmat))
    perm = as.numeric(t(perm))
    vcov. = vcov.[perm, perm]
    emm_basis.multinom(object, ..., vcov. = vcov.)
}

