##############################################################################
#    Copyright (c) 2012-2020 Russell V. Lenth                                #
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

# lqmm and lqm support

#' @exportS3Method recover_data lqmm
recover_data.lqmm = function(object, data = object$mfArgs$data, ...) {
    fcall = object$call
    trms = delete.response(terms(eval(fcall$fixed)))
    recover_data(fcall, trms, object$mfArgs$na.action, data = data, ...)
}

#' @exportS3Method emm_basis lqmm         
emm_basis.lqmm = function(object, trms, xlev, grid, tau = 0.5, ...) {
    taudiff = abs(object$tau - tau)
    col = which(taudiff < 0.0001)
    if (length(col) == 0)
        stop("No coefficients available for tau = ", tau)
    bhat = coef(object)
    # Very touchy here because their boot() function doesn't take dots...
    nm = intersect(names(list(...)), c("method", "R", "seed", "startQR"))
    vargs = c(list(object = object, covariance = TRUE), list(...)[nm])
    V = do.call("summary", vargs)$Cov[names(bhat), names(bhat)]
    if (length(taudiff) > 1) {
        bhat = bhat[, col[1]]
        V = V[, , col]
    }
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    nbasis = estimability::all.estble
    dfargs = list(df = object$rdf)
    dffun = function(k, dfargs) dfargs$df
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, 
         dffun = dffun, dfargs = dfargs, misc = list())
}


# Use same functions for lqm objects
#' @exportS3Method recover_data lqm
recover_data.lqm = function(object, ...) {
    recover_data.lm(object, frame = NULL, ...)
}

#' @exportS3Method emm_basis lqm          
emm_basis.lqm = function(object, ...)
    emm_basis.lqmm(object, ...)



#### rq objects (quantreg)

#' @exportS3Method recover_data rq
recover_data.rq = function(object, ...) {
    recover_data.lm(object, frame = object$model, ...)
}
    
#' @exportS3Method emm_basis rq           
emm_basis.rq = function(object, trms, xlev, grid, tau = 0.5, ...) {
    taudiff = abs(object$tau - tau)
    col = which(taudiff < 0.0001)
    if (length(col) == 0)
        stop("No coefficients available for tau = ", tau)
    bhat = object$coefficients
    summ = summary(object, covariance = TRUE, ...)
    if (length(taudiff) == 1) {
        V = summ$cov
        df = summ$rdf
    }
    else {
        bhat = bhat[, col[1]]
        V = summ[[col]] $ cov
        df = summ[[col]] $ rdf
    }
    nm = if(is.null(names(bhat))) row.names(bhat) else names(bhat)
    m = suppressWarnings(model.frame(trms, grid, na.action = na.pass, xlev = xlev))
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    assign = attr(X, "assign")
    X = X[, nm, drop = FALSE]
    bhat = as.numeric(bhat) 
    nbasis = estimability::all.estble
    misc = list()
    dfargs = list(df = df)
    dffun = function(k, dfargs) dfargs$df
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, 
         dffun = dffun, dfargs = dfargs, misc = misc)
}

# we just reroute rqs objects to emm_basis.rq, as pretty similar
#' @exportS3Method recover_data rqs
recover_data.rqs = recover_data.rq

#' @exportS3Method emm_basis rqs          
emm_basis.rqs = function(object, ...)
    emm_basis.rq(object, ...)
