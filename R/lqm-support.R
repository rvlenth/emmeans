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

recover_data.lqmm = function(object, data = object$mfArgs$data, ...) {
    fcall = object$call
    trms = delete.response(terms(eval(fcall$fixed)))
    recover_data(fcall, trms, object$mfArgs$na.action, data = data, ...)
}

emm_basis.lqmm = function(object, trms, xlev, grid, tau = 0.5, ...) {
    bhat = coef(object)
    col = which(abs(tau - as.numeric(colnames(bhat))) < 0.0001)
    if (length(col) == 0)
        stop("No coefficients available for tau = ", tau)
    col = col[1]
    bhat = bhat[, col]
    nm = names(bhat)
    vcv = summary(object, covariance = TRUE)$Cov
    V = vcv[nm, nm, col]
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    nbasis = estimability::all.estble
    dfargs = list(df = object$rdf)
    dffun = function(k, dfargs) dfargs$df
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, 
         dffun = dffun, dfargs = dfargs, misc = list())
}


# Use same functions for lqm objects
recover_data.lqm = function(object, ...) {
    recover_data.lm(object, frame = NULL, ...)
}

emm_basis.lqm = function(object, ...)
    emm_basis.lqmm(object, ...)
