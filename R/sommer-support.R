##############################################################################
#    Copyright (c) 2012-2019 Russell V. Lenth                                #
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

# sommer package support

#' @exportS3Method recover_data mmer
recover_data.mmer = function(object, data, ...) {
    if (is.null(data))
        data = object$data
    fcall = call("mmer", formula = object$call$fixed, data = data)
    emmeans::recover_data(fcall, delete.response(terms(object$call$fixed)), 
                          object$call$na.method.V, ...)
}

#' @exportS3Method emm_basis mmer         
emm_basis.mmer = function(object, trms, xlev, grid, ...) {
    cf = object$Beta
    bhat = cf$Estimate
    m = suppressWarnings(model.frame(trms, grid, na.action = na.pass, xlev = xlev))
    # if we can get contrasts from the object, fix next line
    X = model.matrix(trms, m, contrasts.arg = NULL)
    V = .my.vcov(object, vcov. = function(., ...) .$VarBeta)
    
    if ((k <- length(bhat)) < ncol(X)) {  # we have rank deficiencies
        QR = qr(model.matrix(trms, object$data, contrasts.arg = NULL))
        bhat = rep(NA, ncol(X))
        bhat[QR$pivot[1:k]] = cf$Estimate
        nbasis = estimability::nonest.basis(QR)
    }
    else
        nbasis = estimability::all.estble 
    misc = list()
    # soup-up following if (1) glms allowed or (2) d.f. available
    dfargs = list(df = object$df.residual)
    dffun = function(k, dfargs) Inf
    bas = list(X = X, bhat = bhat, nbasis = nbasis, V = V, 
               dffun = dffun, dfargs = dfargs, misc = misc)
    # check for multiv resp
    k = length(levels(cf$Trait))
    if (k > 1) {
        bas$misc$ylevs = list(Trait = levels(cf$Trait))
        bas$X = kronecker(diag(rep(1, k)), bas$X)
        # reorder coefs to go one trait at a time
        ord = as.integer(matrix(seq_along(bas$bhat), ncol = k, byrow = TRUE))
        bas$bhat = bas$bhat[ord]
        bas$V = bas$V[ord, ord, drop = FALSE]
    }
    bas
}