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

### Support for objects involving several models (e.g., model averaging, multiple imputation)

### MuMIn:: averaging
### This is just bare-bones. Need to provide data,
### and if applicable, df.
### Optional argument tran sets the tran property

# $data is NOT a standard member, but if it's there, we'll use it
# Otherwise, we need to provide data or its name in the call
recover_data.averaging = function(object, data, ...) {
    trms = terms(object$formula)
    if (is.null(data))
        data = eval(object$data)
    fcall = call("model.avg", formula = object$formula, data = data)
    recover_data(fcall, delete.response(trms), na.action = NULL, ...)
}

emm_basis.averaging = function(object, trms, xlev, grid, ...) {
    bhat = coef(object, full = TRUE)
    V = .my.vcov(object, function(.) vcov(., full = TRUE))
    m = suppressWarnings(model.frame(trms, grid, na.action = na.pass, xlev = xlev))
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)

    nbasis = estimability::all.estble
    misc = list()
    dffun = function(k, dfargs) Inf
    dfargs = list()
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, 
         dffun = dffun, dfargs = dfargs, misc = misc)
}




### On second thought, I don't think following is very useful
# # Average together the coefs (weights defauult to equal)
# # Columns of bhatmat corresp to defferent solutions
# # setNA0 determines whether we take NAs as solutions constrained to 0 or left as NA
# #  so in return value, we get NA in these posituions:
# #    setNA0 = TRUE:  those for which all ests are NA
# #    setNA0 = FALSE: those for which any est is NA
# .coef.avg = function(bhatmat, weights, setNA0 = TRUE) {
#     k = ncol(bhatmat)
#     if (missing(weights))
#         weights = rep(1/k, k)
#     nzw = which(weights > 0)
#     bhatmat = bhatmat[ , nzw, drop = FALSE]
#     weights = weights[nzw]
#     if (setNA0) {
#         nas = which(apply(bhatmat, 1, function(.) all(is.na(.))))
#         bhatmat[is.na(bhatmat)] = 0
#     }
#     else
#         nas = integer(0)
#     wbhat = apply(bhatmat, 1, weighted.mean, weights)
#     wbhat[nas] = NA
#     wbhat
# }
# 
# # Average list of vcovs: weighted average of Vi + (bi - wbi)(bi - wbi)'
# # Caller is responsible for ensuring that all elements of vcovlist conform to bhatmat etc.
# # Better to call with wbhat already defined (with same weights an bhatmat), but if missing, we compute it
# .vcov.avg = function(vcovlist, bhatmat, wbhat, weights, setNA0 = TRUE) {
#     k = ncol(bhatmat)
#     if (missing(weights))
#         weights = rep(1/k, k)
#     nzw = which(weights > 0)
#     bhatmat = bhatmat[ , nzw, drop = FALSE]
#     vcovlist = vcovlist[nzw]
#     weights = weights[nzw]
#     if (missing(wbhat))
#         wbhat = .avg.coef(wbhat, weights, setNA0 = setNA0)
#     notna = wbhat[!is.na(wbhat)]
#     xvmat = sapply(seq_along(vcovlist), function(i) {
#         bdev = bhatmat[notna, i] - wbhat[notna]
#         as.numeric(vcovmat[[i]][notna, notna, drop = FALSE] + outer(bdev, bdev))
#     })
#     vavg = matrix(apply(xvmat, 1, weighted.mean, weights), ncol = length(notna))
#     nms = colnames(bhatmat)
#     dimnames(vavg) = list(nms, nms)
#     vavg
# }
