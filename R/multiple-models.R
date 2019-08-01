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

### MuMIn::averaging
### This is just bare-bones. Need to provide data,
### and if applicable, df.
### Optional argument tran sets the tran property

# $data is NOT a standard member, but if it's there, we'll use it
# Otherwise, we need to provide data or its name in the call
recover_data.averaging = function(object, data, ...) {
    ml = attr(object, "modelList")
    if (is.null(ml))
        return(paste0("emmeans support for 'averaging' models requires a 'modelList' attribute.\n",
                      "Re-fit the model from a model list or with fit = TRUE"))
    if (is.null(object$formula)) {
        lhs = as.formula(paste(formula(ml[[1]])[[2]], "~."))
        rhs = sapply(ml, function(m) {f = formula(m); f[[length(f)]]})
        object$formula = update(as.formula(paste("~", paste(rhs, collapse = "+"))), lhs)
    }
    trms = terms(object$formula)
    if (is.null(data))
        data = ml[[1]]$call$data
    fcall = call("model.avg", formula = object$formula, data = data)
    recover_data(fcall, delete.response(trms), na.action = NULL, ...)
}

emm_basis.averaging = function(object, trms, xlev, grid, ...) {
    bhat = coef(object, full = TRUE)
    V = .my.vcov(object, function(.) vcov(., full = TRUE))
    m = suppressWarnings(model.frame(trms, grid, na.action = na.pass, xlev = xlev))
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)

    nbasis = estimability::all.estble
    ml = attr(object, "modelList")
    ml1 = ml[[1]]
    misc = list()
    if (!is.null(fam <- family(ml1)))
        misc = .std.link.labels(fam, misc)
    if (!is.null(df.residual(ml1))) {
        dffun = function(k, dfargs) dfargs$df
        dfargs = list(df = min(sapply(ml, df.residual)))
    }
    else {
        dffun = function(k, dfargs) Inf
        dfargs = list()
    }
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, 
         dffun = dffun, dfargs = dfargs, misc = misc)
}




### minc::mira support -----------------------------------------
# Here we rely on the methods already in place for elements of $analyses

recover_data.mira = function(object, ...) {
    rdlist = lapply(object$analyses, recover_data, ...)
    rd = rdlist[[1]]
    # we'll average the numeric columns...
    numcols = which(sapply(rd, is.numeric))
    for (j in numcols)
        rd[, j] = apply(sapply(rdlist, function(.) .[, j]), 1, mean)
    rd
}

emm_basis.mira = function(object, trms, xlev, grid, ...) {
    bas = emm_basis(object$analyses[[1]], trms, xlev, grid, ...)
    k = length(object$analyses)
    # we just average the V and bhat elements...
    V = 1/k * bas$V
    allb = cbind(bas$bhat, matrix(0, nrow = length(bas$bhat), ncol = k - 1))
    for (i in 1 + seq_len(k - 1)) {
        basi = emm_basis(object$analyses[[i]], trms, xlev, grid, ...)
        V = V + 1/k * basi$V
        allb[, i] = basi$bhat
    }
    bas$bhat = apply(allb, 1, mean)  # averaged coef
    notna = which(!is.na(bas$bhat))
    bas$V = V + (k + 1)/k * cov(t(allb[notna, , drop = FALSE]))   # pooled via Rubin's rules
    bas
}
