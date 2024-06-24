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
#' @exportS3Method recover_data averaging
recover_data.averaging = function(object, data, formula = object$formula, ...) {
    ml = attr(object, "modelList")
    if (is.null(ml))
        return(paste0("emmeans support for 'averaging' models requires a 'modelList' attribute.\n",
                      "Re-fit the model from a model list or with fit = TRUE"))
    if (is.null(formula)) {
        lhs = as.formula(paste(formula(ml[[1]])[[2]], "~."))
        rhs = sapply(ml, function(m) {f = formula(m); f[[length(f)]]})
        formula = update(as.formula(paste("~", paste(rhs, collapse = "+"))), lhs)
    }
    if (is.null(data))
        data = ml[[1]]$call$data
    ### temporary patch -- $formula often excludes the intercept
    formula = update(formula, .~.+1) ####[-2]
    trms = attr(model.frame(formula, data = data), "terms")
    fcall = call("model.avg", formula = formula, data = data)
    recover_data(fcall, trms, na.action = NULL, ...)
}

# optional argument 'subset' is named character "prefix:pfx" or "wrap:wrap"
# or NAMED vector of indexes
#' @exportS3Method emm_basis averaging
emm_basis.averaging = function(object, trms, xlev, grid, subset, ...) {
    bhat = coef(object, full = TRUE)
    bnms = names(bhat)
    V = .my.vcov(object, function(., ...) vcov(., full = TRUE), ...)

    # support functions for 'what' ...
    prefixed_by = function(pfx) {  # picks elements prefixed by pfx
        w = which(startsWith(bnms, pfx))
        names(w) = gsub(pfx, "", bnms[w])
        w
    }
    wrapped_by = function(wrp) {  # picks elements wrapped by wrp(name)
        pat = paste0(wrp, "\\(")
        w = grep(pat, bnms)
        tmp = sub(pat, "", bnms[w])
        names(w) = substr(tmp, 1, nchar(tmp) - 1)
        ### names(w) = sub("([a-z]+\\()([0-9:_() A-Za-z]+)(\\))", "\\2", bnms[w], perl = TRUE)
        w
    }
    if(!missing(subset)) {
        if (is.character(subset)) {  # wrap:stg or prefix:stg
            tmp = strsplit(subset, ":")[[1]]
            if(tmp[[1]] == "wrap")
                subset = wrapped_by(tmp[2])
            else # anything else is interpreted as prefix
                subset = prefixed_by(rev(tmp)[1])
            if(length(subset) == 0)
                stop("subset specification eliminates all coefficients!")
        }
        bhat = bhat[subset]
        bnms = names(subset)
        bnms[bnms == "(Int)"] = "(Intercept)"  # covers a peculiarity in glmmTMB
        V = V[subset, subset]
    }
    
    # we're gonna assemble all the main-effects model-matrix components
    # and build the model matrix to match the names of the coefs
    fcts = rownames(attr(trms, "factors"))
    mmlist = lapply(fcts, \(f) {
        frm = reformulate(c(-1, f))
        suppressWarnings(model.matrix(frm, data = grid, xlev = xlev))
    })
    ME = cbind("(Intercept)" = 1, do.call(cbind, mmlist))
    cols = strsplit(bnms, ":")  # the cols we need
    ext = setdiff(intersect(unlist(cols), names(grid)), colnames(ME)) # extra expressions we need to compute
    if(length(ext) > 0) {
        xcols = sapply(ext, \(expr) eval(parse(text = expr), envir = grid))
        ME = cbind(ME, xcols)
    }
    # now get the factors for each coef
    X = sapply(cols, \(nms)
        apply(ME[, nms, drop = FALSE], 1, prod))
    colnames(X) = bnms
    # Now, sometimes a coefficient is duplicated, so we need to zero those out
    srt = sapply(cols, \(nms) paste(sort(nms), collapse = ":"))
    tmp = table(srt)
    for (nm in names(tmp[tmp > 1])) {  # work on the dupes
        zap = which(srt == nm)[-1]   # indices of all but the first
        bhat[zap] = 0
    }
    
    
    nbasis = estimability::all.estble
    ml = attr(object, "modelList")
    ml1 = ml[[1]]
    misc = list()
    fam <- try(family(ml1), silent = TRUE)
    if (!inherits(fam, "try-error"))
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

#' @exportS3Method recover_data mira
recover_data.mira = function(object, data = NULL, ...) {
    rdlist = lapply(object$analyses, recover_data, data = data, ...)
    rd = rdlist[[1]]
    # we'll average the numeric columns...
    numcols = which(sapply(rd, is.numeric))
    for (j in numcols)
        rd[, j] = apply(sapply(rdlist, function(.) .[, j]), 1, mean)
    rd
}

#' @exportS3Method emm_basis mira         
emm_basis.mira = function(object, trms, xlev, grid, ...) {
    # In case our method did a "pass it on" with the data, we need to add that attribute
    data = list(...)$misc$data
    if(!is.null(data))
        object$analyses = lapply(object$analyses, function(a) {attr(a, "data") = data; a})
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
    bas$dfargs$m = k
    bas$dfargs$df1 = bas$dffun  # this is the dffun for the 1st analysis
    bas$dfargs$B = cov(t(allb[notna, , drop = FALSE]))
    bas$V = bas$dfargs$T = V + (k + 1)/k * bas$dfargs$B   # pooled via Rubin's rules
    bas$dffun = function(a, dfargs) { # pretty much copied from mice:::barnard.rubin
        dfcom = dfargs$df1(a, dfargs)
        with(dfargs, { 
             b = sum(a * (B %*% a))
             t = sum(a * (T %*% a))
             lambda = (1 + 1/m) * b / t
             dfold = (m - 1)/lambda^2
             dfobs = (dfcom + 1)/(dfcom + 3) * dfcom * (1 - lambda)
             ifelse(is.infinite(dfcom), dfold, 
                    dfold * dfobs/(dfold + dfobs))
        }) }
    bas
}
