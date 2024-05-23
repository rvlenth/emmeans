##############################################################################
#    Copyright (c) 2012-2021 Russell V. Lenth                                #
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

# experimental support for nls, nlme objects

#' @exportS3Method recover_data nls
recover_data.nls = function(object, ...) {
    fcall = object$call
    trms = terms(.reformulate(names(object$dataClasses)))
    recover_data(fcall, trms, object$na.action, ...)
}

#' @exportS3Method emm_basis nls          
emm_basis.nls = function(object, trms, xlev, grid, ...) {
    Vbeta = .my.vcov(object, ...)
    env = object$m$getEnv()
    for (nm in names(grid)) env[[nm]] = grid[[nm]]
    pars = object$m$getAllPars()
    DD = deriv(object$m$formula(), names(pars))
    ests = eval(DD, env)
    bhat = as.numeric(ests)
    grad = attr(ests, "gradient")
    V = grad %*% Vbeta %*% t(grad)
    X = diag(1, nrow(grid))
    list(X=X, bhat=bhat, nbasis=all.estble, V=V, 
         dffun=function(k, dfargs) Inf, dfargs=list(), 
         misc=list())
}
    

### For nlme objects, we can do stuff with the fixed part of the model
### Additional REQUIRED argument is 'param' - parameter name to explore
#' @exportS3Method recover_data nlme
recover_data.nlme = function(object, param, ...) {
    if(missing(param))
        return("'param' argument is required for nlme objects")
    fcall = object$call
    if (!is.null(fcall$weights))
        fcall$weights = nlme::varWeights(object$modelStruct)
    fixed = fcall$fixed
    if (is.call(fixed))
        fixed = eval(fixed, envir = parent.frame())
    if(!is.list(fixed))
        fixed = list(fixed)
    form = NULL
    for (x in fixed)
        if (param %in% all.names(x)) form = x
    if (is.null(form))
        return(paste("Can't find '", param, "' among the fixed parameters", sep = ""))
    fcall$weights = NULL
    trms = delete.response(terms(form))
    if (length(.all.vars(trms)) == 0)
        return(paste("No predictors for '", param, "' in fixed model", sep = ""))
    recover_data(fcall, trms, object$na.action, ...)
}

#' @exportS3Method emm_basis nlme         
emm_basis.nlme = function(object, trms, xlev, grid, param, ...) {
    idx = object$map$fmap[[param]]
    V = object$varFix[idx, idx, drop = FALSE]
    bhat = object$coefficients$fixed[idx]
    contr = attr(object$plist[[param]]$fixed, "contrasts")
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contr)
    dfx = object$fixDF$X[idx]
    dfx[1] = min(dfx) # I'm assuming 1st one is intercept
    dffun = function(k, dfargs) { # containment df
        idx = which(abs(k) > 1e-6)
        ifelse(length(idx) > 0, min(dfargs$dfx[idx]), NA)
    }
    list(X = X, bhat = bhat, nbasis = estimability::all.estble, 
         V = V, dffun = dffun, dfargs = list(dfx = dfx), 
         misc = list(estName = param))
}


# Support for gnls - contributed by Fernando Miguez <femiguez@iastate.edu>
#' @exportS3Method recover_data gnls
recover_data.gnls = function(object, param, data, ...) {
    fcall = object$call$params
    if (is.null(fcall))
        return("Models fitted without a 'params' specification are not supported.")
    if(missing(param))
        return("'param' argument is required for gnls objects")

    plist = object$plist
    pnames = names(plist)
    params = eval(object$call$params)
    if (!is.list(params)) params = list(params)
    
    params = unlist(lapply(params, function(pp) {
        if (is.name(pp[[2]])){
            list(pp)
        }
        else { ## multiple parameters on left hand side
            eval(parse(text = paste("list(",
                                    paste(paste(all.vars(pp[[2]]), deparse(pp[[3]]), sep = "~"),
                                          collapse = ","),
                                    ")")))
        }
    }), recursive = FALSE)

    names(params) = pnames
    form = params[[param]]
    
    trms = delete.response(terms(eval(form, envir = environment(formula(object)))))
    if(is.null(data))
        data = eval(object$call$data, envir = environment(formula(object)))
    recover_data(fcall, trms, object$na.action, data = data, ...)
}

#' @exportS3Method emm_basis gnls         
emm_basis.gnls = function(object, trms, xlev, grid, param, ...) {
    idx = object$pmap[[param]]
    V = object$varBeta[idx, idx, drop = FALSE]
    bhat = object$coefficients[idx]
    contr = attr(object$plist[[param]], "contrasts")
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = contr)
    dfx = with(attributes(logLik(object)), nobs - df)
    dffun = function(k, dfargs) dfargs$dfx
    list(X = X, bhat = bhat, nbasis = estimability::all.estble, 
         V = V, dffun = dffun, dfargs = list(dfx = dfx), 
         misc = list(estName = param))
}

