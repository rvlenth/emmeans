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

# Support for MCMCglmm class and possibly more MCMC-based models

# Method to create a coda 'mcmc' or 'mcmc.list' object from a ref.grid
# (dots not supported, unfortunately)
# If sep.chains is TRUE and there is more than one chain, an mcmc.list is returned
#' Support for MCMC-based estimation
#' 
#' When a model is fitted using Markov chain Monte Carlo (MCMC) methods, 
#' its reference grid contains a \code{post.beta} slot. These functions 
#' transform those posterior samples to posterior samples of EMMs or
#' related contrasts. They can then be summarized or plotted using,
#' e.g., functions in the \pkg{coda} package.
#'
#' @rdname mcmc-support
#' @aliases mcmc-support
#' @method as.mcmc emmGrid
#' @param x An object of class \code{emmGrid}
#' @param names Logical scalar or vector specifying whether variable names are
#'   appended to levels in the column labels for the \code{as.mcmc} or
#'   \code{as.mcmc.list} result -- e.g., column names of \code{treat A} and
#'   \code{treat B} versus  just \code{A} and \code{B}. When there is more than
#'   one variable involved, the elements of \code{names} are used cyclically.
#' @param sep.chains Logical value. If \code{TRUE}, and there is more than one
#'   MCMC chain available, an \code{\link[coda]{mcmc.list}} object is returned
#'   by \code{as.mcmc}, with separate EMMs posteriors in each chain.
#' @param ... arguments passed to \code{\link[coda]{mcmc}}
#'
#' @return An object of class \code{\link[coda]{mcmc}} or \code{\link[coda]{mcmc.list}}.
#' 
#' @section Details:
#' When the object's \code{post.beta} slot is non-trivial, \code{as.mcmc} will
#' return an \code{\link[coda]{mcmc}} or \code{\link[coda]{mcmc.list}} object
#' that can be summarized or plotted using methods in the \pkg{coda} package.
#' In these functions, \code{post.beta} is transformed by post-multiplying it by
#' \code{t(linfct)}, creating a sample from the posterior distribution of LS
#' means. In \code{as.mcmc}, if \code{sep.chains} is \code{TRUE} and there is in
#' fact more than one chain, an \code{mcmc.list} is returned with each chain's
#' results. The \code{as.mcmc.list} method is guaranteed to return an
#' \code{mcmc.list}, even if it comprises just one chain. 
#' 
#' @importFrom coda as.mcmc
#' @method as.mcmc emmGrid
#' @export
as.mcmc.emmGrid = function(x, names = TRUE, sep.chains = TRUE, ...) {
    object = x
    if (is.na(x@post.beta[1]))
        stop("No posterior sample -- can't make an 'mcmc' object")
    mat = x@post.beta %*% t(x@linfct)
    if(!is.null(offset <- x@grid[[".offset."]])) {
        n = nrow(mat)
        mat = mat + matrix(rep(offset, each = n), nrow = n)
    }
    nm = setdiff(names(x@grid), c(".wgt.",".offset."))
    if (any(names)) {
        names = rep(names, length(nm))
        for (i in seq_along(nm))
            if(names[i]) x@grid[nm[i]] = paste(nm[i], x@grid[[nm[i]]])
    }
    if(is.null(dimnames(mat)))
        dimnames(mat) = list(seq_len(nrow(mat)), seq_len(ncol(mat)))
    dimnames(mat)[[2]] = do.call(paste, c(x@grid[, nm, drop = FALSE], sep=", "))
    n.chains = attr(x@post.beta, "n.chains")
    if (!sep.chains || is.null(n.chains) || (n.chains == 1))
        coda::mcmc(mat, ...)
    else {
        n = nrow(mat) / n.chains
        seqn = seq_len(n)
        chains = lapply(seq_len(n.chains), function(i) coda::mcmc(mat[n*(i - 1) + seqn, , drop = FALSE]))
        coda::mcmc.list(chains, ...)
    }
}

### as.mcmc.list - guaranteed to return a list
#' @rdname mcmc-support
#' @importFrom coda as.mcmc.list
#' @export
#' @method as.mcmc.list emmGrid
as.mcmc.list.emmGrid = function(x, names = TRUE, ...) {
    result = as.mcmc.emmGrid(x, names = names, sep.chains = TRUE, ...)
    if(!inherits(result, "mcmc.list"))
        result = coda::mcmc.list(result)
    result
}



#' Summarize an emmGrid from a Bayesian model
#' 
#' This function computes point estimates and HPD intervals for each
#' factor combination in \code{object@emmGrid}. While this function
#' may be called independently, it is called utomatically by the S3 method
#' \code{\link{summary.emmGrid}} when the object is based on a Bayesian model.
#' (Note: the \code{level} argument, or its default, is passed as \code{prob}).
#'
#' @param object an \code{emmGrid} object having a non-missing \code{post.beta} slot
#' @param prob numeric probability content for HPD intervals (note: when not specified,
#'   the current \code{level} option is used; see \code{\link{emm_options}})
#' @param by factors to use as \code{by} variables
#' @param type prediction type as in \code{\link{summary.emmGrid}}
#' @param point.est function to use to compute the point estimates from the 
#'   posterior sample for each grid point
#' @param ... required but not used
#'
#' @return an object of class \code{summary_emm}
#' @export
#' @importFrom coda HPDinterval
#'
#' @examples
#' load(system.file("extdata", "cbpp.RData", package = "emmeans"))
#' cbpp.rg <- do.call(emmobj, cbpp.list)
#' hpd.summary(emmeans(cbpp.rg, "period"))
#'
hpd.summary = function(object, prob, by, type,
                       point.est = median, ...) {
    if (!requireNamespace("coda"))
        stop("Bayesian summary requires the 'coda' package")
    ### require("coda") ### Nope this is a CRAN no-no
    
    # Steal some init code from summary.emmGrid:
    opt = get_emm_option("summary")
    if(!is.null(opt)) {
        opt$object = object
        object = do.call("update.emmGrid", opt)
    }
    
    misc = object@misc
    use.elts = if (is.null(misc$display))  
                    rep(TRUE, nrow(object@grid)) 
                else                        
                    misc$display
    grid = object@grid[use.elts, , drop = FALSE]
    
    if(missing(prob))
        prob = misc$level
    if(missing(by))
        by = misc$by.vars
    
    if (missing(type))
        type = .get.predict.type(misc)
    else
        type = .validate.type(type)
    # if there are two transformations and we want response, then we need to undo both
    if ((type == "response") && (!is.null(misc$tran2)))
        object = regrid(object, transform = "mu")
    if ((type %in% c("mu", "unlink")) && (!is.null(t2 <- misc$tran2))) {
        if (!is.character(t2))
            t2 = "tran"
        object = update(object, inv.lbl = paste0(t2, "(resp)"))
    }
    
    link = .get.link(misc)
    inv = (type %in% c("response", "mu", "unlink")) # flag to inverse-transform
    if (inv && is.null(link))
        inv = FALSE
    
    
    ### OK, finally, here is the real stuff
    mesg = misc$initMesg
    mcmc = as.mcmc(object, names = FALSE, sep.chains = FALSE)
    mcmc = mcmc[, use.elts, drop = FALSE]
    if (inv) {
        for (j in seq_along(mcmc[1, ]))
            mcmc[, j] = with(link, linkinv(mcmc[, j]))
        mesg = c(mesg, paste("Results are back-transformed from the", link$name, "scale"))
        
    }
    else if(!is.null(link))
        mesg = c(mesg, paste("Results are given on the", link$name, "(not the response) scale."))
    
    mesg = c(mesg, paste("HPD interval probability:", prob))
    pt.est = data.frame(apply(mcmc, 2, point.est))
    names(pt.est) = object@misc$estName
    summ = as.data.frame(coda::HPDinterval(mcmc, prob = prob))[c("lower","upper")]
    names(summ) = cnm = paste0(names(summ), ".HPD")
    lblnms = setdiff(names(grid), 
                     c(object@roles$responses, ".offset.", ".wgt."))
    lbls = grid[lblnms]
    if (inv) {
        if (!is.null(misc$inv.lbl)) {
            names(pt.est) = misc$inv.lbl
            if (!is.null(misc$log.contrast))  # contrast of logs - relabel as ratios
                for (ell in seq_along(lbls)){
                    lbls[[ell]] = factor(lbls[[ell]])
                    levels(lbls[[ell]]) = gsub(" - ", " / ", levels(lbls[[ell]]))
                }
        }
        else
            names(pt.est) = "response"
    }
    
    summ = cbind(lbls, pt.est, summ)
    attr(summ, "estName") = misc$estName
    attr(summ, "clNames") = cnm
    if (is.null(misc$pri.vars) || length(misc$pri.vars) == 0)
        misc$pri.vars = names(object@levels)
    attr(summ, "pri.vars") = setdiff(union(misc$pri.vars, misc$by.vars), by)
    attr(summ, "by.vars") = by
    attr(summ, "mesg") = unique(mesg)
    class(summ) = c("summary_emm", "data.frame")
    summ
}




# Currently, data is required, as call is not stored
recover_data.MCMCglmm = function(object, data, ...) {    
    # if a multivariate response, stack the data with `trait` variable
    yvars = .all.vars(update(object$Fixed$formula, ". ~ 1"))
    if(length(yvars) > 1) {
#        for (v in yvars) data[[v]] = NULL
        dat = data
        for (i in seq_len(length(yvars) - 1))
            data = rbind(data, dat)
        data$trait = factor(rep(yvars, each = nrow(dat)))
    }
    attr(data, "call") = object$Fixed
    attr(data, "terms") = trms = delete.response(terms(object$Fixed$formula))
    attr(data, "predictors") = .all.vars(delete.response(trms))
    attr(data, "responses") = yvars
    data
}

emm_basis.MCMCglmm = function(object, trms, xlev, grid, vcov., ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = NULL)
    Sol = as.matrix(object$Sol)[, seq_len(object$Fixed$nfl)] # toss out random effects if included
    bhat = apply(Sol, 2, mean)
    if (missing(vcov.))
        V = cov(Sol)
    else
        V = .my.vcov(object, vcov.)
    misc = list()
    list(X = X, bhat = bhat, nbasis = matrix(NA), V = V, 
         dffun = function(k, dfargs) Inf, dfargs = list(), 
         misc = misc, post.beta = Sol)
}


### Support for MCMCpack , maybe others that produce mcmc objects
### Whether it works depends on:
###    1. if there is a "call" attribute with a formula or fixed member
###    2. if it's right, even then
### Alternatively, maybe providing formula and data will do the trick

recover_data.mcmc = function(object, formula, data, ...) {
    if (missing(formula)) {
        cl = attr(object, "call")
        if (is.null(cl$formula))
            cl$formula = cl$fixed
        if (is.null(cl$formula))
            return("No fixed-effects formula found")
        data = NULL
    }
    else {
        if (missing(formula) || missing(data))
            return("Requires both formula and data to proceed")
        cl = call("mcmc.proxy", formula = formula, data = quote(data))
    }
    trms = delete.response(terms(eval(cl$formula, parent.frame())))
    recover_data(cl, trms, NULL, data, ...)
}

emm_basis.mcmc = function(object, trms, xlev, grid, vcov., ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = NULL)
    samp = as.matrix(object)[, seq_len(ncol(X)), drop = FALSE]
    bhat = apply(samp, 2, mean)
    if (missing(vcov.))
        V = cov(samp)
    else
        V = .my.vcov(object, vcov.)
    misc = list()
    list(X = X, bhat = bhat, nbasis = matrix(NA), V = V, 
         dffun = function(k, dfargs) Inf, dfargs = list(), 
         misc = misc, post.beta = samp)
}


### Support for mcmc.list
recover_data.mcmc.list = function(object, formula, data, ...) {
    recover_data.mcmc(object[[1]], formula, data, ...)
}

emm_basis.mcmc.list = function(object, trms, xlev, grid, vcov., ...) {
    result = emm_basis.mcmc(object[[1]], trms, xlev, grid, vcov, ...)
    cols = seq_len(ncol(result$post.beta))
    for (i in 2:length(object))
        result$post.beta = rbind(result$post.beta, 
            as.matrix(object[[i]])[, cols, drop = FALSE])
    attr(result$post.beta, "n.chains") = length(object)
    result
}


### support for CARBayes package - currently MUST supply data and have
### default contrasts matching what was used in fitting the mdoel
recover_data.carbayes = function(object, data, ...) {
    if(is.null(data)) # Try to recover data from parent frame
        data = model.frame(object$formula, data = parent.frame())
    cl = call("carbayes.proxy", formula = object$formula, data = quote(data))
    trms = delete.response(terms(eval(object$formula, parent.frame())))
    recover_data(cl, trms, NULL, data, ...)
}

emm_basis.carbayes = function(object, trms, xlev, grid, ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = attr(object$X, "contrasts"))
    samp = as.matrix(object$samples$beta)
    bhat = apply(samp, 2, mean)
    V = cov(samp)
    misc = list()
    list(X = X, bhat = bhat, nbasis = matrix(NA), V = V, 
         dffun = function(k, dfargs) Inf, dfargs = list(), 
         misc = misc, post.beta = samp)
}



### Support for the rstanarm package (stanreg objects)
###
recover_data.stanreg = function(object, ...) {
    recover_data.lm(object, ...)
}

# note: mode and rescale are ignored for some models
emm_basis.stanreg = function(object, trms, xlev, grid, mode, rescale, ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    if(is.null(contr <- object$contrasts))
        contr = attr(model.matrix(object), "contrasts")
    X = model.matrix(trms, m, contrasts.arg = contr)
    bhat = fixef(object)
    V = vcov(object)
    misc = list()
    if (!is.null(object$family)) {
        if (is.character(object$family)) # work around bug for stan_polr
            misc$tran = object$method
        else
            misc = .std.link.labels(object$family, misc)
    }
    if(!is.null(object$zeta)) {   # Polytomous regression model
        if (missing(mode))
            mode = "latent"
        else
            mode = match.arg(mode, 
                             c("latent", "linear.predictor", "cum.prob", "exc.prob", "prob", "mean.class"))
        
        xint = match("(Intercept)", colnames(X), nomatch = 0L)
        if (xint > 0L) 
            X = X[, -xint, drop = FALSE]
        k = length(object$zeta)
        if (mode == "latent") {
            if (missing(rescale)) 
                rescale = c(0,1)
            X = rescale[2] * cbind(X, matrix(- 1/k, nrow = nrow(X), ncol = k))
            bhat = c(bhat, object$zeta - rescale[1] / rescale[2])
            misc = list(offset.mult = rescale[2])
        }
        else {
            bhat = c(bhat, object$zeta)
            j = matrix(1, nrow=k, ncol=1)
            J = matrix(1, nrow=nrow(X), ncol=1)
            X = cbind(kronecker(-j, X), kronecker(diag(1,k), J))
            link = object$method
            if (link == "logistic") link = "logit"
            misc = list(ylevs = list(cut = names(object$zeta)), 
                        tran = link, inv.lbl = "cumprob", offset.mult = -1)
            if (mode != "linear.predictor") {
                misc$mode = mode
                misc$postGridHook = ".clm.postGrid" # we probably need to adapt this
            }
        }
        
        misc$respName = as.character(terms(object))[2]
    }
    samp = as.matrix(object$stanfit)[, names(bhat), drop = FALSE]
    attr(samp, "n.chains") = object$stanfit@sim$chains
    list(X = X, bhat = bhat, nbasis = estimability::all.estble, V = V, 
         dffun = function(k, dfargs) Inf, dfargs = list(), 
         misc = misc, post.beta = samp)
}

