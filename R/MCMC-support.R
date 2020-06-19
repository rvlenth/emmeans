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

# NOTE: S3 registration of as.mcmc and as.mcmc.list is done dynamically in zzz.R

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
#' @param x An object of class \code{emmGrid}
#' @param names Logical scalar or vector specifying whether variable names are
#'   appended to levels in the column labels for the \code{as.mcmc} or
#'   \code{as.mcmc.list} result -- e.g., column names of \code{treat A} and
#'   \code{treat B} versus  just \code{A} and \code{B}. When there is more than
#'   one variable involved, the elements of \code{names} are used cyclically.
#' @param sep.chains Logical value. If \code{TRUE}, and there is more than one
#'   MCMC chain available, an \code{\link[coda]{mcmc.list}} object is returned
#'   by \code{as.mcmc}, with separate EMMs posteriors in each chain.
#' @param likelihood Character value or function. If given, simulations are made from 
#'   the corresponding posterior predictive distribution. If not given, we obtain
#'   the posterior distribution of the parameters in \code{object}. See Prediction
#'   section below.
#' @param NE.include Logical value. If \code{TRUE}, non-estimable columns are
#'   kept but returned as columns of \code{NA} values (this may create errors or
#'   warnings in subsequent analyses using, say, \pkg{coda}). If \code{FALSE},
#'   non-estimable columns are dropped, and a warning is issued. (If all are
#'   non-estimable, an error is thrown.)
#' @param ... arguments passed to other methods
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
#' @section Prediction:
#' When \code{likelihood} is specified, it is used to simulate values from the
#' posterior predictive distribution corresponding to the given likelihood and
#' the posterior distribution of parameter values. Denote the likelihood 
#' function as \eqn{f(y|\theta,\phi)}, where \eqn{y} is a response, \eqn{\theta}
#' is the parameter estimated in \code{object}, and \eqn{\phi} comprises zero or
#' more additional parameters to be specified. If \code{likelihood} is a 
#' function, that function should take as its first argument a vector of 
#' \eqn{\theta} values (each corresponding to one row of \code{object@grid}).
#' Any \eqn{\phi} values should be specified as additional named function
#' arguments, and passed to \code{likelihood} via \code{...}. This function should 
#' simulate values of \eqn{y}.
#' 
#' A few standard likelihoods are available by specifying \code{likelihood} as
#' a character value. They are:
#' \describe{
#'   \item{\code{"normal"}}{The normal distribution with mean \eqn{\theta} and
#'   standard deviation specified by additional argument \code{sigma}}
#'   \item{\code{"binomial"}}{The binomial distribution with success probability 
#'     \eqn{theta}, and number of trials specified by \code{trials}}
#'   \item{\code{"poisson"}}{The Poisson distribution with mean \eqn{theta} 
#'     (no additional parameters)}
#'   \item{\code{"gamma"}}{The gamma distribution with scale parameter \eqn{\theta}
#'     and shape parameter specified by \code{shape}}
#' }
#' 
#' @method as.mcmc emmGrid
#' @export as.mcmc.emmGrid
#' @examples
#' require("coda")
#'
#' ### A saved reference grid for a mixed logistic model (see lme4::cbpp)
#' cbpp.rg <- do.call(emmobj, 
#'     readRDS(system.file("extdata", "cbpplist", package = "emmeans")))
#' # Predictive distribution for herds of size 20
#' # (perhaps a bias adjustment should be applied; see "sophisticated" vignette)
#' pred.incidence <- as.mcmc(regrid(cbpp.rg), likelihood = "binomial", trials = 20)
as.mcmc.emmGrid = function(x, names = TRUE, sep.chains = TRUE, 
                           likelihood, NE.include = FALSE, ...) {
    if (is.na(x@post.beta[1])) {
        stop("No posterior sample -- can't make an 'mcmc' object")
    }
# notes on estimabilityn issues:
# 1. Use @bhat to determine which coefs to use
# 2. @nabasis as in freq models
# 3. @post.beta we will EXCLUDE cols corresp to NAs in @bhat
# See stanreg support for hints/details
    use = which(!is.na(x@bhat))
    est = estimability::is.estble(x@linfct, x@nbasis)
    if (!any(est))
        stop("Aborted -- No estimates in the grid are estimable")
    else if(!all(est) && !NE.include) {
        rows = paste(which(!est), collapse = ", ")
        warning("Cases  ", rows, "  were dropped due to non-estimability", call. = FALSE)
    }
    mat = x@post.beta %*% t(x@linfct[, use, drop = FALSE])
    if (NE.include)
        mat[, !est] = NA
    else {
        mat = mat[, est, drop = FALSE]
        x@grid = x@grid[est, , drop = FALSE]
    }
    if(!is.null(offset <- x@grid[[".offset."]])) {
        n = nrow(mat)
        mat = mat + matrix(rep(offset, each = n), nrow = n)
    }
    if (!missing(likelihood)) {
        if (is.character(likelihood)) {
            likelihood = match.arg(likelihood, c("normal", "binomial", "poisson", "gamma"))
            likelihood = switch(likelihood,
                normal = function(theta, sigma, ...) 
                    rnorm(length(theta), mean = theta, sd = sigma),
                binomial = function(theta, trials, ...) 
                    rbinom(length(theta), size = trials, prob = theta),
                poisson = function(theta, ...)
                    rpois(length(theta), lambda = theta),
                gamma = function(theta, shape, ...)
                    rgamma(length(theta), scale = theta, shape = shape)
                #, stop("There is no predefined likelihood named '", likelihood, "'")
                )
        }
        mat = apply(mat, 2, likelihood, ...)
##! TODO: Add "multinomial" support. This will require a flag to observe
##! the 'by' variable(s), then we get parameter values from the columns
##! corresponding to each 'by' group 
    }
    nm = setdiff(names(x@grid), c(".wgt.",".offset."))
    if (any(names)) {
        names = rep(names, length(nm))
        for (i in seq_along(nm))
            if(names[i]) x@grid[nm[i]] = paste(nm[i], x@grid[[nm[i]]])
    }
    if(is.null(dimnames(mat)))
        dimnames(mat) = list(seq_len(nrow(mat)), seq_len(ncol(mat)))
    dimnames(mat)[[2]] = do.call(paste, c(unname(x@grid[, nm, drop = FALSE]), sep=", "))
    n.chains = attr(x@post.beta, "n.chains")
    if (!sep.chains || is.null(n.chains) || (n.chains == 1))
        coda::mcmc(mat)
    else {
        n = nrow(mat) / n.chains
        seqn = seq_len(n)
        chains = lapply(seq_len(n.chains), function(i) coda::mcmc(mat[n*(i - 1) + seqn, , drop = FALSE]))
        coda::mcmc.list(chains)
    }
}


### as.mcmc.list - guaranteed to return a list
#' @rdname mcmc-support
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
#' may be called independently, it is called automatically by the S3 method
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
#' @param bias.adjust Logical value for whether to adjust for bias in
#'   back-transforming (\code{type = "response"}). This requires a value of 
#'   \code{sigma} to exist in the object or be specified.
#' @param sigma Error SD assumed for bias correction (when 
#'   \code{type = "response"}. If not specified,
#'   \code{object@misc$sigma} is used, and an error is thrown if it is not found.
#'   \emph{Note:} \code{sigma} may be a vector, as long as it conforms to the 
#'   number of observations in the posterior sample.
#' @param ... required but not used
#'
#' @return an object of class \code{summary_emm}
#' 
#' @seealso summary.emmGrid
#' 
#' @export
#'
#' @examples
#' if(require("coda")) {
#'   # Create an emmGrid object from a system file
#'   cbpp.rg <- do.call(emmobj, 
#'       readRDS(system.file("extdata", "cbpplist", package = "emmeans")))
#'   hpd.summary(emmeans(cbpp.rg, "period"))
#' }
#' 
hpd.summary = function(object, prob, by, type, point.est = median, 
                       bias.adjust = get_emm_option("back.bias.adj"), sigma, 
                       ...) {
    if(!is.null(object@misc$.predFlag))
        stop("Prediction intervals for MCMC models should be done using 'frequentist = TRUE'\n",
             "or using 'as.mcmc(object, ..., likelihood = ...)'")
    
    .requireNS("coda", "Bayesian summary requires the 'coda' package")
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
    pe.lbl = as.character(substitute(point.est))
    if(length(pe.lbl) > 1) 
        pe.lbl = "user-supplied function"
    mesg = c(misc$initMesg, paste("Point estimate displayed:", pe.lbl))
    mcmc = as.mcmc.emmGrid(object, names = FALSE, sep.chains = FALSE, 
                           NE.include = TRUE, ...)
    mcmc = mcmc[, use.elts, drop = FALSE]
    if (inv) {
        if (bias.adjust) {
            if (missing(sigma))
                sigma = object@misc@sigma
            link = .make.bias.adj.link(link, sigma)
        }
        
        for (j in seq_along(mcmc[1, ]))
            mcmc[, j] = with(link, linkinv(mcmc[, j]))
        mesg = c(mesg, paste("Results are back-transformed from the", link$name, "scale"))
        if(bias.adjust)
            mesg = c(mesg, paste("Bias adjustment applied based on sigma =",
                                 .fmt.sigma(sigma)))
    }
    else if(!is.null(link))
        mesg = c(mesg, paste("Results are given on the", link$name, "(not the response) scale."))
    
    est = !is.na(mcmc[1, ])
    mcmc[, !est] = 0 # temp so we don't get errors
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
            names(pt.est) = misc$estName = misc$inv.lbl
            if (!is.null(misc$log.contrast))  # contrast of logs - relabel as ratios
                for (ell in seq_along(lbls)){
                    lbls[[ell]] = factor(lbls[[ell]])
                    levels(lbls[[ell]]) = gsub(" - ", " / ", levels(lbls[[ell]]))
                }
        }
        else
            names(pt.est) = misc$estName = "response"
    }
    
    summ[!est, ] = NA
    pt.est[!est, ] = NA
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
recover_data.MCMCglmm = function(object, data, trait, ...) {
    if (is.null(data) && !is.null(object$data)) # allow for including data in object
        data = eval(object$data)
    # if a multivariate response, stack the data with `trait` variable
    yvars = .all.vars(update(object$Fixed$formula, ". ~ 1"))
    if ("trait" %in% names(data)) {
        # don't do anything, just use what's provided
    }
    else if(length(yvars) > 1) {
#        for (v in yvars) data[[v]] = NULL
        dat = data
        for (i in seq_len(length(yvars) - 1))
            data = rbind(data, dat)
        data$trait = factor(rep(yvars, each = nrow(dat)))
    }
    else if(!missing(trait)) {
        # we'll create a fake "trait" variable with specified variable
        n = nrow(data)
        levs = levels(data[[trait]])
        attr(data, "misc") = list(resp.levs = levs, trait = trait)
        data$trait = rep(levs[-1], n)[1:n] # way overkill, but easy coding
    }
    attr(data, "call") = object$Fixed
    attr(data, "terms") = trms = delete.response(terms(object$Fixed$formula))
    attr(data, "predictors") = .all.vars(delete.response(trms))
    data
}

# misc may be NULL or a list generated by trait spec
emm_basis.MCMCglmm = function(object, trms, xlev, grid, vcov., 
                              mode = c("default", "multinomial"), misc, ...) {
    nobs.MCMCglmm = function(object, ...) 1   # prevents warning about nobs
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = NULL)
    Sol = as.matrix(object$Sol)[, seq_len(object$Fixed$nfl)] # toss out random effects if included
    bhat = apply(Sol, 2, mean)
    if (missing(vcov.))
        V = cov(Sol)
    else
        V = .my.vcov(object, vcov.)
    if (is.null(misc))
        misc = list()
    mode = match.arg(mode)
    if (mode == "multinomial") {
        misc$postGridHook = .MCMCglmm.multinom.postGrid
    }
    else { # try to figure out the link
        fam = unique(object$family)
        if (length(fam) > 1)
            stop("There is more than one 'family' in this model - too complex for emmeans support")
        link = switch(fam,
                      poisson = "log",
                      multinomial = "log",
                      categorical = "logit",
                      ordinal = "logit") # maybe more later?
        if (!is.null(link))
            misc = .std.link.labels(list(link = link), misc)
    }
    list(X = X, bhat = bhat, nbasis = matrix(NA), V = V, 
         dffun = function(k, dfargs) Inf, dfargs = list(), 
         misc = misc, post.beta = Sol)
}



.MCMCglmm.multinom.postGrid = function(object, ...) {
    linfct = object@linfct
    misc = object@misc
    post.lp = object@post.beta %*% t(linfct)
    sel = .find.by.rows(object@grid, "trait")
    k = length(sel)
    cols = unlist(sel)
    scal = sqrt(1 + 2 * (16 * sqrt(3) / (15 * pi))^2 / (k + 1))  # scaling const for logistic
    # I'm assuming here that diag(IJ) = 2 / (k + 1)
    object@post.beta = post.p = t(apply(post.lp, 1, function(l) {
        expX = exp(cbind(0, matrix(l[cols], ncol = k)) / scal)
        as.numeric(apply(expX, 1, function(z) z / sum(z)))
    })) # These results come out with response levels varying the fastest.
    
    object@bhat = apply(post.p, 2, mean)
    object@V = cov(post.p)
    preds = c(misc$trait, object@roles$predictors)
    object@roles$predictors  = preds[preds != "trait"]
    object@levels[["trait"]] = NULL
    object@levels = c(list(misc$resp.levs), object@levels)
    names(object@levels)[1] = misc$trait
    object@grid = do.call(expand.grid, object@levels)
    
    misc$postGridHook = misc$tran = misc$inv.lbl = 
        misc$trait = misc$resp.levs = NULL
    misc$display = object@model.info$nesting = NULL
    misc$estName = "prob"
    object@linfct = diag(1, ncol(post.p))
    object@misc = misc
    
    object
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
    misc = list()
    if (!is.null(object$family)) {
        if (is.character(object$family)) # work around bug for stan_polr
            misc$tran = object$method
        else
            misc = .std.link.labels(object$family, misc)
    }
    # Previous code...
    ### m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    ### if(is.null(contr <- object$contrasts))
    ###     contr = attr(model.matrix(object), "contrasts")
    ### X = model.matrix(trms, m, contrasts.arg = contr)
    ### bhat = rstanarm::fixef(object)
    ### nms = intersect(colnames(X), names(bhat))
    ### bhat = bhat[nms]
    ### V = vcov(object)[nms, nms, drop = FALSE]

    # Instead, use internal routine in rstanarm to get the model matrix
    # Later, we'll get bhat and V from the posterior sample because
    # the vcov(object) doesn't always jibe with fixef(object)
    pp_data = get("pp_data", envir = getNamespace("rstanarm"))
    X = pp_data(object, newdata = grid, re.form = ~0, ...)[[1]]
    nms = colnames(X)
    
    if(!is.null(object$zeta)) {   # Polytomous regression model
        if (missing(mode))
            mode = "latent"
        else
            mode = match.arg(mode, 
                             c("latent", "linear.predictor", "cum.prob", "exc.prob", "prob", "mean.class"))
        
        xint = match("(Intercept)", nms, nomatch = 0L)
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
        
        misc$respName = as.character.default(terms(object))[2]
    }
    samp = as.matrix(object$stanfit)[, nms, drop = FALSE]
    attr(samp, "n.chains") = object$stanfit@sim$chains

    bhat = apply(samp, 2, mean)
    V = cov(samp)
    
    # estimability...
    nbasis = estimability::all.estble
    all.nms = colnames(X)
    if (length(nms) < length(all.nms)) {
        if(is.null(contr <- object$contrasts))
            contr = attr(model.matrix(object), "contrasts")
        coef = NA * X[1, ]
        coef[names(bhat)] = bhat
        bhat = coef
        mmat = model.matrix(trms, object$data, contrasts.arg = contr)
        nbasis = estimability::nonest.basis(mmat)
    }
    
    
    list(X = X, bhat = bhat, nbasis = nbasis, V = V, 
         dffun = function(k, dfargs) Inf, dfargs = list(), 
         misc = misc, post.beta = samp)
}

