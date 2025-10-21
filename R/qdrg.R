##############################################################################
#    Copyright (c) 2018-2024 Russell V. Lenth                                     #
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

### Quick-and-dirty support for otherwise unsupported models

#' Quick and dirty reference grid
#' 
#' This function may make it possible to compute a reference grid for a model 
#' object that is otherwise not supported.
#' 
#' Usually, you need to provide either \code{object}; or
#' \code{formula}, \code{coef}, \code{vcov}, \code{data}, and perhaps other
#' parameters. It is often fairly straightforward to figure out how to get
#' these from the model \code{object}; see the documentation for the model class that
#' was fitted. Sometimes one or more of these quantities contains extra parameters,
#' and if so, you may need to subset them to make everything conformable. For a given \code{formula} and \code{data},
#' you can find out what is needed via \code{colnames(model.matrix(formula, data))}.
#' (However, for an ordinal model, we expect the first \code{ordinal.dim - 1} coefficients
#' to replace \code{(Intercept)}. And for a multivariate model, we expect \code{coef} 
#' to be a matrix with these row names, and \code{vcov} to have as many rows and columns as
#' the total number of elements of \code{coef}.)
#' 
#' If \code{object} is specified, this function serves as a generic for dispatching
#' a method for \code{object}'s class. See the arguments for \code{qdrg.default} to
#' see how the arguments are determined by default. For many \code{lm}- or \code{glm}-like
#' models, it may suffice override just one or two of these
#' defaults. The possibility of a custom \code{qdrg} method also provides a
#' minimal way for package developers to provide \pkg{emmeans} support: it doesn't
#' allow directly applying \code{emmeans()} on the model, but at least the user
#' can obtain a reference grid, and then go from there. Note that when using a `qdrg` 
#' method, it is best for the user to specify `object =` explicitly in the
#' call, since `object` is not the first argument of `qdrg()`. However, if the
#' first argument is not a formula, the function is retried with it as \code{object}.
#' 
#' The functions \code{\link{qdrg}} and \code{emmobj} are close cousins, in that
#' they both produce \code{emmGrid} objects. When starting with summary
#' statistics for an existing grid, \code{emmobj} is more useful, while
#' \code{qdrg} is more useful when starting from a fitted model.
#'
#' @param formula Formula for the fixed effects
#' @param data Dataset containing the variables in the model
#' @param coef Fixed-effect regression coefficients (must conform to formula)
#' @param vcov Variance-covariance matrix of the fixed effects
#' @param df Error degrees of freedom
#' @param mcmc Posterior sample of fixed-effect coefficients
#' @param object Optional model object. \emph{This rarely works!}; 
#'        but if provided, we try to set 
#'        other arguments based on an expectation that `object` has a similar
#'        structure to `lm` objects. See Details.
#' @param subset Subset of \code{data} used in fitting the model
#' @param weights Weights used in fitting the model
#' @param contrasts List of contrasts specified in fitting the model
#' @param link Link function (character or list) used, if a generalized linear model.
#'     (Note: response transformations are auto-detected from \code{formula})
#' @param qr QR decomposition of the model matrix; used only if there are \code{NA}s
#'     in \code{coef}.
#' @param ordinal \code{list} with elements \code{dim} and \code{mode}.
#'     \code{ordinal$dim} (integer) is the number of levels in an ordinal response. If 
#'     \code{ordinal} is provided, the intercept terms are modified appropriate to predicting 
#'     an ordinal response, as described in \code{vignette("models")}, Section O,
#'     using \code{ordinal$mode} as the \code{mode} argument (if not
#'     provided, \code{"latent"} is assumed).
#'     (All modes are supported except `scale`)
#'     For this to work, we expect
#'     the first \code{ordinal$dim - 1} elements of \code{coef} to be the
#'     estimated threshold parameters, followed by the coefficients for the
#'     linear predictor. Also, if \code{mode} requires back-transforming (e.g.,
#'     \code{"prob"} or \code{"mean.class"}), the user may need to supply the \code{link}
#'     for it to work correctly.
#' @param ... Optional arguments passed to \code{\link{ref_grid}}
#'
#' @return An \code{emmGrid} object constructed from the arguments
#' 
#' @section Rank deficiencies:
#' Different model-fitting packages take different approaches when the model
#' matrix is singular, but \code{qdrg} tries to reconcile them by comparing the
#' linear functions created by \code{formula} to \code{coefs} and \code{vcov}.
#' We may then use the \pkg{estimability} package to determine what quantities
#' are estimable. For reconciling to work properly, \code{coef} should be named
#' and \code{vcov} should have dimnames. To disable this name-matching
#' action, remove the names from \code{coef}, e.g., by calling \code{unname()}.
#' No reconciliation is attempted in multivariate-response cases. For more
#' details on estimability, see the documentation in the \pkg{estimability}
#' package.
#' 
#' @seealso \code{\link{emmobj}} for an alternative way to construct an \code{emmGrid}.
#' 
#' @note For backwards compatibility, an argument \code{ordinal.dim} is invisibly 
#' supported as part of \code{...}, and if present, sets 
#' \code{ordinal = list(dim = ordinal.dim, mode = "latent")}
#' 
#' @export
#' @examples
#' # In these examples, use emm_example(..., list = TRUE) # to see just the code
#' 
#' if (require(biglm, quietly = TRUE)) 
#'     emm_example("qdrg-biglm")
#'     
#' if(require(coda, quietly = TRUE) && require(lme4, quietly = TRUE)) 
#'     emm_example("qdrg-coda")
#'     
#' if(require(ordinal, quietly = TRUE)) 
#'     emm_example("qdrg-ordinal")
#'
qdrg = function(formula, data, coef, vcov, df, mcmc, object,
                subset, weights, contrasts, link, qr, ordinal, ...) {
    
    if (!missing(object)) {
        UseMethod("qdrg", object = object)
    }
    else { 
        if(!inherits(formula, "formula")) {
            cl = match.call()
            cl$object = formula
            cl$formula = NULL
            ### message("Retrying with first argument taken as 'object'...")
            return(eval(cl))
        }
        # back-compatible access to old ordinal.dim arg...
        od = (\(ordinal, ordinal.dim = NULL, ...) {
            if(!missing(ordinal) && is.numeric(ordinal)) ordinal.dim = ordinal
            ordinal.dim
        })(ordinal, ...)
        if(!is.null(od)) ordinal = list(dim = od, mode = "latent")
        
        result = match.call()
        if(missing(formula))
            stop("When 'object' is missing, must at least provide 'formula'")
        result$formula = formula
        if(missing(data))
            result$data = parent.frame()
        else
            result$data = data
        if (!missing(coef)) result$coef = coef
        if (!missing(vcov)) result$vcov = vcov
        if(!missing(df)) result$df = df
        if(missing(contrasts))
            contrasts = attr(model.matrix(result$formula, data = data), "contrasts")
        
        if(!missing(df)) result$df = df
        if(is.null(result$df))
            result$df = Inf
        if(!missing(mcmc)) result$mcmc = mcmc
        if(!missing(subset)) result$subset = subset
        if(!missing(weights)) result$weights = weights
        if(!missing(contrasts)) result$contrasts = contrasts
        if(!missing(link)) result$link = link
        if(!missing(qr) && any(is.na(result$coef))) result$qr = qr
        if(!missing(ordinal)) result$ordinal = ordinal
        
        # make sure "formula" exists, has a LHS and is is 2nd element so that 
        # response transformation can be found
        if (is.null(result$formula))
            stop("No formula; cannot construct a reference grid")
        if(length(result$formula) < 3)
            result$formula = update.formula(result$formula, response ~ .)
        fpos = grep("formula", names(result))[1]
        result = result[c(1, fpos, seq_along(result)[-c(1, fpos)])]
        
        class(result) = c("qdrg", "call")
        ref_grid(result, ...)
    }
}

#' @rdname qdrg
#' @exportS3Method qdrg default
qdrg.default = function(formula = stats::formula(object), 
                        data = try(recover_data.lm(object), silent = TRUE), 
                        coef = stats::coef(object),
                        vcov = stats::vcov(object), 
                        df = stats::df.residual(object), 
                        mcmc, 
                        object,
                        subset, 
                        weights = stats::weights(object), 
                        contrasts = object$contrasts, 
                        link = ifelse(!is.null(lnk<-object$family$link), lnk, object$link),
                        qr = object$qr, 
                        ordinal, 
                        ...) 
{
    
    if(inherits(data, "try-error")) {
        if(is.null(data <- object$data))
            stop("Unable to recover data. You must specify it explicitly in the 'data' argument.")
    }
    if(missing(mcmc)) mcmc = NULL   # for some weird reason, this is needed
    if(missing(subset)) subset = NULL
    if(missing(ordinal)) ordinal = NULL
    qdrg(formula = formula, data = data, coef = coef, vcov = vcov, df = df, mcmc = mcmc,
         subset = subset, weights = weights, contrasts = contrasts, link = link,
         qr = qr, ordinal = ordinal, ...)
}



#' @exportS3Method recover_data qdrg
recover_data.qdrg = function(object, ...) {
    recover_data.call(object, delete.response(terms(object$formula)), object$na.action, ...)
}

#' @exportS3Method vcov qdrg
vcov.qdrg = function(object, ...) 
    object$vcov

#' @exportS3Method emm_basis qdrg         
emm_basis.qdrg = function(object, trms, xlev, grid, ...) {
    m = suppressWarnings(model.frame(trms, grid, na.action = na.pass, xlev = xlev))
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    if (is.null(object$mcmc)) {
        bhat = object$coef
        V = .my.vcov(object, ...)
    }
    else {
        if (is.null(object$coef)) bhat = apply(object$mcmc, 2, mean)
        if (is.null(object$vcov)) V = cov(object$mcmc)
    }
    
    misc = list()
    
    # If ordinal, add extra avgd, subtracted intercepts -- for latent mode
    if(!is.null(ordinal <- object$ordinal)) {
        if(is.null(ordinal$mode)) ordinal$mode = "latent"
        ordinal$mode = match.arg(ordinal$mode, c("latent", "linear.predictor", "cum.prob", "exc.prob", "prob", "mean.class"))
        if(is.null(od <- ordinal$dim)) stop ("'ordinal' MUST have a 'dim' element", call. = FALSE)
        if(ordinal$mode == "latent") {
            intcpt = matrix(-1 / (od - 1), nrow = nrow(X), ncol = od - 1)
            colnames(intcpt) = names(bhat)[1:(od - 1)]
            X = cbind(intcpt, X[, -1, drop = FALSE])
        }
        else {
            levs = seq_len(od - 1)
            misc$ylevs = list(cut = paste(levs, levs+1, sep = "|"))
            misc$inv.lbl = "cumprob"
            misc$offset.mult = -1
            J = matrix(1, nrow = od - 1)
            X = cbind(kronecker(diag(od - 1), matrix(1, nrow = nrow(X))), 
                      kronecker(-J, X[, -1, drop = FALSE]))
            if (ordinal$mode != "linear.predictor") {
                misc$mode = ordinal$mode
                misc$respName = as.character(object$formula)[2]
                misc$postGridHook = ".clm.postGrid"
            }
        }
    }
    
    bhat = .impute.NAs(bhat, X) # make coefs lm-compatible
    nbasis = estimability::all.estble
    if (sum(is.na(bhat)) > 0) {
        if(!is.null(object$qr))
            nbasis = estimability::nonest.basis(object$qr)
        else {
            if (is.name(object$data))
                object$data = eval(object$data)
            mm = suppressWarnings(model.frame(trms, object$data, na.action = na.pass, xlev = xlev))
            XX = model.matrix(trms, mm, contrasts.arg = object$contrasts)
            nbasis = estimability::nonest.basis(XX)
        }
        if (nrow(V) == length(bhat)) {
            ii = which(!is.na(bhat))
            V = V[ii, ii, drop = FALSE]
        }
    }
    
    # check multivariate situation
    if (is.matrix(bhat)) {
        X = kronecker (diag(ncol(bhat)), X)
        nbasis = kronecker(rep(1, ncol(bhat)), nbasis)
        nms = colnames(bhat)
        if (is.null(nms))
            nms = seq_len(ncol(bhat))
        misc$ylevs = list(rep.meas = nms)
        bhat = as.numeric(bhat)
    }
    
    if (!is.null(object$link))
        misc = .std.link.labels(eval(list(link = object$link)), misc)
    dfargs = list(df = object$df)
    dffun = function(k, dfargs) dfargs$df
    
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, 
         misc=misc, post.beta=object$mcmc)
    }
    