##############################################################################
#    Copyright (c) 2012-2017 Russell V. Lenth                                #
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

### Code for emtrends


### emtrends function
#' Estimated marginal means of linear trends
#'
#' The \code{emtrends} function is useful when a fitted model involves a
#' numerical predictor \eqn{x}  interacting with another predictor \code{a}
#' (typically a factor). Such models specify that \eqn{x} has a different trend
#' depending on \eqn{a}; thus, it may be of interest to estimate and compare
#' those trends. Analogous to the \code{\link{emmeans}} setting, we construct a
#' reference grid of these predicted trends, and then possibly average them over
#' some of the predictors in the grid.
#' 
#' @param model A supported model object (\emph{not} a reference grid)
#' @param specs Specifications for what marginal trends are desired -- as in
#'   \code{\link{emmeans}}
#' @param var Character value giving the name of a variable with respect to 
#'   which a difference quotient of the linear predictors is computed. In order
#'   for this to be useful, \code{var} should be a numeric predictor that
#'   interacts with at least one factor in \code{specs}. Then instead of
#'   computing EMMs, we compute and compare the slopes of the \code{var} trend
#'   over levels of the specified other predictor(s). As in EMMs, marginal
#'   averages are computed for the predictors in \code{specs} and \code{by}.
#'   See also the \dQuote{Generalizations} section below.
#' @param delta.var The value of \emph{h} to use in forming the difference
#'   quotient \eqn{(f(x+h) - f(x))/h}. Changing it (especially changing its
#'   sign) may be necessary to avoid numerical problems such as logs of negative
#'   numbers. The default value is 1/1000 of the range of \code{var} over the
#'   dataset.
#' @param data As in \code{\link{ref_grid}}, you may use this argument to supply
#'   the dataset used in fitting the model, for situations where it is not
#'   possible to reconstruct the data. Otherwise, leave it missing.
#' @param transform If \code{object} has a response
#'   transformation or link function, then specifying 
#'   \code{transform = "response"} will cause
#'   \code{emtrends} to calculate the trends after back-transforming to the
#'   response scale. This is done using the chain rule, and standard errors are
#'   estimated via the delta method. With \code{transform = "none"} (the
#'   default), the trends are calculated on the scale of the linear predictor,
#'   without back-transforming it. This argument works similarly to the
#'   \code{transform} argument in \code{\link{ref_grid}}, in that the returned
#'   object is re-gridded to the new scale (see also \code{\link{regrid}}).
#' @param max.degree Integer value. The maximum degree of trends to compute (this
#'   is capped at 5). If greater than 1, an additional factor \code{degree} is
#'   added to the grid, with corresponding numerical derivatives of orders
#'   \code{1, 2, ..., max.degree} as the estimates.
#' @param ... Additional arguments passed to other methods or to 
#'   \code{\link{ref_grid}}
#'
#' @section Generalizations:
#' Instead of a single predictor, the user may specify some monotone function of
#' one variable, e.g., \code{var = "log(dose)"}. If so, the chain rule is
#' applied. Note that, in this example, if \code{model} contains
#' \code{log(dose)} as a predictor, we will be comparing the slopes estimated by
#' that model, whereas specifying \code{var = "dose"} would perform a
#' transformation of those slopes, making the predicted trends vary depending on
#' \code{dose}.
#' 
#' @return An \code{emmGrid} or \code{emm_list} object, according to \code{specs}.
#' See \code{\link{emmeans}} for more details on when a list is returned.
#' 
#' @seealso \code{link{emmeans}}, \code{\link{ref_grid}}
#' @export
#'
#' @examples
#' fiber.lm <- lm(strength ~ diameter*machine, data=fiber)
#' # Obtain slopes for each machine ...
#' ( fiber.emt <- emtrends(fiber.lm, "machine", var = "diameter") )
#' # ... and pairwise comparisons thereof
#' pairs(fiber.emt)
#' 
#' # Suppose we want trends relative to sqrt(diameter)...
#' emtrends(fiber.lm, ~ machine | diameter, var = "sqrt(diameter)", 
#'          at = list(diameter = c(20, 30)))
#'
emtrends = function(model, specs, var, delta.var=.001*rng, data, 
                    transform = c("none", "response"), max.degree = 1, ...) {
    estName = paste(var, "trend", sep=".") # Do now as I may replace var later
    
    if (missing(data)) {
        data = try(recover_data (model, data = NULL))
        if (inherits(data, "try-error"))
            stop("Possible remedy: Supply the data used in the 'data' argument")
    }
    else # attach needed attributes to given data
        data = recover_data(model, data = data)
    
    x = data[[var]]
    fcn = NULL   # differential
    if (is.null(x)) {
        fcn = var
        var = .all.vars(as.formula(paste("~",var)))
        if (length(var) > 1)
            stop("Can only support a function of one variable")
        else {
            x = data[[var]]
            if (is.null(x)) stop("Variable '", var, "' is not in the dataset")      
        }
    }
    rng = diff(range(x))
    if (delta.var == 0)  stop("Provide a nonzero value of 'delta.var'")
    
    RG = orig.rg = ref_grid(model, data = data, ...)
    
    max.degree = max(1, min(5, as.integer(max.degree + .1)))
    transform = match.arg(transform)
    if ((max.degree > 1) && (transform == "response") && hasName(RG@misc, "tran")) {
        max.degree = 1
        warning("Higher-degree trends are not supported with 'transform = 'response'.\n",
        "'max.degree' changed to 1")
    }
    
    # create a vector of delta values, such that a middle one has value 0
    delts = delta.var * (0:max.degree)
    idx.base = as.integer((2 + max.degree)/2)
    delts = delts - delts[idx.base]

    grid = lapply(delts, function(h) {
        g = RG@grid
        g[[var]] = g[[var]] + h
        g})
    if (!is.null(mr <- RG@roles$multresp)) {
        # RG@grid is expanded to mult levels, but we need to unexpand grid accordingly
        if (length(mr) > 0) {
            mri = which(grid[[1]][[mr]] == RG@levels[[mr]][1])
            grid = lapply(grid, function(.) .[mri, , drop = FALSE])
        }
    }
    options = list(...)
    linfct = lapply(grid, function(g) 
        emm_basis(model, attr(data, "terms"), RG@model.info$xlev, g, 
                  misc = attr(data, "misc"), options = options, ...)$X)
    
    if (!is.null(fcn)) { # need a different "h" when diff wrt a function
        tmp = sapply(grid, function(g) 
            eval(parse(text = fcn), envir = g))
        delta.var = apply(tmp, 1, function(.) mean(diff(.)))
    }
    
    newlf = numeric(0)
    h = 1
    for (i in 1:max.degree) { # successively difference linfct
        linfct = lapply(seq_along(linfct)[-1], function(j) linfct[[j]] - linfct[[j-1]])
        h = h * delta.var * i
        what = as.integer((length(linfct) + 1) / 2) # pick out one in the middle
        newlf = rbind(newlf, linfct[[what]] / h)
    }
    
    # Now replace linfct w/ difference quotient
    RG@linfct = newlf
    RG@roles$trend = var
    
    if(hasName(RG@misc, "tran")) {
        tran = RG@misc$tran
        if (is.list(tran)) tran = tran$name
        if (transform == "response") {
            prd = .est.se.df(orig.rg, do.se = FALSE)
            lnk = attr(prd, "link")
            deriv = lnk$mu.eta(prd[[1]])
            RG@linfct = diag(deriv) %*% RG@linfct
            RG@misc$initMesg = paste("Trends are obtained after back-transforming from the", tran, "scale")
        }
        else
            RG@misc$initMesg = paste("Trends are based on the", tran, "(transformed) scale")
    }
    
    # args for emmeans calls
    args = list(object = NULL, specs = specs, ...)
    args$at = args$cov.reduce = args$mult.levs = args$vcov. = args$data = args$trend = NULL
    
#    run_emm = TRUE
    if (max.degree > 1) {
        degnms = c("linear", "quadratic", "cubic", "quartic", "quintic")
        RG@grid$degree = degnms[1]
        g = RG@grid
        for (j in 2:max.degree) {
            g$degree = degnms[j]
            RG@grid = rbind(RG@grid, g)
        }
        RG@roles$predictors = c(RG@roles$predictors, "degree")
        RG@levels$degree = degnms[1:max.degree]
        chk = union(all.vars(specs), args$by)
        if (!("degree" %in% chk))
            args$by = c("degree", args$by)
    }
    RG@grid$.offset. = NULL   # offset never applies after differencing
    RG@misc$tran = RG@misc$tran.mult = NULL
    RG@misc$estName = estName
    RG@misc$methDesc = "emtrends"
    
    .save.ref_grid(RG)  # save in .Last.ref_grid, if enabled
    
    args$object = RG
    do.call("emmeans", args)
}

