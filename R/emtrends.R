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
#' The function works by constructing reference grids for \code{object} with 
#' various values of \code{var}, and then calculating difference quotients of predictions
#' from those reference grids. Finally, \code{\link{emmeans}} is called with
#' the given \code{specs}, thus computing marginal averages as needed of
#' the difference quotients. Any \code{...} arguments are passed to the
#' \code{ref_grid} and \code{\link{emmeans}}; examples of such optional
#' arguments include optional arguments (often \code{mode}) that apply to
#' specific models; \code{ref_grid} options such as \code{data}, \code{at},
#' \code{cov.reduce}, \code{mult.names}, \code{nesting}, or \code{transform};
#' and \code{emmeans} options such as \code{weights} (but please avoid
#' \code{trend} or \code{offset}.
#' 
#' 
#' @param object A supported model object (\emph{not} a reference grid)
#' @param specs Specifications for what marginal trends are desired -- as in
#'   \code{\link{emmeans}}. If \code{specs} is missing or \code{NULL},
#'   \code{emmeans} is not run and the reference grid for specified trends
#'   is returned.
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
#' @param max.degree Integer value. The maximum degree of trends to compute (this
#'   is capped at 5). If greater than 1, an additional factor \code{degree} is
#'   added to the grid, with corresponding numerical derivatives of orders
#'   \code{1, 2, ..., max.degree} as the estimates.
#' @param ... Additional arguments passed to \code{\link{ref_grid}} or 
#'   \code{\link{emmeans}} as appropriate. See Details.
#'
#' @section Generalizations:
#' Instead of a single predictor, the user may specify some monotone function of
#' one variable, e.g., \code{var = "log(dose)"}. If so, the chain rule is
#' applied. Note that, in this example, if \code{object} contains
#' \code{log(dose)} as a predictor, we will be comparing the slopes estimated by
#' that model, whereas specifying \code{var = "dose"} would perform a
#' transformation of those slopes, making the predicted trends vary depending on
#' \code{dose}.
#' 
#' @return An \code{emmGrid} or \code{emm_list} object, according to \code{specs}.
#' See \code{\link{emmeans}} for more details on when a list is returned.
#' 
#' @seealso \code{link{emmeans}}, \code{\link{ref_grid}}
#' 
#' @note
#' In earlier versions of \code{emtrends}, the first argument was named
#' \code{model} rather than \code{object}. (The name was changed because of
#' potential mis-matching with a \code{mode} argument, which is an option for
#' several types of models.) For backward compatibility, \code{model} still works
#' \emph{provided all arguments are named}.
#' 
#' @note
#' It is important to understand that trends computed by \code{emtrends} are
#' \emph{not} equivalent to polynomial contrasts in a parallel model where
#' \code{var} is regarded as a factor. That is because the model \code{object}
#' here is assumed to fit a smooth function of \code{var}, and the estimated
#' trends reflect \emph{local} behavior at particular value(s) of \code{var};
#' whereas when \code{var} is modeled as a factor and polynomial contrasts are
#' computed, those contrasts represent the \emph{global} pattern of changes over
#' \emph{all} levels of \code{var}. 
#' 
#' See the \code{pigs.poly} and \code{pigs.fact} examples below for an
#' illustration. The linear and quadratic trends depend on the value of 
#' \code{percent}, but the cubic trend is constant (because that is true of
#' a cubic polynomial, which is the underlying model). The cubic contrast
#' in the factorial model has the same P value as for the cubic trend,
#' again because the cubic trend is the same everywhere.
#' 
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
#' # Obtaining a reference grid
#' mtcars.lm <- lm(mpg ~ poly(disp, degree = 2) * (factor(cyl) + factor(am)), data = mtcars)
#' 
#' Center trends at mean disp for each no. of cylinders
#' mtcTrends.rg <- emtrends(mtcars.lm, var = "disp", 
#'                           cov.reduce = disp ~ factor(cyl))
#' summary(mtcTrends.rg)  # estimated trends at grid nodes
#' emmeans(mtcTrends.rg, "am", weights = "prop")
#' 
#' 
#' ### Higher-degree trends ...
#' 
#' pigs.poly <- lm(conc ~ poly(percent, degree = 3), data = pigs)
#' emt <- emtrends(pigs.poly, ~ degree | percent, "percent", max.degree = 3,
#'                 at = list(percent = c(9, 13.5, 18)))
#'        # note: 'degree' is an extra factor created by 'emtrends'
#'        
#' summary(emt, infer = c(TRUE, TRUE))
#' 
#' # Compare above results with poly contrasts when 'percent' is modeled as a factor ...
#' pigs.fact <- lm(conc ~ factor(percent), data = pigs)
#' emm <- emmeans(pigs.fact, "percent")
#' 
#' contrast(emm, "poly")
#' # Some P values are comparable, some aren't! See Note in documentation
emtrends = function(object, specs, var, delta.var=.001*rng,
                    max.degree = 1, ...) {
    estName = paste(var, "trend", sep=".") # Do now as I may replace var later
    cl = match.call()
    cl[[1]] = quote(ref_grid)
    cl$var = cl$specs = NULL
    if (is.null(cl$options))
        cl$options = list()
    
    # backward compatibility for when 1st argument was "model"
    if(missing(object) && ("model" %in% names(cl))) {
        names(cl)[names(cl) == "model"] = "object"
        object = eval(cl$object)
    }
    
    # Get data via hook in ref_grid
    cl$options$just.data = TRUE
    data = eval(cl)
    cl$options$just.data = NULL
    
    x = data[[var]]
    fcn = NULL   # differential
    if (is.null(x)) {
        fcn = var
        var = .all.vars(as.formula(paste("~", var)))
        if (length(var) > 1)
            stop("Can only support a function of one variable")
        else {
            x = data[[var]]
            if (is.null(x)) stop("Variable '", var, "' is not in the dataset")      
        }
    }
    rng = diff(range(x))
    if (delta.var == 0)  stop("Provide a nonzero value of 'delta.var'")
    
    max.degree = max(1, min(5, as.integer(max.degree + .1)))

        # create a vector of delta values, such that a middle one has value 0
    delts = delta.var * (0:max.degree)
    idx.base = as.integer((2 + max.degree)/2)
    delts = delts - delts[idx.base]
    
    # set up call for ref_grid
    cl$data = quote(data)
    cl$options$var = var
    cl$options$delts = delts   # ref_grid hook -- expand grid by these increments
    bigRG = eval(cl)
    
    var.subs = as.list(as.data.frame(matrix(seq_len(nrow(bigRG@grid)), ncol = length(delts))))
    RG = orig.rg = bigRG[var.subs[[idx.base]]]  # the subset that corresponds to reference values
    row.names(RG@grid) = seq_along(RG@grid[[1]])

    linfct = lapply(seq_along(delts), function(i) bigRG@linfct[var.subs[[i]], , drop = FALSE])
    
    if (!is.null(fcn)) { # need a different "h" when diff wrt a function
        tmp = sapply(seq_along(delts), function(i)
            eval(parse(text = fcn), envir = bigRG@grid[var.subs[[i]], , drop = FALSE]))
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
    
    # # args for emmeans calls
    # args = list(object = NULL, specs = specs, ...)
    # args$at = args$cov.reduce = args$mult.levs = args$vcov. = args$data = args$trend = args$transform = NULL
    
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
    }
    RG@grid$.offset. = NULL   # offset never applies after differencing
    RG@misc$tran = RG@misc$tran.mult = NULL
    RG@misc$estName = estName
    RG@misc$methDesc = "emtrends"
    
    .save.ref_grid(RG)  # save in .Last.ref_grid, if enabled
    
    if (missing(specs) || is.null(specs))
        return (RG)
    
    # args for emmeans calls
    args = list(object = NULL, specs = specs, ...)
    args$at = args$cov.reduce = args$mult.levs = args$vcov. = args$data = args$trend = args$transform = NULL
    if(max.degree > 1) {
        chk = union(all.vars(specs), args$by)
        if (!("degree" %in% chk))
            args$by = c("degree", args$by)
    }
    
    
    args$object = RG
    do.call("emmeans", args)
}

