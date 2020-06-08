##############################################################################
#    Copyright (c) 2012-2020 Russell V. Lenth                                #
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

# Reference grid code


# Change to cov.reduce specification: can be...
#     a function: is applied to all covariates
#     named list of functions: applied to those covariates (else mean is used)
#     TRUE - same as mean
#     FALSE - same as function(x) sort(unique(x))

#' Create a reference grid from a fitted model
#'
#' Using a fitted model object, determine a reference grid for which estimated
#' marginal means are defined. The resulting \code{ref_grid} object encapsulates
#' all the information needed to calculate EMMs and make inferences on them.
#'
#' To users, the \code{ref_grid} function itself is important because most of
#' its arguments are in effect arguments of \code{\link{emmeans}} and related
#' functions, in that those functions pass their \code{...} arguments to
#' \code{ref_grid}.
#'
#' The reference grid consists of combinations of independent variables over
#' which predictions are made. Estimated marginal means are defined as these
#' predictions, or marginal averages thereof. The grid is determined by first
#' reconstructing the data used in fitting the model (see
#' \code{\link{recover_data}}), or by using the \code{data.frame} provided in
#' \code{data}. The default reference grid is determined by the observed levels
#' of any factors, the ordered unique values of character-valued predictors, and
#' the results of \code{cov.reduce} for numeric predictors. These may be
#' overridden using \code{at}. See also the section below on
#' recovering/overriding model information.
#'
#'
#' @param object An object produced by a supported model-fitting function, such
#'   as \code{lm}. Many models are supported. See
#'   \href{../doc/models.html}{\code{vignette("models", "emmeans")}}.
#' @param at Optional named list of levels for the corresponding variables
#' @param cov.reduce A function, logical value, or formula; or a named list of
#'   these. Each covariate \emph{not} specified in \code{cov.keep} or \code{at}
#'   is reduced according to these specifications. See the section below on
#'   \dQuote{Using \code{cov.reduce} and \code{cov.keep}}.
#' @param cov.keep Character vector: names of covariates that are \emph{not}
#'   to be reduced; these are treated as factors and used in weighting calculations.
#'   \code{cov.keep} may also include integer value(s), and if so, the maximum
#'   of these is used to set a threshold such that any covariate having no more
#'   than that many unique values is automatically included in \code{cov.keep}.
#' @param mult.names Character value: the name(s) to give to the
#'   pseudo-factor(s) whose levels delineate the elements of a multivariate
#'   response. If this is provided, it overrides the default name(s) used for
#'   \code{class(object)} when it has a multivariate response (e.g., the default
#'   is \code{"rep.meas"} for \code{"mlm"} objects).
#' @param mult.levs A named list of levels for the dimensions of a multivariate
#'   response. If there is more than one element, the combinations of levels are
#'   used, in \code{\link{expand.grid}} order. The (total) number of levels must
#'   match the number of dimensions. If \code{mult.name} is specified, this
#'   argument is ignored.
#' @param options If non-\code{NULL}, a named \code{list} of arguments to pass
#'   to \code{\link{update.emmGrid}}, just after the object is constructed.
#' @param data A \code{data.frame} to use to obtain information about the
#'   predictors (e.g. factor levels). If missing, then
#'   \code{\link{recover_data}} is used to attempt to reconstruct the data. See
#'   the note with \code{\link{recover_data}} for an important precaution.
#' @param df Numeric value. This is equivalent to specifying \code{options(df =
#'   df)}. See \code{\link{update.emmGrid}}.
#' @param type Character value. If provided, this is saved as the
#'   \code{"predict.type"} setting. See \code{\link{update.emmGrid}} and the
#'   section below on prediction types and transformations.
#' @param transform Character, logical, or list. If non-missing, the reference
#'   grid is reconstructed via \code{\link{regrid}} with the given
#'   \code{transform} argument. See the section below on prediction types and
#'   transformations.
#' @param nesting If the model has nested fixed effects, this may be specified
#'   here via a character vector or named \code{list} specifying the nesting
#'   structure. Specifying \code{nesting} overrides any nesting structure that
#'   is automatically detected. See the section below on Recovering or Overriding 
#'   Model Information.
#' @param offset Numeric scalar value (if a vector, only the first element is
#'   used). This may be used to add an offset, or override offsets based on the
#'   model. A common usage would be to specify \code{offset = 0} for a Poisson
#'   regression model, so that predictions from the reference grid become rates
#'   relative to the offset that had been specified in the model.
#' @param sigma Numeric value to use for subsequent predictions or
#'   back-transformation bias adjustments. If not specified, we use
#'   \code{sigma(object)}, if available, and \code{NULL} otherwise.
#' @param ... Optional arguments passed to \code{\link{summary.emmGrid}},
#'   \code{\link{emm_basis}}, and
#'   \code{\link{recover_data}}, such as \code{params}, \code{vcov.} (see
#'   \bold{Covariance matrix} below), or options such as \code{mode} for
#'   specific model types (see \href{../doc/models.html}{vignette("models",
#'   "emmeans")}).
#'
#' @section Using \code{cov.reduce} and \code{cov.keep}: 
#'   The \code{cov.keep} argument was not available in \pkg{emmeans} versions
#'   1.4.1 and earlier. Any covariates named in this list are treated as if they
#'   are factors: all the unique levels are kept in the reference grid. The user
#'   may also specify an integer value, in which case any covariate having no more
#'   than that number of unique values is implicitly included in \code{cov.keep}.
#'   The default for \code{cove.keep} is set and retrieved via the 
#'   \code{\link{emm_options}} framework, and the system default is \code{"2"},
#'   meaning that covariates having only two unique values are automatically
#'   treated as two-level factors. See also the Note below on backward compatibility.
#'   
#'   There is a subtle distinction between including a covariate in \code{cov.keep}
#'   and specifying its values manually in \code{at}: Covariates included in 
#'   \code{cov.keep} are treated as factors for purposes of weighting, while
#'   specifying levels in \code{at} will not include the covariate in weighting.
#'   See the \code{mtcars.lm} example below for an illustration.
#'   
#'   \code{cov.reduce} may be a function,
#'   logical value, formula, or a named list of these.
#'   If a single function, it is applied to each covariate.
#'   If logical and \code{TRUE}, \code{mean} is used. If logical and
#'   \code{FALSE}, it is equivalent to including all covariates in
#'   \code{cov.keep}. Use of \samp{cov.reduce = FALSE} is inadvisable because it
#'   can result in a huge reference grid; it is far better to use
#'   \code{cov.keep}.
#'
#'   If a formula (which must be two-sided), then a model is fitted to that
#'   formula using \code{\link{lm}}; then in the reference grid, its response
#'   variable is set to the results of \code{\link{predict}} for that model,
#'   with the reference grid as \code{newdata}. (This is done \emph{after} the
#'   reference grid is determined.) A formula is appropriate here when you think
#'   experimental conditions affect the covariate as well as the response.
#'
#'   If \code{cov.reduce} is a named list, then the above criteria are used to
#'   determine what to do with covariates named in the list. (However, formula
#'   elements do not need to be named, as those names are determined from the
#'   formulas' left-hand sides.) Any unresolved covariates are reduced using
#'   \code{"mean"}.
#'
#'   Any \code{cov.reduce} of \code{cov.keep} specification for a covariate 
#'   also named in \code{at} is ignored.
#'   
#'   
#' @section Interdependent covariates: Care must be taken when covariate values
#'   depend on one another. For example, when a polynomial model was fitted
#'   using predictors \code{x}, \code{x2} (equal to \code{x^2}), and \code{x3}
#'   (equal to \code{x^3}), the reference grid will by default set \code{x2} and
#'   \code{x3} to their means, which is inconsistent. The user should instead
#'   use the \code{at} argument to set these to the square and cube of
#'   \code{mean(x)}. Better yet, fit the model using a formula involving
#'   \code{poly(x, 3)} or \code{I(x^2)} and \code{I(x^3)}; then there is only
#'   \code{x} appearing as a covariate; it will be set to its mean, and the
#'   model matrix will have the correct corresponding quadratic and cubic terms.
#'
#' @section Matrix covariates: Support for covariates that appear in the dataset
#'   as matrices is very limited. If the matrix has but one column, it is
#'   treated like an ordinary covariate. Otherwise, with more than one column,
#'   each column is reduced to a single reference value -- the result of
#'   applying \code{cov.reduce} to each column (averaged together if that
#'   produces more than one value); you may not specify values in \code{at}; and
#'   they are not treated as variables in the reference grid, except for
#'   purposes of obtaining predictions.
#'
#' @section Recovering or overriding model information: Ability to support a
#'   particular class of \code{object} depends on the existence of
#'   \code{recover_data} and \code{emm_basis} methods -- see
#'   \link{extending-emmeans} for details. The call
#'   \code{methods("recover_data")} will help identify these.
#'
#'   \bold{Data.} In certain models, (e.g., results of
#'   \code{\link[lme4]{glmer.nb}}), it is not possible to identify the original
#'   dataset. In such cases, we can work around this by setting \code{data}
#'   equal to the dataset used in fitting the model, or a suitable subset. Only
#'   the complete cases in \code{data} are used, so it may be necessary to
#'   exclude some unused variables. Using \code{data} can also help save
#'   computing, especially when the dataset is large. In any case, \code{data}
#'   must represent all factor levels used in fitting the model. It
#'   \emph{cannot} be used as an alternative to \code{at}. (Note: If there is a
#'   pattern of \code{NAs} that caused one or more factor levels to be excluded
#'   when fitting the model, then \code{data} should also exclude those levels.)
#'
#'   \bold{Covariance matrix.} By default, the variance-covariance matrix for
#'   the fixed effects is obtained from \code{object}, usually via its
#'   \code{\link{vcov}} method. However, the user may override this via a
#'   \code{vcov.} argument, specifying a matrix or a function. If a matrix, it
#'   must be square and of the same dimension and parameter order of the fixed
#'   effects. If a function, must return a suitable matrix when it is called
#'   with \code{object} as its only argument.
#'
#'   \bold{Nested factors.} Having a nesting structure affects marginal
#'   averaging in \code{emmeans} in that it is done separately for each level
#'   (or combination thereof) of the grouping factors. \code{ref_grid} tries to
#'   discern which factors are nested in other factors, but it is not always
#'   obvious, and if it misses some, the user must specify this structure via
#'   \code{nesting}; or later using \code{\link{update.emmGrid}}. The
#'   \code{nesting} argument may be a character vector, a named \code{list}, 
#'   or \code{NULL}.
#'   If a \code{list}, each name should be the name of a single factor in the
#'   grid, and its entry a character vector of the name(s) of its grouping
#'   factor(s). \code{nested} may also be a character value of the form
#'   \code{"factor1 \%in\% (factor2*factor3)"} (the parentheses are optional).
#'   If there is more than one such specification, they may be appended
#'   separated by commas, or as separate elements of a character vector. For
#'   example, these specifications are equivalent: \code{nesting = list(state =
#'   "country", city = c("state", "country")}, \code{nesting = "state \%in\%
#'   country, city \%in\% (state*country)"}, and \code{nesting = c("state \%in\%
#'   country", "city \%in\% state*country")}.
#'
#' @section Predictors with subscripts and data-set references: When the fitted
#'   model contains subscripts or explicit references to data sets, the
#'   reference grid may optionally be post-processed to simplify the variable
#'   names, depending on the \code{simplify.names} option (see
#'   \code{\link{emm_options}}), which by default is \code{TRUE}. For example,
#'   if the model formula is \code{data1$resp ~ data1$trt + data2[[3]] +
#'   data2[["cov"]]}, the simplified predictor names (for use, e.g., in the
#'   \code{specs} for \code{\link{emmeans}}) will be \code{trt},
#'   \code{data2[[3]]}, and \code{cov}. Numerical subscripts are not simplified;
#'   nor are variables having simplified names that coincide, such as if
#'   \code{data2$trt} were also in the model.
#'
#'   Please note that this simplification is performed \emph{after} the
#'   reference grid is constructed. Thus, non-simplified names must be used in
#'   the \code{at} argument (e.g., \code{at = list(`data2["cov"]` = 2:4)}.
#'
#'   If you don't want names simplified, use \code{emm_options(simplify.names =
#'   FALSE)}.
#'
#'
#'
#' @section Prediction types and transformations:
#'   Transformations can exist because of a link function in a generalized linear model, 
#'   or as a response transformation, or even both. In many cases, they are auto-detected,
#'   for example a model formula of the form \code{sqrt(y) ~ ...}. Even transformations
#'   containing multiplicative or additive constants, such as \code{2*sqrt(y + pi) ~ ...},
#'   are auto-detected. A response transformation of \code{y + 1 ~ ...} is \emph{not}
#'   auto-detected, but \code{I(y + 1) ~ ...} is interpreted as \code{identity(y + 1) ~ ...}.
#'   A warning is issued if it gets too complicated.
#'   Complex transformations like the Box-Cox transformation are not auto-detected; but see 
#'   the help page for \code{\link{make.tran}} for information on some advanced methods.
#'   
#'   There is a subtle difference
#'   between specifying \samp{type = "response"} and \samp{transform =
#'   "response"}. While the summary statistics for the grid itself are the same,
#'   subsequent use in \code{\link{emmeans}} will yield different results if
#'   there is a response transformation or link function. With \samp{type =
#'   "response"}, EMMs are computed by averaging together predictions on the
#'   \emph{linear-predictor} scale and then back-transforming to the response
#'   scale; while with \samp{transform = "response"}, the predictions are
#'   already on the response scale so that the EMMs will be the arithmetic means
#'   of those response-scale predictions. To add further to the possibilities,
#'   \emph{geometric} means of the response-scale predictions are obtainable via
#'   \samp{transform = "log", type = "response"}. See also the help page for 
#'   \code{\link{regrid}}.
#'
#' @section Side effect: The most recent result of \code{ref_grid}, whether
#'   called directly or indirectly via \code{\link{emmeans}},
#'   \code{\link{emtrends}}, or some other function that calls one of these, is
#'   saved in the user's environment as \code{.Last.ref_grid}. This facilitates
#'   checking what reference grid was used, or reusing the same reference grid
#'   for further calculations. This automatic saving is enabled by default, but
#'   may be disabled via \samp{emm_options(save.ref_grid = FALSE)}, and
#'   re-enabled by specifying \code{TRUE}.
#'
#' @return An object of the S4 class \code{"emmGrid"} (see
#'   \code{\link{emmGrid-class}}). These objects encapsulate everything needed
#'   to do calculations and inferences for estimated marginal means, and contain
#'   nothing that depends on the model-fitting procedure.
#'
#' @seealso Reference grids are of class \code{\link[=emmGrid-class]{emmGrid}},
#'   and several methods exist for them -- for example
#'   \code{\link{summary.emmGrid}}. Reference grids are fundamental to
#'   \code{\link{emmeans}}. Supported models are detailed in
#'   \href{../doc/models.html}{\code{vignette("models", "emmeans")}}.
#'   See \code{\link{update.emmGrid}} for details of arguments that can be in 
#'   \code{options} (or in \code{...}).
#'   
#' @note The system default for \code{cov.keep} causes models
#'   containing indicator variables to be handled differently than in
#'   \pkg{emmeans} version 1.4.1 or earlier. To replicate older
#'   analyses, change the default via 
#'   \samp{emm_options(cov.keep = character(0))}.
#'   
#' @note Some earlier versions of \pkg{emmeans} offer a \code{covnest} argument.
#'   This is now obsolete; if \code{covnest} is specified, it is harmlessly
#'   ignored. Cases where it was needed are now handled appropriately via the
#'   code associated with \code{cov.keep}.
#'
#' @export
#'
#' @examples
#' fiber.lm <- lm(strength ~ machine*diameter, data = fiber)
#' ref_grid(fiber.lm)
#' summary(.Last.ref_grid)
#'
#' ref_grid(fiber.lm, at = list(diameter = c(15, 25)))
#'
#' \dontrun{
#' # We could substitute the sandwich estimator vcovHAC(fiber.lm)
#' # as follows:
#' summary(ref_grid(fiber.lm, vcov. = sandwich::vcovHAC))
#' }
#'
#' # If we thought that the machines affect the diameters
#' # (admittedly not plausible in this example), then we should use:
#' ref_grid(fiber.lm, cov.reduce = diameter ~ machine)
#' 
#' ### Model with indicator variables as predictors:
#' mtcars.lm <- lm(mpg ~ disp + wt + vs * am, data = mtcars)
#' (rg.default <- ref_grid(mtcars.lm))
#' (rg.nokeep <- ref_grid(mtcars.lm, cov.keep = character(0)))
#' (rg.at <- ref_grid(mtcars.lm, at = list(vs = 0:1, am = 0:1)))
#' 
#' # Two of these have the same grid but different weights:
#' rg.default@grid
#' rg.at@grid
#'
#' # Multivariate example
#' MOats.lm = lm(yield ~ Block + Variety, data = MOats)
#' ref_grid(MOats.lm, mult.names = "nitro")
#' # Silly illustration of how to use 'mult.levs' to make comb's of two factors
#' ref_grid(MOats.lm, mult.levs = list(T=LETTERS[1:2], U=letters[1:2]))
#' 
#' # Using 'params'
#' require("splines")
#' my.knots = c(2.5, 3, 3.5)
#' mod = lm(Sepal.Length ~ Species * ns(Sepal.Width, knots = my.knots), data = iris)
#' ## my.knots is not a predictor, so need to name it in 'params'
#' ref_grid(mod, params = "my.knots") 
#' 
ref_grid <- function(object, at, cov.reduce = mean, cov.keep = get_emm_option("cov.keep"),
                     mult.names, mult.levs, 
                     options = get_emm_option("ref_grid"), data, df, type, 
                     transform, nesting, offset, sigma, ...) 
{
    ### transform = match.arg(transform)
    if (!missing(df)) {
        if(is.null(options)) options = list()
        options$df = df
    }
    
    if(missing(sigma)) { # Get 'sigma(object)' if available, else NULL
        sigma = suppressWarnings(try(stats::sigma(object), silent = TRUE))
        if (inherits(sigma, "try-error"))
            sigma = NULL
    }
    
    # recover the data
    if (missing(data)) {
        data = try(recover_data (object, data = NULL, ...))
        if (inherits(data, "try-error"))
            stop("Perhaps a 'data' or 'params' argument is needed")
    }
    else # attach needed attributes to given data
        data = recover_data(object, data = as.data.frame(data), ...)
    
    if(is.character(data)) # 'data' is in fact an error message
        stop(data)

    ## undocumented hook to use ref_grid as slave to get data
    if (!is.null(options$just.data))
        return(data)
    
    
    trms = attr(data, "terms")
    
    # find out if any variables are coerced to factors or vice versa
    coerced = .find.coerced(trms, data) # now list with members 'factors' and 'covariates'
    
    # convenience function
    sort.unique = function(x) sort(unique(x))
    
    # Get threshold for max #levels to keep a covariate
    if (is.null(cov.keep))
        cov.keep = character(0)
    cov.thresh = max(c(1, suppressWarnings(as.integer(cov.keep))), na.rm = TRUE)
    if (is.logical(cov.reduce)) {
        if (!cov.reduce)
            cov.thresh = 99 # big enough!
        cov.reduce = mean
    }
    
    
    
    
    # Ensure cov.reduce is a function or list thereof
    dep.x = list() # list of formulas to fit later
    fix.cr = function(cvr) {
        if (inherits(cvr, "formula")) {
            if (length(cvr) < 3)
                stop("Formulas in 'cov.reduce' must be two-sided")
            lhs = .all.vars(cvr)[1]
            dep.x[[lhs]] <<- cvr
            cvr = mean 
        }
        else if (!inherits(cvr, c("function","list")))
            stop("Invalid 'cov.reduce' argument")
        cvr
    }
    
    # IMPORTANT: following stmts may also affect x.dep
    if (is.list(cov.reduce))
        cov.reduce = lapply(cov.reduce, fix.cr)
    else
        cov.reduce = fix.cr(cov.reduce)
    
    # zap any formulas that are also in 'at'
    if (!missing(at))
        for (xnm in names(at)) dep.x[[xnm]] = NULL
    
    
    # local cov.reduce function that works with function or named list
    cr = function(x, nm) {
        if (is.function(cov.reduce))
            cov.reduce(x)
        else if (hasName(cov.reduce, nm))
            cov.reduce[[nm]](x)
        else
            mean(x)
    }
    
    # initialize empty lists
    ref.levels = matlevs = xlev = chrlev = list()
    
    for (nm in attr(data, "responses")) {
        y = data[[nm]]
        if (is.matrix(y))
            matlevs[[nm]] = apply(y, 2, mean)
        else
            ref.levels[[nm]] = mean(y)
    }
    
    for (nm in attr(data, "predictors")) {
        x = data[[nm]]
        if (is.matrix(x) && ncol(x) == 1)   # treat 1-column matrices as covariates
            x = as.numeric(x)
        
        # Save the original levels of factors, no matter what
        if (is.factor(x) && !(nm %in% coerced$covariates))
            xlev[[nm]] = levels(factor(x))
            # (applying factor drops any unused levels)
        else if (is.character(x)) # just like a factor
            xlev[[nm]] = sort(unique(x))
    
        # Now go thru and find reference levels...
        # mentioned in 'at' list but not coerced factor
        if (!(nm %in% coerced$factors) && !missing(at) && (hasName(at, nm)))
            ref.levels[[nm]] = at[[nm]]
        # factors not in 'at'
        else if (is.factor(x) && !(nm %in% coerced$covariates))
            ref.levels[[nm]] = levels(factor(x))
        else if (is.character(x) || is.logical(x))
            ref.levels[[nm]] = chrlev[[nm]] = sort.unique(x)
        # matrices
        else if (is.matrix(x)) {
            # Matrices -- reduce columns thereof, but don't add to baselevs
            matlevs[[nm]] = apply(x, 2, cr, nm)
            # if cov.reduce returns a vector, average its columns
            if (is.matrix(matlevs[[nm]]))
                matlevs[[nm]] = apply(matlevs[[nm]], 2, mean)
        }
        # covariate coerced, or not mentioned in 'at'
        else {
            # single numeric pred but coerced to a factor - use unique values
            # even if in 'at' list. We'll fix this up later
            if (nm %in% coerced$factors)            
                ref.levels[[nm]] = sort.unique(x)
            
            # Ordinary covariates - summarize based on cov.keep
            else {
                if ((length(uval <- sort.unique(x)) > cov.thresh) && !(nm %in% cov.keep))
                    ref.levels[[nm]] = cr(as.numeric(x), nm)
                else {
                    ref.levels[[nm]] = uval
                    cov.keep = c(cov.keep, nm)
                }
            }
        }
    }
    

    # Now create the reference grid
    grid = do.call(expand.grid, ref.levels)
    
    # undocumented hook to expand grid by increments of 'var' (needed by emtrends)
    if (!is.null(delts <- options$delts)) {
        var = options$var
        n.orig = nrow(grid) # remember how many rows we had
        grid = grid[rep(seq_len(n.orig), length(delts)), , drop = FALSE]
        options$var = options$delts = NULL
    }
    
    # add any matrices
    for (nm in names(matlevs)) {
        tmp = matrix(rep(matlevs[[nm]], each=nrow(grid)), nrow=nrow(grid))
        dimnames(tmp) = list(NULL, names(matlevs[[nm]]))
        grid[[nm]] = tmp
    }

    # resolve any covariate formulas
    for (xnm in names(dep.x)) {
        if (!all(.all.vars(dep.x[[xnm]]) %in% names(grid)))
            stop("Formulas in 'cov.reduce' must predict covariates actually in the model")
        xmod = lm(dep.x[[xnm]], data = data)
        grid[[xnm]] = predict(xmod, newdata = grid)
        ref.levels[[xnm]] = NULL
    }
    
    # finish-up our hook for expanding the grid
    if (!is.null(delts)) # add increments if any
        grid[[var]] = grid[[var]] + rep(delts, each = n.orig)
    
    if (!is.null(attr(data, "pass.it.on")))   # a hook needed by emm_basis.gamlss
        attr(object, "data") = data
    
    # we've added args `misc` and `options` so emm_basis methods can access and use these if they want
    basis = emm_basis(object, trms, xlev, grid, misc = attr(data, "misc"), options = options, ...)
    if(length(basis$bhat) != ncol(basis$X))
        stop("Something went wrong:\n",
             " Non-conformable elements in reference grid.",
             call. = TRUE)
    
    misc = basis$misc
    
    ### Figure out if there is a response transformation...
    # next stmt assumes that model formula is 1st argument (2nd element) in call.
    # if not, we probably get an error or something that isn't a formula
    # and it is silently ignored
    lhs = try(eval(as.formula(attr(data, "call")[[2]])[-3]), silent = TRUE)
    if (inherits(lhs, "formula")) { # response may be transformed
        tran = setdiff(.all.vars(lhs, functions = TRUE), c(.all.vars(lhs), "~", "cbind", "+", "-", "*", "/", "^", "%%", "%/%"))
        if(length(tran) > 0) {
            if (tran[1] == "scale") { # we'll try to set it up based on terms component
                pv = try(attr(terms(object), "predvars"), silent = TRUE)
                if (!inherits(pv, "try-error")) {
                    pv = c(lapply(pv, as.character), "foo") # make sure it isn't empty
                    scal = which(sapply(pv, function(x) x[1] == "scale"))
                    if(length(scal) > 0) {
                        par = as.numeric(pv[[scal[1]]][3:4]) 
                        tran = make.tran("scale", y = 0, center = par[1], scale = par[2])
                    }
                }
                if (is.character(tran)) { # didn't manage to find params
                    tran = NULL
                    message("NOTE: Unable to recover scale() parameters. See '? make.tran'")
                }
            }
            else if (tran[1] == "linkfun")
                tran = as.list(environment(trms))[c("linkfun","linkinv","mu.eta","valideta","name")]
            else {
                if (tran[1] == "I") 
                    tran = "identity"
                tran = paste(tran, collapse = ".")  
                # length > 1: Almost certainly unsupported, but facilitates a more informative error message
                const.warn = "There are unevaluated constants in the response formula\nAuto-detection of the response transformation may be incorrect"
                # Look for a multiplier, e.g. 2*sqrt(y)
                tst = strsplit(strsplit(as.character(lhs[2]), "\\(")[[1]][1], "\\*")[[1]]
                if(length(tst) > 1) {
                    mul = try(eval(parse(text = tst[1])), silent = TRUE)
                    if(!inherits(mul, "try-error")) {
                        misc$tran.mult = mul
                        tran = gsub("\\*\\.", "", tran)
                    }
                    else
                        warning(const.warn)
                }
                
                # look for added constant, e.g. log(y + 1)
                tst = strsplit(as.character(lhs[2]), "\\(|\\)|\\+")[[1]]
                if (length(tst) > 2) {
                    const = try(eval(parse(text = tst[3])), silent = TRUE)
                    if (!inherits(const, "try-error") && (length(tst) == 3))
                        misc$tran.offset = const
                    else
                        warning(const.warn)
                }
            }
            if(is.null(misc[["tran"]]))
                misc$tran = tran
            else
                misc$tran2 = tran
            misc$inv.lbl = "response"
        }
    }
    
    # Take care of multivariate response
    multresp = character(0) ### ??? was list()
    ylevs = misc$ylevs
    if(!is.null(ylevs)) { # have a multivariate situation
        if (missing(mult.levs))
            mult.levs = ylevs
        if(!missing(mult.names)) {
            k = seq_len(min(length(ylevs), length(mult.names)))
            names(mult.levs)[k] = mult.names[k]
        } 
        if (length(ylevs) > 1)
            ylevs = list(seq_len(prod(sapply(mult.levs, length))))
        
        k = prod(sapply(mult.levs, length))
        if (k != length(ylevs[[1]])) 
            stop("supplied 'mult.levs' is of different length ",
                 "than that of multivariate response")
        for (nm in names(mult.levs))
            ref.levels[[nm]] = mult.levs[[nm]]
        multresp = names(mult.levs)
        MF = do.call("expand.grid", mult.levs)
        grid = merge(grid, MF)
        
    }
            
    # add any matrices
    for (nm in names(matlevs))
        grid[[nm]] = matrix(rep(matlevs[[nm]], 
                                each=nrow(grid)), nrow=nrow(grid))
    
# Here's a complication. If a numeric predictor was coerced to a factor, we had to
# include all its levels in the reference grid, even if altered in 'at'
# Moreover, whatever levels are in 'at' must be a subset of the unique values
# So we now need to subset the rows of the grid and linfct based on 'at'

    problems = if (!missing(at)) 
        intersect(c(multresp, coerced$factors), names(at)) 
    else character(0)
    if (length(problems) > 0) {
        incl.flags = rep(TRUE, nrow(grid))
        for (nm in problems) {
            if (is.numeric(ref.levels[[nm]])) {
                at[[nm]] = round(at[[nm]], 3)
                ref.levels[[nm]] = round(ref.levels[[nm]], 3)
            }
            # get only "legal" levels
            at[[nm]] = at[[nm]][at[[nm]] %in% ref.levels[[nm]]]
            # Now which of those are left out?
            excl = setdiff(ref.levels[[nm]], at[[nm]])
            for (x in excl)
                incl.flags[grid[[nm]] == x] = FALSE
            ref.levels[[nm]] = at[[nm]]
        }
        if (!any(incl.flags))
            stop("Reference grid is empty due to mismatched levels in 'at'")
        grid = grid[incl.flags, , drop=FALSE]
        basis$X = basis$X[incl.flags, , drop=FALSE]
    }

    # Any offsets??? (misc$offset.mult might specify removing or reversing the offset)
    if (!missing(offset)) {  # For safety, we always treat it as scalar
        if (offset[1] != 0)
            grid[[".offset."]] = offset[1]
    }
    else if(!is.null(attr(trms,"offset"))) {
        om = 1
        if (!is.null(misc$offset.mult))
            om = misc$offset.mult
        if (any(om != 0))
            grid[[".offset."]] = om * .get.offset(trms, grid)
    }

    ### --- Determine weights for each grid point
    if (!hasName(data, "(weights)"))
        data[["(weights)"]] = 1
    cov.keep = intersect(unique(cov.keep), names(ref.levels))
    nms = union(union(union(names(xlev), names(chrlev)), coerced$factors), cov.keep)
    #### Old code...
    # if (!covnest)
    #     nms = union(union(names(xlev), names(chrlev)), coerced$factors) # only factors, no covariates or mult.resp
    # else
    #     nms = setdiff(names(ref.levels)[sapply(ref.levels, length) > 1], multresp) # all names (except multiv) for which there is > 1 level
    if (length(nms) == 0)
        wgt = rep(1, nrow(grid))  # all covariates; give each weight 1
    else {
        id = plyr::id(data[, nms, drop = FALSE], drop = TRUE)
        uid = !duplicated(id)
        key = do.call(paste, unname(data[uid, nms, drop = FALSE]))
        key = key[order(id[uid])]
        tgt = do.call(paste, unname(grid[, nms, drop = FALSE]))
        wgt = rep(0, nrow(grid))
        for (i in seq_along(key))
            wgt[tgt == key[i]] = sum(data[["(weights)"]][id==i])
    }
    grid[[".wgt."]] = wgt
    
    model.info = list(call = attr(data,"call"), terms = trms, xlev = xlev)
    # Detect any nesting structures
    nst = .find_nests(grid, trms, coerced$orig, ref.levels)
    if (length(nst) > 0)
        model.info$nesting = nst

    misc$is.new.rg = TRUE
    misc$ylevs = NULL # No longer needed
    misc$estName = "prediction"
    misc$estType = "prediction"
    misc$infer = c(FALSE,FALSE)
    misc$level = .95
    misc$adjust = "none"
    misc$famSize = nrow(grid)
    misc$avgd.over = character(0)
    misc$sigma = sigma

    post.beta = basis$post.beta
    if (is.null(post.beta))
        post.beta = matrix(NA)
    
    predictors = attr(data, "predictors")
    
    simp.tbl = environment(trms)$.simplify.names.
    if (! is.null(simp.tbl)) {
        names(grid) = .simplify.names(names(grid), simp.tbl)
        predictors = .simplify.names(predictors, simp.tbl)
        names(ref.levels) = .simplify.names(names(ref.levels), simp.tbl)
        if (!is.null(post.beta)) names(post.beta) = .simplify.names(names(post.beta), simp.tbl)
        if (!is.null(model.info$nesting)) {
            model.info$nesting = lapply(model.info$nesting, .simplify.names, simp.tbl)
            names(model.info$nesting) = .simplify.names(names(model.info$nesting), simp.tbl)
        }
        
        environment(trms)$.simplify.names. = NULL
    }
    
    result = new("emmGrid",
         model.info = model.info,
         roles = list(predictors = predictors, 
                      responses = attr(data, "responses"), 
                      multresp = multresp),
         grid = grid, levels = ref.levels, matlevs = matlevs,
         linfct = basis$X, bhat = basis$bhat, nbasis = basis$nbasis, V = basis$V,
         dffun = basis$dffun, dfargs = basis$dfargs, 
         misc = misc, post.beta = post.beta)
        
    if (!missing(type)) {
        if (is.null(options)) options = list()
        options$predict.type = type
    }
    
    if (!missing(nesting)) {
        result@model.info$nesting = lst = .parse_nest(nesting)
        if(!is.null(lst)) {
            nms = union(names(lst), unlist(lst))
            if(!all(nms %in% names(result@grid)))
                stop("Nonexistent variables specified in 'nesting'")
            result@misc$display = .find.nonempty.nests(result, nms)
        }
    }
        
    else if (!is.null(nst <- result@model.info$nesting)) {
        result@misc$display = .find.nonempty.nests(result)
        if (get_emm_option("msg.nesting"))
            message("NOTE: A nesting structure was detected in the ",
                    "fitted model:\n    ", .fmt.nest(nst))
    }

    result = .update.options(result, options, ...)

    if(!is.null(hook <- misc$postGridHook)) {
        if (is.character(hook))
            hook = get(hook)
        result@misc$postGridHook = NULL
        result = hook(result, ...)
    }
    if(!missing(transform))
        result = regrid(result, transform = transform, sigma = sigma, ...)
    
    .save.ref_grid(result)
    result
}


#### End of ref_grid ------------------------------------------

# local utility to save each newly constructed ref_grid, if enabled
# Goes into global environment unless .Last.ref_grid is found further up
.save.ref_grid = function(object) {
    if (is.logical(isnewrg <- object@misc$is.new.rg))
        if(isnewrg && get_emm_option("save.ref_grid"))
            assign(".Last.ref_grid", object, inherits = TRUE)
}



# This function figures out which covariates in a model 
# have been coerced to factors. And also which factors have been coerced
# to be covariates
.find.coerced = function(trms, data) {
    if (ncol(data) == 0) 
        return(list(factors = integer(0), covariates = integer(0)))
    isfac = sapply(data, function(x) inherits(x, "factor"))
    
    # Character vectors of factors and covariates in the data...
    facs.d = names(data)[isfac]
    covs.d = names(data)[!isfac]
    
    lbls = attr(trms, "term.labels")
    M = model.frame(trms, utils::head(data, 2)) #### just need a couple rows
    isfac = sapply(M, function(x) inherits(x, "factor"))
    
    # Character vector of terms in the model frame that are factors ...
    facs.m = names(M)[as.logical(isfac)]
    covs.m = setdiff(names(M), facs.m)
    
    # Exclude the terms that are already factors
    # What's left will be things like "factor(dose)", "interact(dose,treat)", etc
    # we're saving these in orig
    orig = cfac = setdiff(facs.m, facs.d)
    if(length(cfac) != 0) {
        cvars = lapply(cfac, function(x) .all.vars(stats::reformulate(x))) # Strip off the function calls
        cfac = intersect(unique(unlist(cvars)), covs.d) # Exclude any variables that are already factors
    }
    
    # Do same with covariates
    ccov = setdiff(covs.m, covs.d)
    orig = c(orig, ccov)
    if(length(ccov) > 0) {
        cvars = lapply(ccov, function(x) .all.vars(stats::reformulate(x)))
        ccov = intersect(unique(unlist(cvars)), facs.d)
    }
    
    list(factors = cfac, covariates = ccov, orig = orig)
}



