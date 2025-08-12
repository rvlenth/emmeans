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
#' which predictions are made. 
#' Estimated marginal means are defined as these
#' predictions, or marginal averages thereof. The grid is determined by first
#' reconstructing the data used in fitting the model (see
#' \code{\link{recover_data}}), or by using the \code{data.frame} provided in
#' \code{data}. 
#' 
#' By \dQuote{independent variables,} we mean (in most cases) the results of
#' \code{all.vars()} applied to the fixed-effects part of the right-hand side of
#' the model formula. Any random effects are excluded. Thus, if the model
#' formula in an \code{lme4::lmer} call is \code{yield ~ fert +
#' seed*density + log(rain) + (1|block/plot)}, the independent variables are
#' \code{fert}, \code{seed}, \code{density}, and \code{rain} (not
#' \code{log(rain)}). In multivariate models, the dimension of the multivariate
#' response is also considered an independent variable.
#' 
#' 
#' The default reference grid is determined by the observed levels
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
#' @param regrid Character, logical, or list. If non-missing, the reference
#'   grid is reconstructed via \code{\link{regrid}} with the argument
#'   \code{transform = regrid}. See the section below on prediction types and
#'   transformations. \emph{Note:} This argument was named \code{transform} in
#'   version 1.7.2 and earlier. For compatibility with old code, \code{transform}
#'   is still accepted if found among \code{...}, 
#'   as long as it doesn't match \code{tran}.
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
#'   Note: This applies only when the family is \code{"gaussian"}; for other families,
#'   \code{sigma} is set to \code{NA} and cannot be overridden.
#' @param counterfactuals \code{counterfactuals} specifies character
#'   names of counterfactual factors. If this is non-missing, a reference grid
#'   is created consisting of combinations of counterfactual levels 
#'   and the actual levels of those same factors.
#'   This grid is always converted to the response transformation scale
#'   and averaged over the actual factor levels. See the section below 
#'   on counterfactuals. 
#' @param nuisance,non.nuisance,wt.nuis If \code{nuisance} is a vector of predictor names,
#'   those predictors are omitted from the reference grid. Instead, the result 
#'   will be as if we had averaged over the levels of those factors, with either 
#'   equal or proportional weights as specified in \code{wt.nuis} (see the 
#'   \code{weights} argument in \code{\link{emmeans}}). The factors in 
#'   \code{nuisance} must not interact with other factors, not even other
#'   nuisance factors. Specifying nuisance factors can save considerable
#'   storage and computation time, and help avoid exceeding the maximum
#'   reference-grid size (\code{get_emm_option("rg.limit")}). 
#'   (\emph{Note:} For certain models where the \code{emm_basis} method returns a
#'   re-gridded parameterization, nuisance factors cannot be used, and an error
#'   is thrown.)
#' @param rg.limit Integer limit on the number of reference-grid rows to allow
#'   (checked before any multivariate responses are included).
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
#'   The default for \code{cov.keep} is set and retrieved via the 
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
#'   To allow for situations where a simple \code{lm()} call as described above won't
#'   be adequate, a formula of the form \code{ext ~ fcnname} is also supported,
#'   where the left-hand side may be \code{ext}, \code{extern}, or
#'   \code{external} (and must \emph{not} be a predictor name) and the
#'   right-hand side is the name of an existing function. The function is called
#'   with one argument, a data frame with columns for each variable in the
#'   reference grid. The function is expected to use that frame as new data to
#'   be used to obtain predictions for one or more models; and it should return
#'   a named list or data frame with replacement values for one or more of the
#'   covariates.
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
#'   with arguments \code{(object, ...)}. Be careful with possible 
#'   unintended conflicts with arguments in \code{...}; for example, 
#'   \code{sandwich::vcovHAC()} has optional arguments \code{adjust} and \code{weights}
#'   that may be intended for \code{emmeans()} but will also be passed to \code{vcov.()}.
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
#'   between specifying \samp{type = "response"} and \samp{regrid =
#'   "response"}. While the summary statistics for the grid itself are the same,
#'   subsequent use in \code{\link{emmeans}} will yield different results if
#'   there is a response transformation or link function. With \samp{type =
#'   "response"}, EMMs are computed by averaging together predictions on the
#'   \emph{linear-predictor} scale and then back-transforming to the response
#'   scale; while with \samp{regrid = "response"}, the predictions are
#'   already on the response scale so that the EMMs will be the arithmetic means
#'   of those response-scale predictions. To add further to the possibilities,
#'   \emph{geometric} means of the response-scale predictions are obtainable via
#'   \samp{regrid = "log", type = "response"}. See also the help page for 
#'   \code{\link{regrid}}.
#'   
#'   \emph{Order-of-processing issues:} 
#'   The \code{regrid} argument, if present, is acted on immediately after the reference 
#'   grid is constructed, while some of the \code{...} arguments may be used to
#'   update the object at the very end. Thus, code like
#'   \code{ref_grid(mod, tran = "sqrt", regrid = "response")} will not work correctly
#'   if the intention was to specify the response transformation, because the re-grid 
#'   is done \emph{before} it processes \code{tran = "sqrt"}. To get the intended
#'   result, do
#'   \code{regrid(ref_grid(mod, tran = "sqrt"), transform = "response")}.
#'
#' @section Counterfactuals:
#'   If \code{counterfactuals} is specified, the rows of the entire dataset
#'   become part of the reference grid, and the other reference levels are
#'   confined to those named in \code{counterfactuals}. In this type of analysis
#'   (called G-computation), we substitute (or impute) each combination of counterfactual
#'   levels into the entire dataset. Thus, predictions from this grid are those
#'   of each observation under each of the counterfactual levels. For this to
#'   make sense, we require an assumption of exchangeability of these levels.
#'   
#'   This grid is always converted to the response scale, as G-computation on
#'   the linear-predictor scale produces the same results as ordinary weighted EMMs.
#'   If we have counterfactual factors \code{A, B}, the reference grid also includes
#'   factors \code{actual_A, actual_B} which are used to track which observations
#'   originally had the \code{A, B} levels before they were changed by the 
#'   counterfactuals code. We average the response-scale predictions for each
#'   combination of actual levels and imputed levels (and multivariate levels, 
#'   if any). See additional discussion of how \code{\link{emmeans}} handles
#'   counterfactuals under that documentation.
#'   
#'   Currently, counterfactuals are not supported when the reference grid
#'   requires post-processing (e.g., ordinal models with \code{mode = "prob"}).
#'   Cases where we have nested factor levels can be complicated if mixed-in with counterfactuals, 
#'   and we make no guarantees. 
#'   Note that past implementations included arguments \code{wt.counter} and \code{avg.counter},
#'   which are now deprecated and are just ignored if specified.
#' 
#' @section Optional side effect: If the \code{save.ref_grid} option is set to
#'   \code{TRUE} (see \code{\link{emm_options}}),
#'   The most recent result of \code{ref_grid}, whether
#'   called directly or indirectly via \code{\link{emmeans}},
#'   \code{\link{emtrends}}, or some other function that calls one of these, is
#'   saved in the user's environment as \code{.Last.ref_grid}. This facilitates
#'   checking what reference grid was used, or reusing the same reference grid
#'   for further calculations. This automatic saving is disabled by default, but
#'   may be enabled via \samp{emm_options(save.ref_grid = TRUE)}.
#'   
#' 
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
#' ### Using cov.reduce formulas...
#' # Above suggests we can vary disp indep. of other factors - unrealistic
#' rg.alt <- ref_grid(mtcars.lm, at = list(wt = c(2.5, 3, 3.5)),
#'     cov.reduce = disp ~ vs * wt)
#' rg.alt@grid
#' 
#' # Alternative to above where we model sqrt(disp)
#' disp.mod <- lm(sqrt(disp) ~ vs * wt, data = mtcars)
#' disp.fun <- function(dat)
#'     list(disp = predict(disp.mod, newdata = dat)^2)
#' rg.alt2 <- ref_grid(mtcars.lm, at = list(wt = c(2.5, 3, 3.5)),
#'     cov.reduce = external ~ disp.fun)
#' rg.alt2@grid
#' 
#'
#' # Multivariate example
#' MOats.lm = lm(yield ~ Block + Variety, data = MOats)
#' ref_grid(MOats.lm, mult.names = "nitro")
#' # Silly illustration of how to use 'mult.levs' to make comb's of two factors
#' ref_grid(MOats.lm, mult.levs = list(T=LETTERS[1:2], U=letters[1:2]))
#' 
#' # Comparing estimates with and without counterfactuals
#' neuralgia.glm <- glm(Pain ~ Treatment + Sex + Age + Duration, 
#'                      family = binomial(), data = neuralgia)
#' emmeans(neuralgia.glm, "Treatment", type = "response")
#' 
#' emmeans(neuralgia.glm, "Treatment", counterfactuals = "Treatment")
#' 
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
                     regrid, nesting, offset, sigma, 
                     counterfactuals, ## wt.counter, avg.counter = TRUE,
                     nuisance = character(0), non.nuisance, wt.nuis = "equal",
                      rg.limit = get_emm_option("rg.limit"), ...) 
{
    if(!missing(counterfactuals)) { # route this to a different routine
        cl = match.call()
        cl[[1]] = as.name(".cf.refgrid") # internal function for counterfactuals
        return(eval(cl))
    }
    
    # hack to ignore 'tran' in dots arguments and interpret 'transform' as `regrid` :
    .foo = function(t,tr,tra,tran, transform = NULL, ...) transform
    .bar = .foo(...)
    if (!is.null(.bar)) {
        regrid = .bar
        message("In 'ref_grid()', use 'regrid = ...' rather than 'transform = ...' ",
        "to avoid this message.")
    }
    
    if (!missing(df)) {
        if(is.null(options)) options = list()
        options$df = df
    }
    
    # recover the data
    if (missing(data)) {
        data = try(recover_data (object, data = NULL, ...))
        if (inherits(data, "try-error"))
            stop("Perhaps a 'data' or 'params' argument is needed")
    }
    # note if options$delts non-null, we called 2nd time from emtrends and data is ok as-is
    else if (is.null(options$delts)) # attach needed attributes to given data
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
            xlev[[nm]] = levels(.chk.fac(x))
            # (drops any used levels unless allow.na.levs option is TRUE)
        else if (is.character(x)) # just like a factor
            xlev[[nm]] = sort(unique(x))
    
        # Now go thru and find reference levels...
        # mentioned in 'at' list but not coerced factor
        if (!(nm %in% coerced$factors) && !missing(at) && (hasName(at, nm)))
            ref.levels[[nm]] = at[[nm]]
        # factors not in 'at'
        else if (is.factor(x) && !(nm %in% coerced$covariates))
            ref.levels[[nm]] = levels(.chk.fac(x))
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
    

    if (!missing(non.nuisance))
        nuisance = setdiff(names(ref.levels), non.nuisance)

    # Now create the reference grid
    if(no.nuis <- (length(nuisance) == 0)) {
        .check.grid(ref.levels, rg.limit)
        grid = do.call(expand.grid, ref.levels)
    }
    else {
        nuis.info = .setup.nuis(nuisance, ref.levels, trms, rg.limit)
        grid = nuis.info$grid
    }
    
    # undocumented hook to expand grid by increments of 'var' (needed by emtrends)
    if (!is.null(delts <- options$delts)) {
        var = options$var
        n.orig = nrow(grid) # remember how many rows we had
        grid = grid[rep(seq_len(n.orig), length(delts)), , drop = FALSE]
        options$var = options$delts = NULL
        # grid[[var]] = grid[[var]] + rep(delts, each = n.orig)
        #     (we have to wait until after covariate calcs to do this)
    }
    
    # add any matrices
    for (nm in names(matlevs)) {
        tmp = matrix(rep(matlevs[[nm]], each=nrow(grid)), nrow=nrow(grid))
        dimnames(tmp) = list(NULL, names(matlevs[[nm]]))
        grid[[nm]] = tmp
    }

    # resolve any covariate formulas
    for (xnm in names(dep.x)) {
        if ((xnm %in% c("ext", "extern", "external")) && !(xnm %in% names(grid))) { # use external fcn
            fun = get(as.character(dep.x[[xnm]][[3]]), inherits = TRUE)
            rslts = fun(grid)   # should be some thing that supports names() and [[]], e.g., a list or d.f.
            for (nm in  intersect(names(rslts), names(grid))) {
                grid[[nm]] = rslts[[nm]]
                ref.levels[[nm]] = NULL
            }
        }
        else if (!all(.all.vars(dep.x[[xnm]]) %in% names(grid)))
            stop("Formulas in 'cov.reduce' must predict covariates actually in the model")
        else { # Use lm() to predict this covariate
            xmod = lm(dep.x[[xnm]], data = data)
            grid[[xnm]] = predict(xmod, newdata = grid)
            ref.levels[[xnm]] = NULL
        }
    }
    
    # finish-up our hook for expanding the grid
    if (!is.null(delts)) # add increments if any
        grid[[var]] = grid[[var]] + rep(delts, each = n.orig)
    
    if (!is.null(attr(data, "pass.it.on")))   # a hook needed by emm_basis.gamlss et al.
        attr(object, "data") = data
    
    ###!! Prevent a warning like in https://stackoverflow.com/questions/68969384/emmeans-warning-in-model-frame-defaultformula-data-data-variable-gr/68990172#68990172
    xl = xlev
    modnm = rownames(attr(trms, "factors"))
    chk = sapply(modnm, function(mn) mn %in% names(xl))
    for(i in which(!chk)) { # replace names in xl - e.g., as.factor(trt) where trt already a factor
        fn = all.vars(reformulate(modnm[i]))
        if (length(fn) == 1)
            names(xl)[names(xl) == fn] = modnm[i]
    }
    ###!! If we remove this code, also need to change xl back to xlev in 'basis =' call below
    
    # we've added args `misc` and `options` so emm_basis methods can access and use these if they want
    basis = emm_basis(object, trms, xl, grid, misc = attr(data, "misc"), options = options, ...)

    environment(basis$dffun) = baseenv()   # releases unnecessary storage
    if(length(basis$bhat) != ncol(basis$X))
        stop("Something went wrong:\n",
             " Non-conformable elements in reference grid.",
             call. = TRUE)
    
    if(!no.nuis) {
        basis = .basis.nuis(basis, nuis.info, wt.nuis, ref.levels, data, grid, ref.levels)
        grid = basis$grid
        nuisance = ref.levels[nuis.info$nuis] # now nuisance has the levels info
        ref.levels = basis$ref.levels
    }
    
    misc = basis$misc
    
    ### Figure out if there is a response transformation...
    # next stmt assumes that model formula is 1st argument (2nd element) in call.
    # if not, we probably get an error or something that isn't a formula
    # and it is silently ignored
    frm = try(formula(eval(attr(data, "call")[[2]], environment(trms))), silent = TRUE)
    if (inherits(frm, "formula")) { # response may be transformed
        lhs = if(length(frm) == 2) NULL
              else frm[-3]
        tran = setdiff(.all.vars(lhs, functions = TRUE), c(.all.vars(lhs), "~", "cbind", "+", "-", "*", "/", "^", "%%", "%/%"))
        if(length(tran) > 0) {
            # we are now supporting scale() as well as some functions from datawizard
            if (tran[1] %in% c("scale", "center", "centre", "standardize", "standardise")) { # we'll try to set it up based on terms component
                pv = try(attr(terms(object), "predvars"), silent = TRUE)
                if (!inherits(pv, "try-error") && !is.null(pv)) {
                    scal = which(sapply(c(sapply(pv, as.character), "foo"), 
                                        function(x) x[1]) == tran[1])
                    if(length(scal) > 0) {
                        pv = pv[[scal[1]]]
                        ctr = ifelse(is.null(pv$center), 0, ifelse(pv$center, pv$center, 0))
                        scl = ifelse(is.null(pv$scale), 1, ifelse(pv$scale, pv$scale, 1))
                        tran = make.tran("scale", y = 0, center = ctr, scale = scl)
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
                dig = 3 - log10(max(abs(ref.levels[[nm]])))
                at[[nm]] = round(at[[nm]], digits = dig)
                ref.levels[[nm]] = round(ref.levels[[nm]], digits = dig)
                grid[[nm]] = round(grid[[nm]], digits = dig)
            }
            # get only "legal" levels
            at[[nm]] = ref.levels[[nm]] = at[[nm]][at[[nm]] %in% ref.levels[[nm]]]
            rows = numeric(0)
            for(x in at[[nm]])
                rows = c(rows, which(grid[[nm]] == x))
            grid = grid[rows, , drop = FALSE]
            grid[[nm]] = factor(grid[[nm]], levels = at[[nm]])
            basis$X = basis$X[rows, , drop = FALSE]
        }
    }

    # Any offsets??? (misc$offset.mult might specify removing or reversing the offset)
    om = ifelse(is.null(misc$offset.mult), 1, misc$offset.mult)
    oval = 0
    if (!missing(offset)) {  # For safety, we always treat it as scalar
        if (offset[1] != 0)
            oval = offset[1]
        if(".static.offset." %in% names(grid))
            grid$.static.offset. = ref.levels$.static.offset. = oval
    }
    else {
        if(".static.offset." %in% names(grid)) { # static offset is available
            oval = om * grid[[".static.offset."]]
        }
        #else implied
        if(!is.null(attr(trms,"offset"))) {
            if (any(om != 0))
                oval = om * (oval + .get.offset(trms, grid))
        }
    }
    if(any(oval != 0))
        grid[[".offset."]] = oval
    
    
    ### --- Determine weights for each grid point
    if (!hasName(data, "(weights)"))
        data[["(weights)"]] = 1
    cov.keep = intersect(unique(cov.keep), names(ref.levels))
    nms = union(union(union(names(xlev), names(chrlev)), coerced$factors), cov.keep)
    nms = intersect(nms, names(grid))
    #### Old code...
    # if (!covnest)
    #     nms = union(union(names(xlev), names(chrlev)), coerced$factors) # only factors, no covariates or mult.resp
    # else
    #     nms = setdiff(names(ref.levels)[sapply(ref.levels, length) > 1], multresp) # all names (except multiv) for which there is > 1 level
    if (length(nms) == 0)
        wgt = rep(1, nrow(grid))  # all covariates; give each weight 1
    else {
        id = .my.id(data[, nms, drop = FALSE])
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
    if (!is.null(mm <- basis$model.matrix)) { # submodel support
        attr(mm, "factors") = .smpFT(trms)
        model.info$model.matrix = mm  
    }
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
    if(is.null(misc$avgd.over))
       misc$avgd.over = character(0)
    
    # emm_basis may have provided a sigma. If not, use sigma value from model or args
    if(is.null(misc$sigma) && missing(sigma)) { # Get 'sigma(object)' if available, else NULL
        sigma = suppressWarnings(try(stats::sigma(object), silent = TRUE))
        if (inherits(sigma, "try-error"))
            sigma = NULL
        misc$sigma = sigma
    }
    # ... but override it only if sigma != NA
    if(is.null(misc$sigma) || (length(misc$sigma) == 0) || !is.na(misc$sigma[1]))
        misc$sigma = sigma
    

    post.beta = basis$post.beta
    if (is.null(post.beta))
        post.beta = matrix(NA)
    
    predictors = intersect(attr(data, "predictors"), names(grid))
    # if(!missing(counterfactuals)) predictors = c(predictors, ".obs.no.")
    
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
                      multresp = multresp,
                      nuisance = nuisance),
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
    if(!missing(regrid)) {
        # if(missing(wt.counter)) wt.counter = 1
        result = regrid(result, transform = regrid, sigma = sigma, ...)
    }
    
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

# My replacement for plyr::id(, drop = TRUE)
.my.id = function(data){
    p = do.call(paste, data)
    u = unique(p)
    match(p, u)
}

# Utility to error-out when potential reference grid size is too big
.check.grid = function(levs, limit = get_emm_option("rg.limit")) {
    size = prod(sapply(levs, length))
    if (size > limit)
        stop("The rows of your requested reference grid would be ", size, ", which exceeds\n",
             "the limit of ", limit, " (not including any multivariate responses).\n",
             "Your options are:\n",
             "  1. Specify some (or more) nuisance factors using the 'nuisance' argument\n",
             "     (see ?ref_grid). These must be factors that do not interact with others.\n",
             "  2. Add the argument 'rg.limit = <new limit>' to the call. Be careful,\n",
             "     because this could cause excessive memory use and performance issues.\n",
             "     Or, change the default via 'emm_options(rg.limit = <new limit>)'.\n",
             call. = FALSE)
}

# Utility to set up the grid for nuisance factors. This consists of two or more
# data.frames rbinded together:
#   * the expanded grid for all factors *not* in nuis, with the nuis factors set 
#     at their first levels
#   * for each factor f in nuis, a set of rows for (levs$f), with the other 
#     factors at their first levels (this is arbitrary)
# In addition, we return a character vector 'row.assign' corresponding to the rows 
# in the grid, telling which is what: ".main.grid." for the first part of the grid,
# otherwise factor names from nuis.
# We also return 'nuis' itself - which may be reduced since illegal 
# entries are silently removed
###
# Update July 2024...
# All this trickery assumes that columns of X are associated with model terms!
# So this does not work for anything that's been re-gridded.
# Those (almost?) always have X = I. Unfortunately we won't see X until
# later, but we can still stop it before we inflict damage...
# Note: we also provide a simpler check: misc$regrid.flag is non-NULL
.setup.nuis = function(nuis, levs, trms, rg.limit) {
    firsts =  args = lapply(levs, function(x) x[1])
    nuis = intersect(nuis, names(levs))
    # sanity checks on terms, and term indexes
    fsum = rep(99, length(nuis))
    tbl = attr(trms, "factors")
    rn = row.names(tbl) = sapply(row.names(tbl), function(nm)
        paste(all.vars(reformulate(nm)), collapse = ","))
    for (i in seq_along(nuis)) {
        f = nuis[i]
        if(f %in% rn)
            fsum[i] = sum(tbl[f, ])
    }
    nuis = nuis[fsum == 1]   # silently remove any unfound or interacting factors

    # top part...
    non.nuis = setdiff(names(levs), nuis)
    for (n in non.nuis)
        args[[n]] = levs[[n]]
    .check.grid(args, rg.limit)
    grid = do.call(expand.grid, args)
    ra = rep(".main.grid.", nrow(grid))
    # bottom parts
    for (f in nuis) {
        args = firsts
        args[[f]] = levs[[f]]
        grid = rbind(grid, do.call(expand.grid, args))
        ra = c(ra, rep(f, length(levs[[f]])))
    }
    list(grid = grid, row.assign = ra, nuis = nuis)
}

# Do the required post-processing for nuisance factors...
# After we get the model matrix for this grid, we'll average each set of rows in the
# bottom part, and substitute those averages in the required columns in the top part 
# of the model matrix.
.basis.nuis = function(basis, info, wt, levs, data, grid, ref.levels) {
    X = basis$X
    if(!is.null(basis$misc$regrid.flag) || all(apply(X, 1, \(x) sum(x != 0)) == 1))   # each row has 1 nonzero element
        stop("Sorry, 'nuisance' specs are not allowed for this situation.",
             " Revise the call accordingly.", call. = FALSE)
    ra = info$row.assign
    r. = rep(".", length(ra))  # fillers
    n = sum(ra == ".main.grid.")
    k = nrow(X) / length(ra)   # multivariate dimension
    nuis = info$nuis
    wts = lapply(nuis, function(f) {
        if (wt == "equal")
            w = rep(1, length(levs[[f]]))
        else {
            x = data[[f]]
            w = sapply(levs[[f]], function(lev) sum(x == lev))
        }
        w / sum(w)
    })
    names(wts) = nuis
    
    # In a multivariate case, we have to repeat the same operations for each block of X rows
    for (m in 1:k) {
        RA = c(rep(r., m - 1), ra, rep(r., k - m))
        for (f in nuis) {
            subX = X[RA == f, , drop = FALSE]
            cols = which(apply(subX, 2, function(x) diff(range(x)) > 0))
            subX = sweep(subX[, cols, drop = FALSE], 1, wts[[f]], "*")
            avg = apply(subX, 2, sum)
            avg = matrix(rep(avg, each = n), nrow = n) # now several copies
            X[RA == ".main.grid.", cols] = avg
        }
    }
    basis$misc$nuis = nuis
    basis$misc$avgd.over = paste(length(nuis), "nuisance factors")
    RA = rep(ra, k)
    basis$X = X[RA == ".main.grid.", , drop = FALSE]
    non.nuis = setdiff(names(ref.levels), info$nuis)
    basis$ref.levels = ref.levels[non.nuis]
    basis$grid = grid[ra == ".main.grid.", non.nuis, drop = FALSE]
    basis
}

# Internal function to do reference grid for counterfactuals
.cf.refgrid = function(object, counterfactuals, data, options = list(), ...) {
    if(missing(data))
        data = recover_data(object, ...)
    if(!hasName(data, "(weights)"))
        pwts = rep(1, nrow(data))
    else
        pwts = data[["(weights)"]]
    
    # Start with just the ordinary reference grid
    rg = ref_grid(object, data = data, ...)
    cfac = intersect(counterfactuals, names(rg@levels))
    clevs = rg@levels[cfac]
    cgrid = do.call(expand.grid, clevs)
    # handle nasty fact that character predictors don't act like factors
    for (j in cfac)
        if(is.character(data[[j]])) 
            cgrid[[j]] = as.character(cgrid[[j]])
    
    # Get the stuff we need for each main dataset step
    link = .get.link(rg@misc)
    if(is.null(link))
        link = make.link("identity")
    trms = rg@model.info$terms
    xlev = rg@model.info$xlev
    misc = list()
    offset = .get.offset(object, grid = data)
    k = ifelse(length(mr <- rg@roles$multresp) == 0, 1, length(rg@levels[[mr]])) # grid expansion factor
    
    # Index sets for combinations of factors
    pg = apply(cgrid, 1, paste, collapse = ".")
    pd = apply(data[cfac], 1, paste, collapse = ".")
    cidx = lapply(pg, \(x) which(pd == x))

    # special case for covariates with no matches
    if(all(sapply(cidx, length) == 0))
        cidx = list(seq_len(nrow(data)))
    
    # account for any NAs in bhat
    nonNA = !is.na(rg@bhat)
    # ensure we include all levels of cfacs with data
    all.active = sort(unlist(cidx))
    deadrows = sapply(cidx, function(x) x[1])
    offset = c(offset, rep(mean(offset), length(deadrows)))
    pwts = c(pwts, rep(mean(pwts), length(deadrows)))
    data = rbind(data, data[deadrows, ])
    n = nrow(data)
    mymean = function(x) ifelse(is.null(x), NA, mean(x))
    
    # get means of groups of prior weights
    mpwt = sapply(cidx, \(i) mean(pwts[i]))
    
    ## Compile the averaged results for delta method
    DL = matrix(nrow = 0, ncol = sum(nonNA))
    bh = numeric(0)
    embmeth = .find_method(object, "emm_basis")
    for (i in seq_len(nrow(cgrid))) {
        g = data
        for(j in cfac)
            g[all.active, j] = cgrid[i,j]
        bas = embmeth(object, trms, xlev, g, ...)
        if(!is.null(bas$misc$postGridHook))
            stop("Sorry, we do not support counterfactuals for this situation.")
        X = bas$X[, nonNA, drop = FALSE]
        eta = X %*% bas$bhat[nonNA]
        yhat = link$linkinv(eta + offset)
        d = link$mu.eta(eta) * rep(pwts, k)   # includes prior weights
        X = sweep(X, 1, d, "*")
        
        pos = 0
        XX = matrix(nrow = 0, ncol = ncol(X))
        for(I in 1:k) {
            XX = sapply(cidx, \(i) apply(X[pos + i, , drop = FALSE], 2, mymean)) / mpwt
            DL = rbind(DL, t(XX))
            yy = sapply(cidx, \(i) ifelse(length(i) == 0, NA, mean(yhat[i + pos])))
            bh = c(bh, yy)
            pos = pos + n
        }
    }
    RG = rg
    RG@bhat = bh
    nonNA = !is.na(bh)
    RG@linfct = diag(nrow(DL))
    RG@V = DL %*% rg@V %*% t(DL)
    levs = rg@levels
    levs[cfac] = NULL
    alevs = clevs
    names(alevs) = paste0("actual_", cfac)
    levs = c(alevs, clevs)
    if (k > 1)
        levs = c(levs, rg@levels[length(rg@levels)])
    RG@levels = levs
    wts = sapply(cidx, length) * mpwt
    RG@grid = do.call("expand.grid", levs)
    RG@grid$.wgt. = rep(wts, length(bh)/length(wts))
    misc = rg@misc
    if(!is.null(misc$inv.lbl))
        misc$estName = misc$inv.lbl
    misc[c("tran", "inv.lbl", "sigma")] = NULL
    RG@misc = c(misc, famSize = length(bh), cf.grid = TRUE)
    RG@model.info$model.matrix = NULL
    RG@roles$predictors = c(names(alevs), names(clevs))
    RG@roles$counterfactuals = cfac
    if (all(nonNA))
        RG@nbasis = estimability::all.estble
    else {
        RG@nbasis = matrix(0, ncol = sum(!nonNA), nrow = length(bh))
        idx = which(!nonNA)
        for (j in seq_len(ncol(RG@nbasis)))
            RG@nbasis[idx[j], j] = 1
        RG@V = RG@V[nonNA, nonNA]
    }
    RG
}


## OLD CODE...
# Set up grid for counterfactuals - i.e., copies of whole dataset with
# counterfactual levels substituted, with obs index varying the slowest
# .setup.cf = function(levs, data) {
#     lv = arg = levs[-length(levs)] # omit .obs.no.
#     arg$stringsAsFactors = FALSE
#     g = do.call(expand.grid, arg)
#     idx = rep(seq_len(nrow(data)), each = nrow(g))
#     xdata = data[idx, , drop = FALSE]
#     idx = matrix(seq_len(nrow(xdata)), nrow = nrow(g))
#     for (i in seq_along(g[, 1]))
#         for (j in names(lv))
#             xdata[idx[i, ], j] = g[i, j]
#     xdata
# }
