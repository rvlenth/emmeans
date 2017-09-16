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

# documentation for models

#' Models supported in \pkg{emmeans}
#'
#' Here we document what model objects may be used with \pkg{emmeans}, and some
#' special features of some of them that may be accessed by passing additional
#' arguments through \code{\link{ref_grid}} or \code{\link{emmeans}}.
#'
#' Certain objects are affected by optional arguments to functions that
#' construct \code{emm} objects, including \code{\link{ref_grid}},
#' \code{\link{emmeans}}, \code{\link{emtrends}}, and \code{\link{emmip}}. When
#' \dQuote{arguments} are mentioned in the subsequent quick reference and
#' object-by-object documentation, we are talking about arguments in these
#' constructors.
#'
#' Additional models can be supported by writing appropriate \code{recover_data}
#' and \code{emm_basis} methods. See \code{\link{extending-emmeans}} and
#' \code{vignette("extending")} for details.
#'
#' @section Quick reference for supported objects and options:
#' Here is an alphabetical list of model classes that are supported, and
#' the arguments that apply. Detailed documentation follows, with objects
#' grouped by the code in the \dQuote{Group} column.
#' \tabular{llcl}{
#'   \bold{Object class} \tab \bold{Package} \tab \bold{Group} \tab \bold{Arguments / notes} \cr
#'   \code{aov} \tab \pkg{stats} \tab A \tab  \cr
#'   \code{aovList} \tab \pkg{stats} \tab V \tab Best with balanced designs, orthogonal coding \cr
#'   \code{betareg} \tab \pkg{betareg} \tab B \tab \code{mode = c("link", "precision", "phi.link",)} \cr
#'            \tab \tab \tab \code{"variance", "quantile")} \cr
#'   \code{carbayes} \tab \pkg{CARBayes} \tab S \tab \code{data} is \emph{required} \cr
#'   \code{clm} \tab \pkg{ordinal} \tab O \tab \code{mode = c("latent", "linear.predictor", "cum.prob",} \cr
#'               \tab \tab \tab \code{"exc.prob", "prob", "mean.class", "scale")} \cr
#'   \code{clmm} \tab \pkg{ordinal} \tab O \tab Like \code{clm} but no \code{"scale"} mode \cr
#'   \code{coxme} \tab \pkg{coxme} \tab G \tab \cr
#'   \code{coxph} \tab \pkg{survival} \tab G \cr
#'   \code{gee} \tab \pkg{gee} \tab E \tab \code{vcov.method = c("naive", "robust")} \cr
#'   \code{geeglm} \tab \pkg{geepack} \tab E \tab \code{vcov.method = c("vbeta", "vbeta.naiv", "vbeta.j1s",} \cr
#'                      \tab  \tab \tab \code{"vbeta.fij", "robust", "naive")} or a matrix \cr
#'   \code{geese} \tab \pkg{geepack} \tab E \tab Like \code{geeglm} \cr
#'   \code{glm} \tab \pkg{stats} \tab G \tab  \cr
#'   \code{glm.nb} \tab \pkg{MASS} \tab G \tab Requires \code{data} argument \cr
#'   \code{glmmadmb} \tab \pkg{glmmADMB} \tab G \cr
#'   \code{glmerMod} \tab \pkg{lme4} \tab G \cr
#'   \code{glmmPQL} \tab \pkg{MASS} \tab G \tab inherits \code{lm} support \cr
#'   \code{gls} \tab \pkg{nlme} \tab G \cr
#'   \code{hurdle} \tab \pkg{pscl} \tab C \tab \code{mode = c("response", "count", "zero", "prob0")},\cr
#'               \tab  \tab \tab \code{lin.pred = c(FALSE, TRUE)} \cr
#'   \code{lm} \tab \pkg{stats} \tab A \tab  Several other classes inherit from this and may be supported \cr
#'   \code{lme} \tab \pkg{nlme} \tab A \tab \code{sigmaAdjust = c(TRUE, FALSE)} \cr
#'   \code{lmerMod} \tab \pkg{lme4} \tab L \tab \code{mode = c("kenward-roger", "satterthwaite",} \cr
#'               \tab \tab \tab \code{ "asymptotic")}. Also \code{emm_options(disable.pbkrtest =,} \cr
#'                 \tab  \tab \tab \code{lmer.df =, pbkrtest.limit =)} \cr
#'   \code{manova} \tab \pkg{stats} \tab M \tab  \code{mult.name}, \code{mult.levs} \cr
#'   \code{maov} \tab \pkg{stats} \tab M \tab  \code{mult.name}, \code{mult.levs} \cr
#'   \code{mcmc} \tab \pkg{mcmc} \tab S \tab May require \code{formula}, \code{data} \cr
#'   \code{MCMCglmm} \tab \pkg{MCMCglmm} \tab M,S \tab \code{mult.name}, \code{mult.levs}\cr
#'                \tab \tab \tab \code{data} is required \cr
#'   \code{mer} \tab \pkg{lme4.0} \tab A \tab May become obsolete \cr
#'   \code{mixed} \tab \pkg{afex} \tab P \tab Supported in \pkg{afex} package \cr
#'   \code{mlm} \tab \pkg{stats} \tab M \tab  \code{mult.name}, \code{mult.levs} \cr
#'   \code{multinom} \tab \pkg{nnet} \tab N \tab \code{mode = c("prob", "latent")} \cr
#'         \tab  \tab \tab Always include response in \code{specs} for \code{emmeans} \cr
#'   \code{nauf} \tab \pkg{} \tab P \tab Supported in \pkg{nauf} package \cr
#'   \code{nlme} \tab \pkg{nlme} \tab A \tab Supports \code{fixed} part. Requires \code{param} \cr
#'   \code{polr} \tab \pkg{MASS} \tab O \tab \code{mode = c("latent", "linear.predictor", "cum.prob",} \cr
#'              \tab  \tab \tab \code{"exc.prob", "prob", "mean.class")} \cr
#'   \code{rlm} \tab \pkg{MASS} \tab A \tab inherits \code{lm} support \cr
#'   \code{rms} \tab \pkg{rms} \tab O \tab \code{mode = "middle", "latent", "linear.predictor",}\cr
#'              \tab \tab \tab \code{"cum.prob", "exc.prob", "prob", "mean.class")} \cr
#'   \code{rsm} \tab \pkg{rsm} \tab P \tab Supported in \pkg{rsm} package \cr
#'   \code{stanreg} \tab \pkg{rstanarm} \tab S \tab Args for \code{stanreg_xxx} similar to those for \code{xxx}\cr
#'   \code{survreg} \tab \pkg{survival} \tab A \cr
#'   \code{zeroinf} \tab \pkg{pscl} \tab C \tab \code{mode = c("response", "count", "zero", "prob0")},\cr
#'                 \tab \tab \tab \code{lin.pred = c(FALSE, TRUE)} \cr
#' } % end of tabular
#' 
#' @section Group A -- \dQuote{Standard} models (typically linear and mixed):
#' Models in this group, such as \code{lm}, do not have unusual features that
#' need special support; hence no extra arguments are needed.
#' 
#' @section B -- Beta regression:
#' The additional \code{mode} argument for \code{betareg} objects has possible
#' values of \code{"response"}, \code{"link"}, \code{"precision"},
#' \code{"phi.link"}, \code{"variance"}, and \code{"quantile"}, which have the
#' same meaning as the \code{type} argument in \code{predict.betareg} -- with
#' the addition that \code{"phi.link"} is like \code{"link"}, but for the
#' precision portion of the model. When \code{mode = "quantile"} is specified,
#' the additional argument \code{quantile} (a numeric scalar or vector)
#' specifies which quantile(s) to compute; the default is 0.5 (the median). Also
#' in \code{"quantile"} mode, an additional variable \code{quantile} is added to
#' the reference grid, and its levels are the values supplied.
#' 
#' @section Group C -- Count models:
#' Two optional arguments -- \code{mode} and \code{lin.pred} -- are provided.
#' The \code{mode} argument has possible values \code{"response"} (the default),
#' \code{"count"}, \code{"zero"}, or \code{"prob0"}. \code{lin.pred} is logical
#' and defaults to \code{FALSE}.
#' 
#' With \code{lin.pred = FALSE}, the results are comparable to those returned by
#' \code{predict(..., type = "response")}, \code{predict(..., type = "count")},
#' \code{predict(..., type = "zero")}, or \code{predict(..., type = "prob")[,
#' 1]}. See the documentation for \code{\link[pscl]{predict.hurdle}} and
#' \code{\link[pscl]{predict.zeroinfl}}.
#' 
#' The option \code{lin.pred = TRUE} only applies to \code{mode = "count"} and
#' \code{mode = "zero"}. The results returned are on the linear-predictor scale,
#' with the same transformation as the link function in that part of the model.
#' The predictions for a reference grid with \code{mode = "count"},
#' \code{lin.pred = TRUE}, and \code{type = "response"} will be the same as
#' those obtained with \code{lin.pred = FALSE} and \code{mode = "count"};
#' however, any LS means derived from these grids will be different, because the
#' averaging is done on the log-count scale and the actual count scale,
#' respectively -- thereby producing geometric means versus arithmetic means of
#' the predictions.
#' 
#' If the \code{vcov.} argument is used (see details in \code{\link{ref.grid}}),
#' it must yield a matrix of the same size as would be obtained using
#' \code{\link[pscl]{vcov.hurdle}} or \code{\link[pscl]{vcov.zeroinfl}} with its
#' \code{model} argument set to \code{("full", "count", "zero")} in respective
#' correspondence with \code{mode} of \code{("mean", "count", "zero")}. If
#' \code{vcov.} is a function, it must support the \code{model} argument.
#' 
#' @section Group E -- GEE models:
#' These models all have more than one covariance estimate available, and it may
#' be selected by supplying a string as the \code{vcov.method} argument. It is
#' partially matched with the available choices shown in the quick reference. In
#' \code{geese} and \code{geeglm}, the aliases \code{"robust"} (for
#' \code{"vbeta"}) and \code{"naive"} (for \code{"vbeta.naiv"} are also
#' accepted.
#' 
#' If a matrix or function is supplied as \code{vcov.method}, it is
#' interpreted as a \code{vcov.} specification as described for \code{...}
#' in \code{\link{ref_grid}}.
#'  
#' @section Group G -- Generalized linear models:
#' Models in this group receive only standard support as in Group A, but
#' typically the tests and confidence intervals are asymptotic. Thus the
#' \code{df} column for tabular results will be \code{NA}.
#' 
#' @section Group L -- \code{lmerMod}:
#' There is an optional \code{mode} argument that defaults to 
#' \code{get_EMM_option("lmer.df")} (which in turn defaults to 
#' \code{"kenward-roger"}). The possible values are "kenward-roger",
#' "satterthwaite", and "asymptotic" (these are partially matched and 
#' case-insensitive). With \code{"kenward-roger"}, d.f. are obtained using code
#' from the \pkg{pbkrtest} package, if installed. With \code{"satterthwaite"},
#' d.f. are obtained using code from the \pkg{lmerTest} package, if installed. 
#' With \code{"asymptotic"}, or if the needed package is not installed, d.f. are
#' set to \code{NA}.
#' 
#' A by-product of the Kenward-Roger method is that the covariance matrix is 
#' adjusted using \code{\link[pbkrtest]{vcovAdj}}. This can require considerable
#' computation; so to avoid that overhead, the user should opt for the 
#' Satterthwaite or asymptotic method; or, for backward compatibility, may 
#' disable the use of \pkg{pbkrtest} via \samp{emm_options(disable.pbkrtest =
#' TRUE)} (this does not disable the \pkg{pbkrtest} package entirely, just its
#' use in \pkg{emmeans}). The computation time required depends roughly on the
#' number of observations, \emph{N}, in the design matrix (because a major part
#' of the computation involves inverting an \emph{N x N} matrix). Thus,
#' \pkg{pbkrtest} is automatically disabled if \emph{N} exceeds the value of 
#' \code{get_emm_option("pbkrtest.limit")}. If desired, the user may use 
#' \code{emm_options} to adjust this limit from the default of 3000.
#' 
#' The \code{df} argument may be used to specify some other degrees of freedom.
#' Note that if \code{df} and \code{method = "satterthwaite"} are both
#' specified, the covariance matrix is adjusted but the K-R degrees of freedom
#' are not used.
#' 
#' @section Group M -- Multivariate models:
#' When there is a multivariate response, the different responses are treated as
#' if they were levels of a factor -- named \code{rep.meas} by default. The
#' \code{mult.name} argument may be used to change this name. The
#' \code{mult.levs} argument may specify a named list of one or more sets of
#' levels. If this has more than one element, then the multivariate levels are
#' expressed as combinations of the named factor levels via the function
#' \code{\link{expand.grid}}.
#' 
#' @section N - Multinomial responses:
#' The reference grid includes a pseudo-factor with the same name and levels as
#' the multinomial response. There is an optional \code{mode} argument which
#' should match \code{"prob"} or \code{"latent"}. With \code{mode = "prob"}, the
#' reference-grid predictions consist of the estimated multinomial
#' probabilities. The \code{"latent"} mode returns the linear predictor,
#' recentered so that it averages to zero over the levels of the response
#' variable (similar to sum-to-zero contrasts). Thus each latent variable can be
#' regarded as the log probability at that level minus the average log
#' probability over all levels.
#' 
#' There are two optional arguments: \code{mode} and \code{rescale} (which
#' defaults to \code{c(0, 1)}). \code{polr} does not support scale models.
#' 
#' Please note that, because the probabilities sum to 1 (and the latent values
#' sum to 0) over the multivariate-response levels, all sensible results from
#' \code{emmeans} must involve that response as one of the factors. For example,
#' if \code{resp} is a response with \eqn{k} levels, \code{emmeans(model, ~ resp
#' | trt)} will yield the estimated multinomial distribution for each
#' \code{trt}; but \code{emmeans(model, ~ trt)} will just yield the average
#' probability of \eqn{1/k} for each \code{trt}.
#' 
#' 
#' @section Group O - Ordinal responses:
#' The reference grid for ordinal models will include all variables that appear in
#' the main model
#' as well as those in the \code{scale} or \code{nominal} models (if provided). 
#' There are two optional arguments: \code{mode} (a character string) and 
#' \code{rescale} (which defaults to \samp{c(0, 1)}). \code{mode} should match 
#' one of \code{"latent"} (the default), \code{"linear.predictor"}, 
#' \code{"cum.prob"}, \code{"exc.prob"}, \code{"prob"}, \code{"mean.class"}, or
#' \code{"scale"} -- see the quick reference and note which are supported.
#' 
#' With \samp{mode = "latent"}, the reference-grid predictions are made on the
#' scale of the latent variable implied by the model. The scale and location of
#' this latent variable are arbitrary, and may be altered via \code{rescale}.
#' The predictions are multiplied by \samp{rescale[2]}, then \samp{rescale[1]}
#' is added. Keep in mind that the scaling is related to the link function used
#' in the model; for example, changing from a probit link to a logistic link
#' will inflate the latent values by around \eqn{\pi/\sqrt{3}}{pi/sqrt(3)}, all
#' other things being equal. \code{rescale} has no effect for other values of
#' \code{mode}.
#' 
#' With \samp{mode = "linear.predictor"} \code{mode = "cum.prob"}, and
#' \code{mode = "exc.prob"}, the boundaries between categories (i.e.,
#' thresholds) in the ordinal response are included in  the reference grid as a
#' pseudo-factor named \code{cut}. The reference-grid predictions are then of
#' the cumulative probabilities at each threshold (for \code{mode =
#' "cum.prob"}), exceedance probabilities (one minus cumulative probabilities,
#' for \code{mode = "exc.prob"}), or the link function thereof (for \code{mode =
#' "linear.predictor"}).
#' 
#' With \code{mode = "prob"}, a pseudo-factor with the same name as the model's
#' response variable is created, and the grid predictions are of the
#' probabilities of each class of the ordinal response. With
#' \code{"mean.class"}, the returned results are means of the ordinal response,
#' interpreted as a numeric value from 1 to the number of classes, using the
#' \code{"prob"} results as the estimated probability distribution for each
#' case.
#' 
#' With \code{mode = "scale"}, and the fitted object incorporates a scale model,
#' EMMs are obtained for the factors in the scale model (with a log response)
#' instead of the response model. The grid is constructed using only the factors
#' in the scale model.
#' 
#' Any grid point that is non-estimable by either the location or the scale
#' model (if present) is set to \code{NA}, and any EMMs involving such a
#' grid point will also be non-estimable. A consequence of this is that if there
#' is a rank-deficient \code{scale} model, and then \emph{all} latent responses
#' become non-estimable because the predictions are made using the average
#' log-scale estimate.
#' 
#' \code{rms} models have an additional \code{mode}. With \code{mode = "middle"}
#' (this is the default), the middle intercept is used, comparable to the
#' default for \pkg{rms}'s \code{Predict} function. This is quite similar in
#' concept to \code{mode = "latent"}, where all intercepts are averaged
#' together.
#' 
#' 
#' @section P -- Other packages:
#' Models in this group have their \pkg{emmeans} support provided by the package
#' that implements the model-fitting procedure. Users should refer to the
#' package documentation for details on \pkg{emmeans} support.
#' 
#' @section S -- Sampling methods:
#' Models fitted using MCMC methods contain a sample from the posterior
#' distribution of fixed-effect coefficients. In some cases (e.g., results of
#' \code{MCMCregress} and \code{MCMCpoisson}), the object may include a
#' \code{"call"} attribute that \code{emmeans} can use to reconstruct the data
#' and obtain a basis for the least-squares means. If not, a \code{formula} and
#' \code{data} argument are provided that may help produce the right results. In
#' addition, the \code{contrasts} specifications are not necessarily recoverable
#' from the object, so the system default must match what was actually used in
#' fitting the model.
#' 
#' The usual \code{summary}, \code{test},
#' etc. methods provide frequentist analyses of the results based on the
#' posterior means and covariances. However, an \code{as.mcmc} method is
#' provided that creates an \code{mcmc} object that can be summarized or plotted
#' using the \pkg{coda} package. It provides a posterior sample of EMMs for the
#' given reference grid, based on the posterior sample of the fixed effects from
#' themodel object.
#' 
#' @section Group V -- \code{aovList}:
#' Support for these objects is limited. To avoid strong biases in the
#' predictions, the \code{contrasts} attribute of all factors should be of a
#' type that sums to zero -- for example, \code{"contr.sum"},
#' \code{"contr.poly"}, or \code{"contr.helmert"} but \emph{not}
#' \code{"contr.treatment"}.  Only intra-block estimates of covariances are
#' used. That is, if a factor appears in more than one error stratum, only the
#' covariance structure from its lowest stratum is used in estimating standard
#' errors. Degrees of freedom are obtained using the Satterthwaite method. In
#' general, \code{aovList} support is best with balanced designs, and due
#' caution in the use of contrasts. If a \code{vcov.} argument is supplied, it
#' must yield a single covariance matrix for the unique fixed effects, and the
#' degrees of freedom are set to \code{NA}.
#'
#' @name models
NULL
