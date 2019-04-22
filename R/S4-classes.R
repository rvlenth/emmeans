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

### S4 class definitions for emmeans package


### emmGrid object - For reference grids, emmeans results, etc.

#' The \code{emmGrid} class
#' 
#' The \code{emmGrid} class encapsulates linear functions of regression
#' parameters, defined over a grid of predictors. This includes reference
#' grids and grids of marginal means thereof (aka estimated marginal means).
#' Objects of class `emmGrid` may be used independently of the underlying model
#' object. Instances are created primarily by \code{\link{ref_grid}} and
#' \code{\link{emmeans}}, and several related functions.
#'   
#' @rdname emmGrid-class 
#' @slot model.info list. Contains the elements \code{call} (the call that
#'   produced the model), \code{terms} (its \code{terms} object), and
#'   \code{xlev} (factor-level information)
#' @slot roles list. Contains at least the elements \code{predictors}, 
#'   \code{responses}, and \code{multresp}. Each is a character vector of names 
#'   of these variables.
#' @slot grid data.frame. Contains the combinations of the variables that define
#'   the reference grid. In addition, there is an auxiliary column named
#'   \code{".wgt."} holding the observed frequencies or weights for each factor
#'   combination (excluding covariates). If the model has one or more
#'   \code{\link{offset}()} calls, there is an another auxiliary column named
#'   \code{".offset."}. Auxiliary columns are not considered part of the
#'   reference grid. (However, any variables included in \code{offset} calls
#'   \emph{are} in the reference grid.)
#' @slot levels list. Each entry is a character vector with the distinct levels
#'   of each variable in the reference grid. Note that \code{grid} is obtained
#'   by applying the function \code{\link{expand.grid}} to this list
#' @slot matlevs list. Like \code{levels} but has the levels of any matrices in
#'   the original dataset. Matrix columns are always concatenated and treated as
#'   a single variable for purposes of the reference grid
#' @slot linfct matrix. Each row consists of the linear function of the
#'   regression coefficients for predicting its corresponding element of the
#'   reference grid. The rows of this matrix go in one-to-one correspondence
#'   with the rows of \code{grid}, and the columns with elements of \code{bhat}.
#' @slot bhat numeric. The regression coefficients. If there is a multivariate
#'   response, the matrix of coefficients is flattened to a single vector, and
#'   \code{linfct} and \code{V} redefined appropriately. Important: \code{bhat}
#'   must \emph{include} any \code{NA} values produced as a result of 
#'   collinearity in the predictors. These are taken care of later in the 
#'   estimability check.
#' @slot nbasis matrix. The basis for the non-estimable functions of the
#'   regression coefficients. Every EMM will correspond to a linear combination
#'   of rows of \code{linfct}, and that result must be orthogonal to all the
#'   columns of \code{nbasis} in order to be estimable. If everything is
#'   estimable, \code{nbasis} should be a 1 x 1 matrix of \code{NA}.
#' @slot V matrix. The symmetric variance-covariance matrix of \code{bhat}
#' @slot dffun function having two arguments. \code{dffun(k, dfargs)} should
#'   return the degrees of freedom for the linear function \code{sum(k*bhat)},
#'   or \code{NA} if unavailable
#' @slot dfargs list. Used to hold any additional information needed by
#'   \code{dffun}.
#' @slot misc list. Additional information used by methods. These include at
#'   least the following: \code{estName} (the label for the estimates of linear
#'   functions), and the default values of \code{infer}, \code{level}, and
#'   \code{adjust} to be used in the \code{\link{summary.emmGrid}} method. Elements in
#'   this slot may be modified if desired using the \code{\link{update.emmGrid}}
#'   method.
#' @slot post.beta matrix. A sample from the posterior distribution of the
#'   regression coefficients, if MCMC methods were used; or a 1 x 1 matrix of
#'   \code{NA} otherwise. When it is non-trivial, the \code{\link{as.mcmc.emmGrid}}
#'   method returns \code{post.beta \%*\% t(linfct)}, which is a sample from the
#'   posterior distribution of the EMMs.
#' 
#' @section Methods:
#'   All methods for these objects are S3 methods except for \code{show}. 
#'   They include \code{\link{[.emmGrid}}, \code{\link{as.glht.emmGrid}},
#'   \code{\link{as.mcmc.emmGrid}}, \code{\link{as.mcmc.list.emmGrid}} (see \pkg{coda}),
#'   \code{\link{cld.emmGrid}} (see \pkg{multcomp}),
#'   \code{\link{coef.emmGrid}}, \code{\link{confint.emmGrid}}, 
#'   \code{\link{contrast.emmGrid}}, \code{\link{pairs.emmGrid}},
#'   \code{\link{plot.emmGrid}}, \code{\link{predict.emmGrid}}, \code{\link{print.emmGrid}},
#'   \code{\link{rbind.emmGrid}}, \code{show.emmGrid}, \code{\link{str.emmGrid}}, 
#'   \code{\link{summary.emmGrid}}, \code{\link{test.emmGrid}}, 
#'   \code{\link{update.emmGrid}}, \code{\link{vcov.emmGrid}}, and 
#'   \code{\link{xtable.emmGrid}}
#' 
#' @export
setClass("emmGrid", slots = c(
    model.info = "list",
    roles = "list",
    grid = "data.frame", 
    levels = "list",
    matlevs = "list",
    linfct = "matrix",
    bhat = "numeric",
    nbasis = "matrix",
    V = "matrix",
    dffun = "function",
    dfargs = "list",
    misc = "list",
    post.beta = "matrix"
))
# Note: misc will hold various extra params,
# including at least the following req'd by the summary method
#   estName: column name for the estimate in the summary ["prediction"]
#   infer: booleans (CIs?, tests?)  [(FALSE,FALSE)]
#   level: default conf level [.95]
#   adjust: default adjust method ["none"]
#   famSize: number of means in family

# NOTE: Old ref.grid and lsmobj classes moved to deprecated.R