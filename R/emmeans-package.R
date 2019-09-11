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

#' Estimated marginal means (aka Least-squares means)
#' 
#' This package provides methods for obtaining estimated marginal means (EMMs, also 
#' known as least-squares means) for factor combinations in a variety of models.
#' Supported models include [generalized linear] models, models for counts,
#' multivariate, multinomial and ordinal responses, survival models, GEEs, and
#' Bayesian models. For the latter, posterior samples of EMMs are provided.
#' The package can compute contrasts or linear
#' combinations of these marginal means with various multiplicity adjustments.
#' One can also estimate and contrast slopes of trend lines.
#' Some graphical displays of these results are provided.
#' 
#' @section Overview:
#' \describe{
#' \item{Vignettes}{A number of vignettes are provided to help the user get
#' acquainted with the \pkg{emmeans} package and see some examples.}
#' 
#' \item{Concept}{Estimated marginal means (see Searle \emph{et al.} 1980 are
#' popular for summarizing linear models that include factors. For balanced
#' experimental designs, they are just the marginal means. For unbalanced data,
#' they in essence estimate the marginal means you \emph{would} have observed
#' that the data arisen from a balanced experiment. Earlier developments
#' regarding these techniques were developed in a least-squares context and are
#' sometimes referred to as \dQuote{least-squares means}. Since its early
#' development, the concept has expanded far beyond least-squares settings.}
#' 
#' \item{Reference grids}{ The implementation in \pkg{emmeans} relies on our own
#' concept of a \emph{reference grid}, which is an array of factor and predictor
#' levels. Predictions are made on this grid, and estimated marginal means (or
#' EMMs) are defined as averages of these predictions over zero or more
#' dimensions of the grid. The function \code{\link{ref_grid}} explicitly
#' creates a reference grid that can subsequently be used to obtain
#' least-squares means. The object returned by \code{ref_grid} is of class
#' \code{"emmGrid"}, the same class as is used for estimated marginal means (see
#' below).
#' 
#' Our reference-grid framework expands slightly upon Searle \emph{et al.}'s
#' definitions of EMMs, in that it is possible to include multiple levels of
#' covariates in the grid. }
#' 
#' \item{Models supported}{As is mentioned in the package description, many
#' types of models are supported by the package. 
#' See \href{../doc/models.html}{vignette("models", "emmeans")} for full details. 
#' Some models may require other packages be
#' installed in order to  access all of the available features.}
#' 
#' \item{Estimated marginal means}{
#' The \code{\link{emmeans}} function computes EMMs given a fitted model (or a
#' previously constructed \code{emmGrid} object), using a specification indicating
#' what factors to include. The \code{\link{emtrends}} function creates the same
#' sort of results for estimating and comparing slopes of fitted lines. Both
#' return an \code{emmGrid} object.}
#' 
#' \item{Summaries and analysis}{
#' The \code{\link{summary.emmGrid}}  method may be used to display an \code{emmGrid}
#' object. Special-purpose summaries are available via \code{\link{confint.emmGrid}} and
#' \code{\link{test.emmGrid}}, the latter of which can also do a joint test of several
#' estimates. The user may specify by variables, multiplicity-adjustment
#' methods, confidence levels, etc., and if a transformation or link function is
#' involved, may reverse-transform the results to the response scale.}
#' 
#' \item{Contrasts and comparisons}{
#' The \code{\link{contrast}} method for \code{emmGrid} objects is used to obtain
#' contrasts among the estimates; several standard contrast families are
#' available such as deviations from the mean, polynomial contrasts, and
#' comparisons with one or more controls. Another \code{emmGrid} object is returned,
#' which can be summarized or further analyzed. For convenience, a \code{pairs.emmGrid}
#' method is provided for the case of pairwise comparisons. 
#' }
#' \item{Graphs}{The \code{\link{plot.emmGrid}} method will display
#' side-by-side confidence intervals for the estimates, and/or
#' \dQuote{comparison arrows} whereby the *P* values of pairwise differences
#' can be observed by how much the arrows overlap. The \code{\link{emmip}} function
#' displays estimates like an interaction plot, multi-paneled if there are by
#' variables. These graphics capabilities require the \pkg{lattice} package be
#' installed.}
#' 
#' \item{MCMC support}{When a model is fitted using MCMC methods, the posterior
#' chains(s) of parameter estimates are retained and converted into posterior
#' samples of EMMs or contrasts thereof. These may then be summarized or plotted
#' like any other MCMC results, using tools in, say \pkg{coda} or
#' \pkg{bayesplot}.}
#' 
#' \item{\pkg{multcomp} interface}{The \code{\link{as.glht}} function and
#' \code{glht} method for \code{emmGrid}s provide an interface to the
#' \code{glht} function in the \pkg{multcomp} package, thus
#' providing for more exacting simultaneous estimation or testing. The package
#' also provides an \code{\link{emm}} function that works as an alternative to
#' \code{mcp} in a call to \code{glht}.
#' }
#' } %%% end describe
#' 
#' @import estimability
#' @import mvtnorm
#' @import stats
#' @importFrom graphics pairs plot
#' @importFrom methods as is new slot slot<- slotNames 
#' @importFrom utils getS3method installed.packages methods str
#' @name emmeans-package
NULL
