##############################################################################
#    Copyright (c) 2012-2018 Russell V. Lenth                                #
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

### This file has summary.emmGrid S3 method and related functions


# S3 summary method

#' Summaries, predictions, intervals, and tests for \code{emmGrid} objects
#' 
#' These are the primary methods for obtaining numerical or tabular results 
#' from an \code{emmGrid} object. 
#' \code{summary.emmGrid} is the general function for summarizing \code{emmGrid} objects. 
#' It also serves as the print method for these objects; so for convenience,
#' \code{summary()} arguments may be included in calls to functions such as 
#' \code{\link{emmeans}} and \code{\link{contrast}} that construct \code{emmGrid} 
#' objects. Note that by default, summaries for Bayesian models are
#' diverted to \code{\link{hpd.summary}}.
#' 
#' \code{confint.emmGrid} is equivalent to \code{summary.emmGrid with 
#' infer = c(TRUE, FALSE)}. When called with \code{joint = FALSE}, \code{test.emmGrid}
#' is equivalent to \code{summary.emmGrid} with \code{infer = c(FALSE, TRUE)}. 
#' 
#' With \code{joint = TRUE}, \code{test.emmGrid} calculates the Wald test of the
#' hypothesis \code{linfct \%*\% bhat = null}, where \code{linfct} and
#' \code{bhat} refer to slots in \code{object} (possibly subsetted according to
#' \code{by} or \code{rows}). An error is thrown if any row of \code{linfct} is
#' non-estimable. It is permissible for the rows of \code{linfct} to be linearly
#' dependent, as long as \code{null == 0}, in which case a reduced set of 
#' contrasts is tested. Linear dependence and nonzero \code{null} cause an 
#' error.
#'
#' @param object An object of class \code{"emmGrid"} (see \link{emmGrid-class})
#' @param infer A vector of one or two logical values. The first determines
#'   whether confidence intervals are displayed, and the second determines
#'   whether \emph{t} tests and \emph{P} values are displayed. If only one value
#'   is provided, it is used for both.
#' @param level Numerical value between 0 and 1. Confidence level for confidence
#'   intervals, if \code{infer[1]} is \code{TRUE}.
#' @param adjust Character value naming the method used to adjust \eqn{p} values
#'   or confidence limits; or to adjust comparison arrows in \code{plot}. See
#'   the P-value adjustments section below.
#' @param by Character name(s) of variables to use for grouping into separate 
#'   tables. This affects the family of tests considered in adjusted \emph{P}
#'   values. 
#' @param type Character: type of prediction desired. This only has an effect if
#'   there is a known transformation or link function. \code{"response"} 
#'   specifies that the inverse transformation be applied. \code{"mu"} (or 
#'   equivalently, \code{"unlink"}) is usually the same as \code{"response"},
#'   but in the case where the model has both a link function and a response 
#'   transformation, only the link part is back-transformed. Other valid values 
#'   are \code{"link"}, \code{"lp"}, and \code{"linear.predictor"}; these are
#'   equivalent, and request that results be shown for the linear predictor,
#'   with no back-transformation. The default is \code{"link"}, unless the 
#'   \code{"predict.type"} option is in force; see \code{\link{emm_options}},
#'   and also the section below on transformations and links.
#' @param df Numeric. If non-missing, a constant number of degrees of freedom to
#'   use in constructing confidence intervals and \emph{P} values (\code{NA}
#'   specifies asymptotic results).
#' @param calc Named list of character value(s) or formula(s).
#'   The expressions in \code{char} are evaluated and appended to the
#'   summary, just after the \code{df} column. The expression may include
#'   any names up through \code{df} in the summary, any additional names in 
#'   \code{object@grid} (such as \code{.wgt.} or \code{.offset.}), or any
#'   earlier elements of \code{calc}.
#' @param null Numeric. Null hypothesis value(s), on the linear-predictor scale,
#'   against which estimates are tested. May be a single value used for all, or
#'   a numeric vector of length equal to the number of tests in each family
#'   (i.e., \code{by} group in the displayed table).
#' @param delta Numeric value (on the linear-predictor scale). If zero, ordinary
#'   tests of significance are performed. If positive, this specifies a
#'   threshold for testing equivalence (using the TOST or two-one-sided-test
#'   method), non-inferiority, or non-superiority, depending on \code{side}. See
#'   Details for how the test statistics are defined.
#' @param side Numeric or character value specifying whether the test is
#'   left-tailed (\code{-1}, \code{"-"}, code{"<"}, \code{"left"}, or
#'   \code{"nonsuperiority"}); right-tailed (\code{1}, \code{"+"}, \code{">"},
#'   \code{"right"}, or \code{"noninferiority"}); or two-sided (\code{0},
#'   \code{2}, \code{"!="}, \code{"two-sided"}, \code{"both"},
#'   \code{"equivalence"}, or \code{"="}). See the special section below for
#'   more details.
#' @param frequentist Ignored except if a Bayesian model was fitted. If missing
#'   or \code{FALSE}, the object is passed to \code{\link{hpd.summary}}. Otherwise, 
#'   a logical value of \code{TRUE} will have it return a frequentist summary.
#' @param bias.adjust Logical value for whether to adjust for bias in
#'   back-transforming (\code{type = "response"}). This requires a value of 
#'   \code{sigma} to exist in the object or be specified.
#' @param sigma Error SD assumed for bias correction (when 
#'   \code{type = "response"} and a transformation
#'   is in effect), or for constructing prediction intervals. If not specified,
#'   \code{object@misc$sigma} is used, and an error is thrown if it is not found.
#'   \emph{Note:} \code{sigma} may be a vector, as long as it conforms to the number of rows
#'   of the reference grid.
#' @param ... Optional arguments such as \code{scheffe.rank} 
#'   (see \dQuote{P-value adjustments}). 
#'   In \code{as.data.frame.emmGrid}, \code{confint.emmGrid}, 
#'   \code{predict.emmGrid}, and 
#'   \code{test.emmGrid}, these arguments are passed to
#'   \code{summary.emmGrid}.
#'
#' @return \code{summary.emmGrid}, \code{confint.emmGrid}, and
#'   \code{test.emmGrid} return an object of class \code{"summary_emm"}, which
#'   is an extension of \code{\link{data.frame}} but with a special \code{print}
#'   method that displays it with custom formatting. For models fitted using
#'   MCMC methods, the call is diverted to \code{\link{hpd.summary}} (with 
#'   \code{prob} set to \code{level}, if specified); one may
#'   alternatively use general MCMC summarization tools with the 
#'   results of \code{as.mcmc}.
#'   
#' @section Defaults:
#'   The \code{misc} slot in \code{object} may contain default values for
#'   \code{by}, \code{calc}, \code{infer}, \code{level}, \code{adjust}, 
#'   \code{type}, \code{null}, \code{side}, and \code{delta}. 
#'   These defaults vary depending
#'   on the code that created the object. The \code{\link{update}} method may be
#'   used to change these defaults. In addition, any options set using 
#'   \samp{emm_options(summary = ...)} will trump those stored in the object's 
#'   \code{misc} slot.
#' 
#' @section Transformations and links:
#'   With \code{type = "response"}, the transformation assumed can be found in
#'   \samp{object@misc$tran}, and its label, for the summary is in
#'   \samp{object@misc$inv.lbl}. Any \eqn{t} or \eqn{z} tests are still performed
#'   on the scale of the linear predictor, not the inverse-transformed one.
#'   Similarly, confidence intervals are computed on the linear-predictor scale,
#'   then inverse-transformed.
#'   
#'   When \code{bias.adjust} is \code{TRUE}, then back-transformed estimates
#'   are adjusted by adding 
#'   \eqn{0.5 h''(u)\sigma^2}, where \eqn{h} is the inverse transformation and
#'   \eqn{u} is the linear predictor. This is based on a second-order Taylor
#'   expansion. There are better or exact adjustments for certain specific
#'   cases, and these may be incorporated in future updates.
#' 
#' @section P-value adjustments:
#'   The \code{adjust} argument specifies a multiplicity adjustment for tests or
#'   confidence intervals. This adjustment always is applied \emph{separately}
#'   to each table or sub-table that you see in the printed output (see
#'   \code{\link{rbind.emmGrid}} for how to combine tables). 
#'   
#'   The valid values of \code{adjust} are as follows:
#'   \describe{
#'   \item{\code{"tukey"}}{Uses the Studentized range distribution with the number
#'     of means in the family. (Available for two-sided cases only.)}
#'   \item{\code{"scheffe"}}{Computes \eqn{p} values from the \eqn{F}
#'     distribution, according to the Scheffe critical value of
#'     \eqn{\sqrt{rF(\alpha; r, d)}}{sqrt[r*qf(alpha, r, d)]}, where \eqn{d} is
#'     the error degrees of freedom and \eqn{r} is the rank of the set of linear
#'     functions under consideration. By default, the value of \code{r} is
#'     computed from \code{object@linfct} for each by group; however, if the
#'     user specifies an argument matching \code{scheffe.rank}, its value will
#'     be used instead. Ordinarily, if there are \eqn{k} means involved, then
#'     \eqn{r = k - 1} for a full set of contrasts involving all \eqn{k} means, and
#'     \eqn{r = k} for the means themselves. (The Scheffe adjustment is available
#'     for two-sided cases only.)}
#'   \item{\code{"sidak"}}{Makes adjustments as if the estimates were independent
#'     (a conservative adjustment in many cases).}
#'   \item{\code{"bonferroni"}}{Multiplies \eqn{p} values, or divides significance
#'     levels by the number of estimates. This is a conservative adjustment.}
#'   \item{\code{"dunnettx"}}{Uses our own\emph{ad hoc} approximation to the 
#'     Dunnett distribution for a family of estimates having pairwise
#'     correlations of \eqn{0.5} (as is true when comparing treatments with a
#'     control with equal sample sizes). The accuracy of the approximation
#'     improves with the number of simultaneous estimates, and is much faster
#'     than \code{"mvt"}. (Available for two-sided cases only.)}
#'   \item{\code{"mvt"}}{Uses the multivariate \eqn{t} distribution to assess the
#'     probability or critical value for the maximum of \eqn{k} estimates. This
#'     method produces the same \eqn{p} values and intervals as the default
#'     \code{summary} or \code{confint} methods to the results of
#'     \code{\link{as.glht}}. In the context of pairwise comparisons or comparisons
#'     with a control, this produces \dQuote{exact} Tukey or Dunnett adjustments,
#'     respectively. However, the algorithm (from the \pkg{mvtnorm} package) uses a
#'     Monte Carlo method, so results are not exactly repeatable unless the same
#'     random-number seed is used (see \code{\link[base:Random]{set.seed}}). As the family
#'     size increases, the required computation time will become noticeable or even
#'     intolerable, making the \code{"tukey"}, \code{"dunnettx"}, or others more
#'     attractive.}
#'   \item{\code{"none"}}{Makes no adjustments to the \eqn{p} values.}
#'   } %%%%%%%%%%%%%%%% end \describe {}
#' 
#'   For tests, not confidence intervals, the Bonferroni-inequality-based adjustment
#'   methods in \code{\link{p.adjust}} are also available (currently, these
#'   include \code{"holm"}, \code{"hochberg"}, \code{"hommel"},
#'   \code{"bonferroni"}, \code{"BH"}, \code{"BY"}, \code{"fdr"}, and
#'   \code{"none"}). If a \code{p.adjust.methods} method other than
#'   \code{"bonferroni"} or \code{"none"} is specified for confidence limits, the
#'   straight Bonferroni adjustment is used instead. Also, if an adjustment method
#'   is not appropriate (e.g., using \code{"tukey"} with one-sided tests, or with
#'   results that are not pairwise comparisons), a more appropriate method
#'   (usually \code{"sidak"}) is substituted.
#' 
#'   In some cases, confidence and \eqn{p}-value adjustments are only approximate
#'   -- especially when the degrees of freedom or standard errors vary greatly
#'   within the family of tests. The \code{"mvt"} method is always the correct
#'   one-step adjustment, but it can be very slow. One may use
#'   \code{\link{as.glht}} with methods in the \pkg{multcomp} package to obtain
#'   non-conservative multi-step adjustments to tests.
#'
#' @section Tests of significance, nonsuperiority, noninferiority, or equivalence:
#'   When \code{delta = 0}, test statistics are the usual tests of significance.
#'   They are of the form 
#'   \samp{(estimate - null)/SE}. Notationally: 
#'   \describe{
#'     \item{Significance}{\eqn{H_0: \theta = \theta_0}  versus \cr
#'        \eqn{H_1: \theta < \theta_0} (left-sided), or\cr
#'       \eqn{H_1 \theta > \theta_0} (right-sided), or\cr
#'       \eqn{H_1: \theta \ne \theta_0} (two-sided)\cr
#'       The test statistic is\cr
#'       \eqn{t = (Q - \theta_0)/SE}\cr 
#'       where \eqn{Q} is our estimate of \eqn{\theta};
#'       then left, right, or two-sided \eqn{p} values are produced, 
#'       depending on \code{side}.}
#'   }
#'   When \code{delta} is positive, the test statistic depends on \code{side} as
#'   follows.
#'   \describe{
#'   \item{Left-sided (nonsuperiority)}{\eqn{H_0: \theta \ge \theta_0 + \delta}
#'     versus \eqn{H_1: \theta < \theta_0 + \delta}\cr 
#'     \eqn{t = (Q - \theta_0 - \delta)/SE}\cr 
#'     The \eqn{p} value is the lower-tail probability.}
#'   \item{Right-sided (noninferiority)}{\eqn{H_0: \theta \le \theta_0 - \delta}
#'     versus \eqn{H_1: \theta > \theta_0 - \delta}\cr 
#'     \eqn{t = (Q - \theta_0 + \delta)/SE}\cr
#'     The \eqn{p} value is the upper-tail probability.}
#'   \item{Two-sided (equivalence)}{\eqn{H_0: |\theta - \theta_0| \ge \delta}
#'     versus \eqn{H_1: |\theta - \theta_0| < \delta}\cr
#'     \eqn{t = (|Q - \theta_0| - \delta)/SE}\cr
#'     The \eqn{p} value is the \emph{lower}-tail probability.\cr
#'     Note that \eqn{t} is the maximum of \eqn{t_{nonsup}} and \eqn{-t_{noninf}}. 
#'     This is equivalent to choosing the less 
#'     significant result in the two-one-sided-test (TOST) procedure.}
#'   } %%%%%%%%%%%% end \describe{}
#'
#' 
#' @section Non-estimable cases:
#'   When the model is rank-deficient, each row \code{x} of \code{object}'s
#'   \code{linfct} slot is checked for estimability. If \code{sum(x*bhat)}
#'   is found to be non-estimable, then the string \code{NonEst} is displayed for the
#'   estimate, and associated statistics are set to \code{NA}. 
#'   The estimability check is performed
#'   using the orthonormal basis \code{N} in the \code{nbasis} slot for the null
#'   space of the rows of the model matrix. Estimability fails when
#'   \eqn{||Nx||^2 / ||x||^2} exceeds \code{tol}, which by default is
#'   \code{1e-8}. You may change it via \code{\link{emm_options}} by setting
#'   \code{estble.tol} to the desired value.
#'   
#' @section Warning about potential misuse of P values:
#'   A growing consensus in the statistical and scientific community is that
#'   the term \dQuote{statistical significance} should be completely abandoned, and
#'   that criteria such as \dQuote{p < 0.05} never be used to assess the
#'   importance of an effect. These practices are just too misleading and prone to abuse.
#'   See \href{../doc/basics.html#pvalues}{the \dQuote{basics} vignette} for more
#'   discussion.
#'   
#' 
#' @note In doing testing and a transformation and/or link is in force, any
#'   \code{null} and/or \code{delta} values specified must always be on the
#'   scale of the linear predictor, regardless of the setting for `type`. If
#'   \code{type = "response"}, the null value displayed in the summary table 
#'   will be back-transformed from the value supplied by the user. But the
#'   displayed \code{delta} will not be changed, because there (usually) is
#'   not a natural way to back-transform it.
#' 
#' @note The default \code{show} method for \code{emmGrid} objects (with the
#'   exception of newly created reference grids) is \code{print(summary())}.
#'   Thus, with ordinary usage of \code{\link{emmeans}} and such, it is
#'   unnecessary to call \code{summary} unless there is a need to
#'   specify other than its default options.
#'   
#' @seealso \code{\link{hpd.summary}}
#' 
#' @method summary emmGrid  
#' @export
#'
#' @examples
#' warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
#' warp.emm <- emmeans(warp.lm, ~ tension | wool)
#' warp.emm    # implicitly runs 'summary'
#' 
#' confint(warp.emm, by = NULL, level = .90)
#' 
#' # --------------------------------------------------------------
#' pigs.lm <- lm(log(conc) ~ source + factor(percent), data = pigs)
#' pigs.emm <- emmeans(pigs.lm, "percent", type = "response")
#' summary(pigs.emm)    # (inherits type = "response")
#' summary(pigs.emm, calc = c(n = ".wgt."))  # Show sample size
#' 
#' # For which percents is EMM non-inferior to 35, based on a 10% threshold?
#' # Note the test is done on the log scale even though we have type = "response"
#' test(pigs.emm, null = log(35), delta = log(1.10), side = ">")
#' 
#' con <- contrast(pigs.emm, "consec")
#' test(con)
#' 
#' test(con, joint = TRUE)
#' 
#' # default Scheffe adjustment - rank = 3
#' summary(con, infer = c(TRUE, TRUE), adjust = "scheffe")
#' 
#' # Consider as some of many possible contrasts among the six cell means
#' summary(con, infer = c(TRUE, TRUE), adjust = "scheffe", scheffe.rank = 5)
#'
summary.emmGrid <- function(object, infer, level, adjust, by, type, df, calc,
                        null, delta, side, frequentist, 
                        bias.adjust = get_emm_option("back.bias.adj"),
                        sigma, ...) {
    if(missing(sigma))
        sigma = object@misc$sigma
    if(missing(frequentist) && !is.null(object@misc$frequentist))
        frequentist = object@misc$frequentist
    if(missing(bias.adjust)) {
        if (!is.null(object@misc$bias.adjust)) 
            bias.adjust = object@misc$bias.adjust
        else
            bias.adjust = get_emm_option("back.bias.adj")
    }
    
    if(!is.na(object@post.beta[1]) && (missing(frequentist) || !frequentist))
        return (hpd.summary(object, prob = level, by = by, type = type, 
                            bias.adjust = bias.adjust, sigma = sigma, ...))
    
    # Any "summary" options override built-in
    opt = get_emm_option("summary")
    if(!is.null(opt)) {
        opt$object = object
        object = do.call("update.emmGrid", opt)
    }
    
    if(!missing(sigma))
        object = update(object, sigma = sigma)  # we'll keep sigma in misc
    
    misc = object@misc
    use.elts = .reconcile.elts(object)
    grid = object@grid[use.elts, , drop = FALSE]
    
    ### For missing arguments, get from misc, else default    
    if(missing(infer))
        infer = misc$infer
    if(missing(level))
        level = misc$level
    if(missing(adjust))
        adjust = misc$adjust
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
    
    if(missing(df)) 
        df = misc$df
    if(!is.null(df)) {
        object@dffun = function(k, dfargs) df[1]
        attr(object@dffun, "mesg") = "user-specified"
    }
    
    # for missing args that default to zero unless provided or in misc slot
    .nul.eq.zero = function(val) {
        if(is.null(val)) 0
        else val
    }
    if(missing(null))
        null = .nul.eq.zero(misc$null)
    if(missing(delta))
        delta = .nul.eq.zero(misc$delta)
    if(missing(side))
        side = .nul.eq.zero(misc$side)
    
    
    # reconcile all the different ways we could specify the alternative
    # ... and map each to one of the first 3 subscripts
    side.opts = c("left","both","right","two-sided","noninferiority","nonsuperiority","equivalence","superiority","inferiority","0","2","-1","1","+1","<",">","!=","=")
    side.map =  c( 1,     2,     3,      2,          3,               1,               2,            3,            1,            2,  2,   1,  3,   3,  1,  3,  2,   2)
    side = side.map[pmatch(side, side.opts, 2)[1]] - 2
    delta = abs(delta)
    
    result = .est.se.df(object, ...)
    
    lblnms = setdiff(names(grid), 
                     c(object@roles$responses, ".offset.", ".wgt."))
    lbls = grid[lblnms]
    
    zFlag = (all(is.na(result$df) | is.infinite(result$df)))
    inv = (type %in% c("response", "mu", "unlink")) # flag to inverse-transform
    link = attr(result, "link")
    if (inv && is.null(link))
        inv = FALSE
    
    if ((length(infer) == 0) || !is.logical(infer)) 
        infer = c(FALSE, FALSE)
    if(length(infer == 1)) 
        infer = c(infer,infer)
    
    if(inv && !is.null(misc$tran)) {
        if (!is.null(misc$inv.lbl)) {
            names(result)[1] = misc$inv.lbl
            if (!is.null(misc$log.contrast))  # contrast of logs - relabel as ratios
                for (ell in seq_along(lbls)){
                    lbls[[ell]] = factor(lbls[[ell]])
                    levels(lbls[[ell]]) = gsub(" - ", " / ", levels(lbls[[ell]]))
                }
        }
        else
            names(result)[1] = "response"
    }
    
    attr(result, "link") = NULL
    estName = names(result)[1]
    
    mesg = misc$initMesg
    
    # Look for "mesg" attribute in dffun
    if (!is.null(dfm <- attr(object@dffun, "mesg")))
        mesg = c(mesg, paste("Degrees-of-freedom method:", dfm))
    
    ### Add an annotation when we show results on lp scale and
    ### there is a transformation
    if (!inv && !is.null(link)) {
        mesg = c(mesg, paste("Results are given on the", link$name, "(not the response) scale."))
    }
    if (inv && !is.null(link$unknown)) {
        mesg = c(mesg, paste0('Unknown transformation "', link$name, '": no transformation done'))
        inv = FALSE
        link = NULL
    }
    if(inv && bias.adjust && !is.null(link)) 
        link = .make.bias.adj.link(link, sigma)
    
    # et = 1 if a prediction, 2 if a contrast (or unmatched or NULL), 3 if pairs
    et = pmatch(c(misc$estType, "c"), c("prediction", "contrast", "pairs"), nomatch = 2)[1]
    
    by.size = nrow(grid)
    by.rows = .find.by.rows(grid, by)
    if (!is.null(by)) {
        if (length(unique(sapply(by.rows, length))) > 1) {
            by.size = misc$famSize = "(varies)"
        }
        else for (nm in by)
            by.size = by.size / length(unique(object@levels[[nm]]))
    }
    fam.info = c(misc$famSize, by.size, et)
    cnm = NULL
    adjust = tolower(adjust)
    
    # get vcov matrix only if needed (adjust == "mvt")
    corrmat = sch.rank = NULL
    if (!is.na(pmatch(adjust, "mvt"))) {
        corrmat = cov2cor(vcov(object))
        attr(corrmat, "by.rows") = by.rows
    }
    else if (!is.na(pmatch(adjust, "scheffe"))) {
        if(is.null(sch.rank <- .match.dots("scheffe.rank", ...)))
            sch.rank = sapply(by.rows, function(.) qr(zapsmall(object@linfct[., , drop = FALSE]))$rank)
        if(length(unique(sch.rank)) > 1)
            fam.info[1] = "uneven"   # This forces ragged.by = TRUE in .adj functions
    }
    
    # Add calculated columns
    if(!missing(calc) || !is.null(calc <- misc$calc)) {
        env = c(result, grid[setdiff(names(grid), names(result))])
        for (v in names(calc)) {
            elt = rev(as.character(calc[[v]]))[1] # pick out rhs if a formula
            val = try(eval(parse(text = elt), envir = env), silent = TRUE)
            if(!inherits(val, "try-error"))
                result[[v]] = env[[v]] = val
            else
                warning("The column '", v, "' could not be calculated, ",
                " so it is omitted")
        }
    }

    if(infer[1]) { # add CIs
        acv = .adj.critval(level, result$df, adjust, fam.info, side, corrmat, by.rows, sch.rank)
        ###adjust = acv$adjust # in older versions, I forced same adj method for tests
        cv = acv$cv
        cv = switch(side + 2, cbind(-Inf, cv), cbind(-cv, cv), cbind(-cv, Inf))
        cnm = if (zFlag) c("asymp.LCL", "asymp.UCL") else c("lower.CL","upper.CL")
        if(!is.null(misc$.predFlag)) {
            cnm = c("lower.PL", "upper.PL")
            sigma = misc$sigma
            mesg = c(mesg, paste0(
                "Prediction intervals and SEs are based on an error SD of ", 
                .fmt.sigma(sigma)))
            estName = names(result)[1] = "prediction"
        }
        result[[cnm[1]]] = result[[1]] + cv[, 1]*result$SE
        result[[cnm[2]]] = result[[1]] + cv[, 2]*result$SE
        mesg = c(mesg, paste("Confidence level used:", level), acv$mesg)
        if (inv) {
            clims = with(link, cbind(linkinv(result[[cnm[1]]]), linkinv(result[[cnm[2]]])))
            tmp = apply(clims, 1, function(x) { 
                z = diff(x); ifelse(is.na(z), 0, z) })
            idx = if (all(tmp >= 0)) 1:2 else 2:1
            result[[cnm[1]]] = clims[, idx[1]]
            result[[cnm[2]]] = clims[, idx[2]]
            mesg = c(mesg, paste("Intervals are back-transformed from the", link$name, "scale"))
        }
    }
    if(infer[2]) { # add tests
        if (!all(null == 0)) {
            result[["null"]] = null
            if (inv && !is.null(link))
                result[["null"]] = link$linkinv(null)
        }
        tnm = ifelse (zFlag, "z.ratio", "t.ratio")
        tail = ifelse(side == 0, -sign(abs(delta)), side)
        if (side == 0) {
            if (delta == 0) # two-sided sig test
                t.ratio = result[[tnm]] = (result[[1]] - null) / result$SE
            else
                t.ratio = result[[tnm]] = (abs(result[[1]] - null) - delta) / result$SE
        }
        else {
            t.ratio = result[[tnm]] = (result[[1]] - null + side * delta) / result$SE            
        }
        apv = .adj.p.value(t.ratio, result$df, adjust, fam.info, tail, corrmat, by.rows, sch.rank)
        adjust = apv$adjust   # in case it was abbreviated
        result$p.value = apv$pval
        mesg = c(mesg, apv$mesg)
        if (delta > 0)
            mesg = c(mesg, paste("Statistics are tests of", c("nonsuperiority","equivalence","noninferiority")[side+2],
                                 "with a threshold of", signif(delta, 5)))
        if(tail != 0) 
            mesg = c(mesg, paste("P values are ", ifelse(tail<0,"left-","right-"),"tailed", sep=""))
        if (inv) 
            mesg = c(mesg, paste("Tests are performed on the", link$name, "scale"))
    }
    if (inv) {
        result[["SE"]] = with(link, abs(mu.eta(result[[1]]) * result[["SE"]]))
        result[[1]] = with(link, linkinv(result[[1]]))
        if(bias.adjust)
            mesg = c(mesg, paste("Bias adjustment applied based on sigma =", 
                                 .fmt.sigma(sigma)))
    }
    
    if (length(misc$avgd.over) > 0) {
        qual = attr(misc$avgd.over, "qualifier")
        if (is.null(qual)) qual = ""
        mesg = c(paste0("Results are averaged over", qual, " the levels of: ",
                        paste(misc$avgd.over, collapse = ", ")), mesg)
    }
    summ = cbind(lbls, result)
    attr(summ, "estName") = estName
    attr(summ, "clNames") = cnm  # will be NULL if infer[1] is FALSE
    if (is.null(misc$pri.vars) || length(misc$pri.vars) == 0)
        misc$pri.vars = names(object@levels)
    attr(summ, "pri.vars") = setdiff(union(misc$pri.vars, misc$by.vars), c(by, ".wgt.", ".offset."))
    attr(summ, "by.vars") = by
    attr(summ, "adjust") = adjust
    attr(summ, "side") = side
    attr(summ, "delta") = delta
    attr(summ, "type") = type
    attr(summ, "mesg") = unique(mesg)
    class(summ) = c("summary_emm", "data.frame")
    summ
}

# S3 predict method

#' @rdname summary.emmGrid
#' @method predict emmGrid
#' @param interval Type of interval desired (partial matching is allowed): 
#' \code{"none"} for no intervals,
#'   otherwise confidence or prediction intervals with given arguments, 
#'   via \code{\link{confint.emmGrid}}.
#'   
#' @export
#' @return \code{predict} returns a vector of predictions for each row of \code{object@grid}.
predict.emmGrid <- function(object, type, 
                            interval = c("none", "confidence", "prediction"),
                            level = 0.95,
                            bias.adjust = get_emm_option("back.bias.adj"), sigma, 
                            ...) 
{
    # update with any "summary" options
    opt = get_emm_option("summary")
    if(!is.null(opt)) {
        opt$object = object
        object = do.call("update.emmGrid", opt)
    }

    interval = match.arg(interval)
    if (interval %in% c( "confidence", "prediction")) {
        if (interval == "prediction")
            object@misc$.predFlag = TRUE
        return(confint.emmGrid(object, type = type, level = level, 
                               bias.adjust = bias.adjust, sigma = sigma, ...))
    }
    
    if (missing(type))
        type = .get.predict.type(object@misc)
    else
        type = .validate.type(type)
    
    # if there are two transformations and we want response, then we need to undo both
    if ((type == "response") && (!is.null(object@misc$tran2)))
        object = regrid(object, transform = "mu")
    
    pred = .est.se.df(object, do.se = FALSE, ...)
    result = pred[[1]]
    
    if (type %in% c("response", "mu", "unlink")) {
        link = attr(pred, "link")
        if (!is.null(link)) {
            if (bias.adjust) {
                if(missing(sigma))
                    sigma = object@misc$sigma
                link = .make.bias.adj.link(link, sigma)
            }
            result = link$linkinv(result)
            if (is.logical(link$unknown) && link$unknown)
                warning("Unknown transformation: \"", link$name, "\" -- no transformation applied.")
        }
    }
    result
}

# as.data.frame method

#' @rdname summary.emmGrid
#' @param x object of the given class
#' @param row.names passed to \code{\link{as.data.frame}}
#' @param optional passed to \code{\link{as.data.frame}}
#' @return The \code{as.data.frame} method returns a plain data frame,
#'   equivalent to \code{as.data.frame(summary(.))}.
#' @method as.data.frame emmGrid
#' @export
as.data.frame.emmGrid = function(x, row.names = NULL, optional = FALSE, ...) {
    as.data.frame(summary(x, ...), row.names = row.names, optional = optional)
}


#' @rdname summary.emmGrid
#' @method [ summary_emm
#' @param as.df Logical value. With \code{x[..., as.df = TRUE]}, the result is
#'   object is coerced to an ordinary \code{\link{data.frame}}; otherwise, it is left as a 
#'   \code{summary_emm} object.
#' @export
"[.summary_emm" = function(x, ..., as.df = TRUE) {
    if (as.df)
        as.data.frame(x)[...]
    else
        base::`[.data.frame`(x, ...)
}



# Computes the quadratic form y'Xy after subsetting for the nonzero elements of y
.qf.non0 = function(X, y) {
    ii = (zapsmall(y) != 0)
    if (any(ii))
        sum(y[ii] * (X[ii, ii, drop = FALSE] %*% y[ii]))
    else 0
}


# Workhorse for getting estimates etc.
.est.se.df = function(object, do.se=TRUE, tol = get_emm_option("estble.tol"), ...) {
    if (nrow(object@grid) == 0) {
        result = data.frame(NA, NA, NA)
        names(result) = c(object@misc$estName, "SE", "df")
        return(result[-1, ])
    }
    misc = object@misc
    use.elts = .reconcile.elts(object)

    if (!is.null(hook <- misc$estHook)) {
        if (is.character(hook)) hook = get(hook)
        result = hook(object, do.se=do.se, tol=tol, ...)
    }
    else {
        active = which(!is.na(object@bhat))
        bhat = object@bhat[active]
        result = t(apply(object@linfct[use.elts, , drop = FALSE], 1, function(x) {
            if (!any(is.na(x)) && estimability::is.estble(x, object@nbasis, tol)) {
                x = x[active]
                est = sum(bhat * x)
                if(do.se) {
                    se = sqrt(.qf.non0(object@V, x))
                    df = object@dffun(x, object@dfargs)
                    ### Brute force: if (is.na(df)) df = Inf    ### Added
                }
                else # if these unasked-for results are used, we're bound to get an error!
                    se = df = 0
                c(est, se, df)
            }
            else c(NA,NA,NA)
        }))
            
        if (!is.null(object@grid$.offset.))
            result[, 1] = result[, 1] + object@grid$.offset.[use.elts]
    }
    result[1] = as.numeric(result[1]) # silly bit of code to avoid getting a data.frame of logicals if all are NA
    result = as.data.frame(result)
    names(result) = c(misc$estName, "SE", "df")
    if (!is.null(misc$.predFlag)) {
        if (is.null(misc$sigma))
            stop("No 'sigma' is available for obtaining Prediction intervals.", 
                 call. = FALSE)
        result$SE = sqrt(result$SE^2 + misc$sigma^2)
    }
    if (!is.null(misc$tran)) {
        attr(result, "link") = .get.link(misc)
        if(is.character(misc$tran) && (misc$tran == "none"))
            attr(result, "link") = NULL
    }
    result
}

# workhorse for nailing down the link function
.get.link = function(misc) {
    link = if(is.character(misc$tran))
        .make.link(misc$tran)
    else if (is.list(misc$tran))
        misc$tran
    else 
        NULL
    
    if (is.list(link)) {  # See if multiple of link is requested, or variable is offset
        if (!is.null(misc$tran.mult) || !is.null(misc$tran.offset)) {
            name = link$name
            mult = link$mult = ifelse(is.null(misc$tran.mult), 1, misc$tran.mult)
            off = link$offset = ifelse(is.null(misc$tran.offset), 0, misc$tran.offset)
            link = with(link, list(
                linkinv = function(eta) linkinv(eta / mult) - offset,
                mu.eta = function(eta) mu.eta(eta / mult) / mult))
            if(mult != 1)
                name =  paste0(round(mult, 3), "*", name)        
            if(off != 0)
                name = paste0(name, "(mu + ", round(off, 3), ")")
            link$name = name         
        }
    }
    
    if (!is.null(link) && is.null(link$name))
        link$name = "linear-predictor"
    link
}

# patch-in alternative back-transform stuff for bias adjustment
# Currently, we just use a 2nd-order approx for everybody:
#   E(h(nu + E))  ~=  h(nu) + 0.5*h"(nu)*var(E)
.make.bias.adj.link = function(link, sigma) {
    if (is.null(sigma))
        stop("Must specify 'sigma' to obtain bias-adjusted back transformations", call. = FALSE)
    link$inv = link$linkinv
    link$der = link$mu.eta
    link$sigma22 = sigma^2 / 2
    link$der2 = function(eta) with(link, 1000 * (der(eta + .0005) - der(eta - .0005)))
    link$linkinv = function(eta) with(link, inv(eta) + sigma22 * der2(eta))
    link$mu.eta = function(eta) with(link, der(eta) +
                                         1000 * sigma22 * (der2(eta + .0005) - der2(eta - .0005)))
    link
}
####!!!!! TODO: Re-think whether we are handling Scheffe adjustments correctly
####!!!!!       if/when we shift around 'by' specs, etc.

# utility to compute an adjusted p value
# tail is -1, 0, 1 for left, two-sided, or right
# Note fam.info is c(famsize, ncontr, estTypeIndex)
# lsmeans >= 2.14: added corrmat arg, dunnettx & mvt adjustments
# emmeans > 1.3.4: we have sch.rank of same length as by.rows
# NOTE: corrmat is NULL unless adjust == "mvt"
.adj.p.value = function(t, DF, adjust, fam.info, tail, corrmat, by.rows, sch.rank) {
    fam.size = fam.info[1]
    n.contr = fam.info[2]
    et = as.numeric(fam.info[3])

    if (n.contr == 1) # Force no adjustment when just one test
        adjust = "none"
    
    # do a pmatch of the adjust method, case insensitive
    adj.meths = c("sidak", "tukey", "scheffe", "dunnettx", "mvt", p.adjust.methods)
    k = pmatch(tolower(adjust), adj.meths)
    if(is.na(k))
        stop("Adjust method '", adjust, "' is not recognized or not valid")
    adjust = adj.meths[k]
    if ((tail != 0) && (adjust %in% c("tukey", "scheffe", "dunnettx"))) # meth not approp for 1-sided
        adjust = "sidak"
    if ((et != 3) && adjust == "tukey") # not pairwise
        adjust = "sidak"
    
    ragged.by = (is.character(fam.size) || adjust %in% c("mvt", p.adjust.methods))   # flag that we need to do groups separately
    if (!ragged.by)
        by.rows = list(seq_along(t))       # not ragged, we can do all as one by group
    
    # asymptotic results when df is NA
    DF[is.na(DF)] = Inf
    
    # if estType is "prediction", use #contrasts + 1 as family size
    # (produces right Scheffe CV; Tukey ones are a bit strange)
    # deprecated - used to try to keep track of scheffe.adj = ifelse(et == 1, 0, - 1)
    if (tail == 0)
        p.unadj = 2*pt(abs(t), DF, lower.tail=FALSE)
    else
        p.unadj = pt(t, DF, lower.tail = (tail<0))
    
#    pvals = lapply(by.rows, function(rows) {
    pval = numeric(length(t))
    for(jj in seq_along(by.rows)) { ####(rows in by.rows) {
        rows = by.rows[[jj]]
        unadj.p = p.unadj[rows]
        abst = abs(t[rows])
        df = DF[rows]
        if (ragged.by) {
            n.contr = length(rows)
            fam.size = (1 + sqrt(1 + 8*n.contr)) / 2   # tukey family size - e.g., 6 pairs -> family of 4
        }
        if (adjust %in% p.adjust.methods) {
            if (n.contr == length(unadj.p))
                pval[rows] = p.adjust(unadj.p, adjust, n = n.contr)
            else # only will happen when by.rows is length 1
                pval[rows] = as.numeric(apply(matrix(unadj.p, nrow=n.contr), 2, 
                                        function(pp) p.adjust(pp, adjust, n=sum(!is.na(pp)))))
        }
        else pval[rows] = switch(adjust,
                           sidak = 1 - (1 - unadj.p)^n.contr,
                           # NOTE: tukey, scheffe, dunnettx all assumed 2-sided!
                           tukey = ptukey(sqrt(2)*abst, fam.size, zapsmall(df), lower.tail=FALSE),
                           scheffe = pf(t[rows]^2 / (sch.rank[jj]), sch.rank[jj], 
                                        df, lower.tail = FALSE),
                           dunnettx = 1 - .pdunnx(abst, n.contr, df),
                           mvt = 1 - .my.pmvt(t[rows], df, corrmat[rows,rows,drop=FALSE], -tail) # tricky - reverse the tail because we're subtracting from 1 
        )
    }

    chk.adj = match(adjust, c("none", "tukey", "scheffe"), nomatch = 99)
    
    if (ragged.by) {
        nc = max(sapply(by.rows, length))
        fs = (1 + sqrt(1 + 8*nc)) / 2
        scheffe.dim = "(varies)"
    }
    else {
        nc = n.contr
        fs = fam.size
        scheffe.dim = sch.rank[1]
    }
    do.msg = (chk.adj > 1) && (nc > 1) && !((fs < 3) && (chk.adj < 10)) 
    if (do.msg) {
#         xtra = if(chk.adj < 10) paste("a family of", fam.size, "tests")
#         else             paste(n.contr, "tests")
        xtra = switch(adjust, 
                      tukey = paste("for comparing a family of", fam.size, "estimates"),
                      scheffe = paste("with rank", scheffe.dim),
                      paste("for", n.contr, "tests")
                )
        mesg = paste("P value adjustment:", adjust, "method", xtra)
    }
    else mesg = NULL
    list(pval = pval, mesg = mesg, adjust = adjust)
}

# Code needed for an adjusted critical value
# returns a list similar to .adj.p.value
# lsmeans >= 2.14: Added tail & corrmat args, dunnettx & mvt adjustments
# emmeans > 1.3.4: Added sch.rank
# NOTE: corrmat is NULL unless adjust == "mvt"
.adj.critval = function(level, DF, adjust, fam.info, tail, corrmat, by.rows, sch.rank) {
    mesg = NULL
    
    fam.size = fam.info[1]
    n.contr = fam.info[2]
    et = as.numeric(fam.info[3])
    
    ragged.by = (is.character(fam.size))   # flag that we need to do groups separately
    if (!ragged.by && n.contr == 1) # Force no adjustment when just one interval
        adjust = "none"
    
    adj.meths = c("sidak", "tukey", "scheffe", "dunnettx", "mvt", "bonferroni", "none")
    k = pmatch(tolower(adjust), adj.meths)
    if(is.na(k))
        k = which(adj.meths == "bonferroni") 
    adjust = adj.meths[k]
    
    if (!ragged.by && (adjust != "mvt") && (length(unique(DF)) == 1))
        by.rows = list(seq_along(DF))       # not ragged, we can do all as one by group
    
    if ((et != 3) && adjust == "tukey") # not pairwise
        adjust = "sidak"
    if ((tail != 0) && (adjust %in% c("tukey", "scheffe", "dunnettx"))) # meth not approp for 1-sided
        adjust = "sidak"
    if ((et != 3) && adjust == "tukey") # not pairwise
        adjust = "sidak"
    
    # asymptotic results when df is NA
    DF[is.na(DF)] = Inf
    #### No longer used scheffe.adj = ifelse(et == 1, 0, - 1)
    
    chk.adj = match(adjust, c("none", "tukey", "scheffe"), nomatch = 99)
    if (ragged.by) {
        nc = max(sapply(by.rows, length))
        fs = (1 + sqrt(1 + 8*nc)) / 2
        scheffe.dim = "(varies)"
    }
    else {
        nc = n.contr
        fs = fam.size
        scheffe.dim = sch.rank[1]
    }
    do.msg = (chk.adj > 1) && (nc > 1) && 
        !((fs < 3) && (chk.adj < 10)) 
    
    if (do.msg) {
#        xtra = if(chk.adj < 10) paste("a family of", fam.size, "estimates")
#        else             paste(n.contr, "estimates")
        xtra = switch(adjust, 
                      tukey = paste("for comparing a family of", fam.size, "estimates"),
                      scheffe = paste("with rank", scheffe.dim),
                      paste("for", n.contr, "estimates")
        )
        mesg = paste("Conf-level adjustment:", adjust, "method", xtra)
    }
    
    adiv = ifelse(tail == 0, 2, 1) # divisor for alpha where needed
    
    ###cvs = lapply(by.rows, function(rows) {
    cv = numeric(sum(sapply(by.rows, length)))
    for (jj in seq_along(by.rows)) { ####(rows in by.rows) {
        rows = by.rows[[jj]]
        df = DF[rows]
        if (ragged.by) {
            n.contr = length(rows)
            fam.size = (1 + sqrt(1 + 8*n.contr)) / 2   # tukey family size - e.g., 6 pairs -> family of 4
        }
        cv[rows] = switch(adjust,
               none = -qt((1-level)/adiv, df),
               sidak = -qt((1 - level^(1/n.contr))/adiv, df),
               bonferroni = -qt((1-level)/n.contr/adiv, df),
               tukey = qtukey(level, fam.size, df) / sqrt(2),
               scheffe = sqrt((sch.rank[jj]) * qf(level, sch.rank[jj], df)),
               dunnettx = .qdunnx(level, n.contr, df),
               mvt = .my.qmvt(level, df, corrmat[rows, rows, drop = FALSE], tail)
        )
    }
    
    list(cv = cv, mesg = mesg, adjust = adjust)
}


### My own functions to ease access to mvt functions
### These use one argument at a time and expands each (lower, upper) or p to a k-vector
### Use tailnum = -1, 0, or 1
### NOTE: corrmat needs "by.rows" attribute to tell which rows
###   belong to which submatrix.
.my.pmvt = function(x, df, corrmat, tailnum) {
    lower = switch(tailnum + 2, -Inf, -abs(x), x)
    upper = switch(tailnum + 2, x, abs(x), Inf)
    by.rows = attr(corrmat, "by.rows")
    if (is.null(by.rows)) 
        by.rows = list(seq_len(length(x)))
    by.sel = numeric(length(x))
    for (i in seq_along(by.rows))
        by.sel[by.rows[[i]]] = i
    df = .fix.df(df)
    apply(cbind(lower, upper, df, by.sel), 1, function(z) {
        idx = by.rows[[z[4]]]
        k = length(idx)
        pval = try(mvtnorm::pmvt(rep(z[1], k), rep(z[2], k), 
                        df = as.integer(z[3]), corr = corrmat[idx, idx]), 
                    silent = TRUE)
        if (inherits(pval, "try-error"))   NA
        else                               pval
    })
}

# Vectorized for df but needs p to be scalar
.my.qmvt = function(p, df, corrmat, tailnum) {
    tail = c("lower.tail", "both.tails", "lower.tail")[tailnum + 2] 
    df = .fix.df(df)
    by.rows = attr(corrmat, "by.rows")
    if (is.null(by.rows)) 
        by.rows = list(seq_len(length(df)))
    by.sel = numeric(length(df))
    for (i in seq_along(by.rows))
        by.sel[by.rows[[i]]] = i
    # If df all equal, compute just once for each by group
    eq.df = (diff(range(df)) == 0)
    i1 = if (eq.df)   sapply(by.rows, function(r) r[1])
         else         seq_along(df)
    result = apply(cbind(p, df[i1], by.sel[i1]), 1, function(z) {
        idx = by.rows[[z[3]]]
        cv = try(mvtnorm::qmvt(z[1], tail = tail, 
                    df = as.integer(z[2]), corr = corrmat[idx, idx])$quantile,
                 silent = TRUE)
        if (inherits(cv, "try-error"))     NA
        else                               cv
    })
    if (eq.df) {
        res = result
        result = numeric(length(df))
        for(i in seq_along(by.rows))
            result[by.rows[[i]]] = res[i]
    } 
    result
}

# utility to get appropriate integer df
.fix.df = function(df) {
    sapply(df, function(d) {
        if (d > 0) d = max(1, d)
        if (is.infinite(d) || (d > 9999)) d = 0
        floor(d + .25) # tends to round down
    })
}

### My approximate dunnett distribution 
### - a mix of the Tukey cdf and Sidak-corrected t
.pdunnx = function(x, k, df, twt = (k - 1)/k) {
    tukey = ptukey(sqrt(2)*x, (1 + sqrt(1 + 8*k))/2, df)
    sidak = (pf(x^2, 1, df))^k
    twt*tukey + (1 - twt)*sidak
}

# Uses linear interpolation to get quantile
.qdunnx = function(p, k, df, ...) {
     if (k < 1.005)
         return(qt(1 - .5*(1 - p), df))
    xtuk = qtukey(p, (1 + sqrt(1 + 8*k))/2, df) / sqrt(2)
    xsid = sqrt(qf(p^(1/k), 1, df))
    fcn = function(x, d) 
        .pdunnx(x, k, d, ...) - p
    apply(cbind(xtuk, xsid, df), 1, function(r) {
        if (abs(diff(r[1:2])) < .0005)
            return (r[1])
        x = try(uniroot(fcn, r[1:2], tol = .0005, d = r[3]), silent = TRUE)
        if (inherits(x, "try-error")) {
            warning("Root-finding failed; using qtukey approximation for Dunnett quantile")
            return(xtuk)
        }
        else
            x$root
    })
}



### Support for different prediction types ###

# Valid values for type arg or predict.type option
#   "link", "lp", "linear" are all legal but equivalent
#   "mu" and "response" are usually equivalent -- but in a GLM with a response transformation,
#      "mu" (or "unlink") would back-transform the link only, "response" would do both
.valid.types = c("link","lp","linear", "response", "mu", "unlink")

# get "predict.type" option from misc, and make sure it's legal
.get.predict.type = function(misc) {
    type = misc$predict.type
    if (is.null(type))
        .valid.types[1]
    else
        .validate.type(type)
}

# check a "type" arg to make it legal
# NOTE: if not matched, returns "link", i.e., no back-transformation will be done
.validate.type = function (type) {
    type = .valid.types[pmatch(type, .valid.types, 1)]
    if (length(type) > 1) {
        type = type[1]
        warning("You specified more than one prediction type. Only type = \"", type, "\" was used")
    }
    type
}





# left-or right-justify column labels for m depending on "l" or "R" in just
.just.labs = function(m, just) {
    nm = dimnames(m)[[2]]
    for (j in seq_len(length(nm))) {
        if(just[nm[j]] == "L") 
            nm[j] = format(nm[j], width = nchar(m[1,j]), just="left")
    }
    dimnames(m) = list(rep("", nrow(m)), nm)
    m
}

# Format a data.frame produced by summary.emmGrid
#' @method print summary_emm
#' @export
print.summary_emm = function(x, ..., digits=NULL, quote=FALSE, right=TRUE) {
    x.save = x
    for(i in which(sapply(x, is.matrix))) 
        x[[i]] = NULL   # hide matrices
    for (i in seq_along(names(x)))   # zapsmall the numeric columns
        if (is.numeric(x[[i]]))  
            x[[i]] = zapsmall(x[[i]])
    if (!is.null(x$df)) x$df = round(x$df, 2)
    if (!is.null(x$t.ratio)) 
        x$t.ratio = format(round(x$t.ratio, 3), nsmall = 3, sci = FALSE)
    if (!is.null(x$z.ratio)) 
        x$z.ratio = format(round(x$z.ratio, 3), nsmall = 3, sci = FALSE)
    if (!is.null(x$p.value)) {
        fp = x$p.value = format(round(x$p.value, 4), nsmall = 4, sci = FALSE)
        x$p.value[fp=="0.0000"] = "<.0001"
    }
    estn = attr(x, "estName")
    just = sapply(x, function(col) if(is.numeric(col)) "R" else "L")
    est = x[[estn]]
    if (get_emm_option("opt.digits") && is.null(digits)) {
        if (!is.null(x[["SE"]]))
            tmp = est + x[["SE"]] * cbind(rep(-2, nrow(x)), 0, 2)
        else if (!is.null(x[["lower.HPD"]]))
            tmp = x[, c("lower.HPD", estn, "upper.HPD"), drop = FALSE]
        else tmp = NULL
        if (!is.null(tmp))
            digits = max(apply(tmp, 1, .opt.dig))
    }
    if (any(is.na(est))) {
        x[[estn]] = format(est, digits=digits)
        x[[estn]][is.na(est)] = "nonEst"
    }
    xc = as.matrix(format.data.frame(x, digits=digits, na.encode=FALSE))
    m = apply(rbind(just, names(x), xc), 2, function(x) {
        w = max(sapply(x, nchar))
        if (x[1] == "R") format(x[-seq_len(2)], width = w, justify="right")
        else format(x[-seq_len(2)], width = w, justify="left")
    })
    if(!is.matrix(m)) m = t(as.matrix(m))
    by.vars = attr(x, "by.vars")
    if (is.null(by.vars)) {
        m = .just.labs(m, just)
        print(m, quote=FALSE, right=TRUE)
        cat("\n")
    }
    else { # separate listing for each by variable
        m = .just.labs(m[, setdiff(names(x), by.vars), drop = FALSE], just)
        pargs = unname(as.list(x[,by.vars, drop=FALSE]))
        pargs$sep = ", "
        lbls = do.call(paste, pargs)
        for (lb in unique(lbls)) {
            rows = which(lbls==lb)
            levs = paste(by.vars, "=", xc[rows[1], by.vars])
            cat(paste(paste(levs, collapse=", ")), ":\n", sep="")
            print(m[rows, , drop=FALSE], ..., quote=quote, right=right)
            cat("\n")
        }
    }
    
    msg = unique(attr(x, "mesg"))
    if (!is.null(msg))
        for (j in seq_len(length(msg))) cat(paste(msg[j], "\n"))
    
    invisible(x.save)
}

# determine optimum digits to display based on a conf or cred interval
# (but always at least 3)
.opt.dig = function(x) {
    z = range(x) / max(abs(x))
    dz = zapsmall(c(diff(z), z), digits = 8)[1]
    zz = round(1.51 - log(dz, 10))  # approx 1 - log(diff(z/3))
    zz[is.infinite(zz)] = 3   # we get z = Inf when SE is 0
    max(zz, 3, na.rm = TRUE)  
}

# Utility -- When misc$display present, reconcile which elements to use.
# Needed if we messed with previous nesting
.reconcile.elts = function(object) {
    display = object@misc$display
    nrows = nrow(object@grid)
    use.elts = rep(TRUE, nrows)
    if (!is.null(display) && (length(display) == nrows))  
        use.elts = display
    use.elts
}
