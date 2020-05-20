# Cohen's effect sizes


#' Calculate effect sizes and confidence bounds thereof
#' 
#' Standardized effect sizes are typically calculated using pairwise differences of estimates,
#' divided by the SD of the population providing the context for those effects.
#' This function calculates effect sizes from an \code{emmGrid} object,
#' and confidence intervals for them, accounting for uncertainty in both the estimated
#' effects and the population SD. 
#' 
#' Any \code{by} variables specified in \code{object} will remain in force in the returned
#' effects, unless overridden in the optional arguments.
#' 
#' For models having a single random effect, such as those fitted using
#' \code{\link{lm}}; in that case, the \code{stats::sigma} and
#' \code{stats::df.residual} functions may be useful for specifying \code{sigma}
#' and \code{edf}. For models with more than one random effect, \code{sigma} may
#' be based on some combination of the random-effect variances. 
#' 
#' Specifying \code{edf} can be rather unintuitive but is also relatively
#' uncritical; but the smaller the value, the wider the confidence intervals for
#' effect size. The value of \code{sqrt(2/edf)} can be interpreted as the
#' relative accuracy of \code{sigma}; for example, with \code{edf = 50},
#' \eqn{\sqrt(2/50) = 0.2}, meaning that \code{sigma} is accurate to plus or
#' minus 20 percent. Note in an example below, we tried two different \code{edf}
#' values as kind of a bracketing/sensitivity-analysis strategy. A value of
#' \code{Inf} is allowable, in which case you are assuming that \code{sigma} is
#' known exactly. Obviously, this narrows the confidence intervals for the
#' effect sizes -- unrealistically if in fact \code{sigma} is unknown.
#' 
#'
#' @param object an \code{\link[=emmGrid-class]{emmGrid}} object, 
#' typically one defining the EMMs to 
#' be contrasted. If instead, \code{class(object) == "emm_list"},
#' such as is produced by \code{emmeans(model, pairwise ~ treatment)},
#' a message is displayed; the contrasts already therein are used; and 
#' \code{method} is replaced by \code{"identity"}.
#' @param sigma numeric scalar, value of the population SD. 
#' @param edf numeric scalar that specifies the equivalent degrees of freedom
#'   for the \code{sigma}. This is a way of specifying the uncertainty in \code{sigma},
#'   in that we regard our estimate of \code{sigma^2} as being proportional to
#'   a chi-square random variable with \code{edf} degrees of freedom. (\code{edf} should
#'   not be confused with the \code{df} argument that may be passed via \code{...}
#'   to specify the degrees of freedom to use in \eqn{t} statistics and confidence intervals.)
#' @param method the contrast method to use to define the effects.
#'   This is passed to \code{\link{contrast}} after the elements of \code{object}
#'   are scaled.
#' @param ... Additional arguments passed to \code{contrast}
#'
#' @return an \code{\link[=emmGrid-class]{emmGrid}} object containing the effect sizes
#' 
#' @section Computation:
#' This function uses calls to \code{\link{regrid}} to put the estimated
#' marginal means (EMMs) on the log scale. Then an extra element is added to
#' this grid for the log of \code{sigma} and its standard error (where we assume
#' that \code{sigma} is uncorrelated with the log EMMs). Then a call to
#' \code{\link{contrast}} subtracts \code{log{sigma}} from each of the log EMMs,
#' yielding values of \code{log(EMM/sigma)}.
#' Finally, the results are re-gridded back to the original scale and the
#' desired contrasts are computed using \code{method}. In the log-scaling
#' part, we actually rescale the absolute values and keep track of the signs.
#' 
#' @note
#' The effects are always computed on the scale of the \emph{linear-predictor};
#' any response transformation or link function is completely ignored. If you
#' wish to base the effect sizes on the response scale, it is \emph{not} enough
#' to replace \code{object} with \code{regrid(object)}, because this
#' back-transformation changes the SD required to compute effect sizes. 
#' 
#' @note
#' \strong{Disclaimer:} There is substantial disagreement among practitioners on
#' what is the appropriate \code{sigma} to use in computing effect sizes; or,
#' indeed, whether \emph{any} effect-size measure is appropriate for some
#' situations. The user is completely responsible for specifying 
#' appropriate parameters (or for failing to do so).
#' 
#' @export
#'
#' @examples
#' fiber.lm <- lm(strength ~ diameter + machine, data = fiber)
#' 
#' emm <- emmeans(fiber.lm, "machine")
#' eff_size(emm, sigma = sigma(fiber.lm), edf = df.residual(fiber.lm))
#' 
#' # or equivalently:
#' eff_size(pairs(emm), sigma(fiber.lm), df.residual(fiber.lm), method = "identity")
#' 
#' 
#' ### Mixed model example:
#' if (require(nlme)) {
#' 
#'   Oats.lme <- lme(yield ~ Variety + factor(nitro), 
#'                   random = ~ 1 | Block / Variety,
#'                   data = Oats)
#'                   
#'   # Combine variance estimates
#'   VarCorr(Oats.lme)
#'   totSD <- sqrt(214.4724 + 109.6931 + 162.5590)
#'   # I figure edf is somewhere between 5 (Blocks df) and 51 (Resid df)
#'   
#'   emmV <- emmeans(Oats.lme, ~ Variety)
#'   print(eff_size(emmV, sigma = totSD, edf = 5))
#'   print(eff_size(emmV, sigma = totSD, edf = 51))
#' }
#' 
#' # Multivariate model for the same data:
#'  MOats.lm <- lm(yield ~ Variety, data = MOats)
#'  eff_size(emmeans(MOats.lm, "Variety"), 
#'           sigma = sqrt(mean(sigma(MOats.lm)^2)),   # RMS of sigma()
#'           edf = df.residual(MOats.lm))
#' 
#' # These results illustrate a sobering message that effect sizes are often
#' # not nearly as accurate as you may think.
eff_size = function(object, sigma, edf, method = "pairwise", ...) {
    if (inherits(object, "emm_list") && ("contrasts" %in% names(object))) {
        message("Since 'object' is a list, we are using the contrasts already present.")
        object = object$contrasts
        method = "identity"
    }
    SE.logsigma = sqrt(1 / (2 * edf))
 
    object = update(object, tran = NULL)
    object = regrid(object, transform = "response")
    
    # put on absolute scale
    emm = object@bhat
    sgn = sign(emm)
    emm = abs(emm)
    sgn[emm < 1e-15] = 0
    emm[sgn == 0] = 1
    object@bhat = emm
    
    # put on log scale and incorporate the SD estimate
    logobj = regrid(object, transform = "log")
    logobj@bhat = c(logobj@bhat, log(sigma))
    V = rbind(cbind(logobj@V, 0), 0)
    n = nrow(V)
    V[n,n] = SE.logsigma^2
    logobj@V = V
    logobj@levels = list(dummy = 1:n)
    logobj@grid = data.frame(dummy = 1:n)
    logobj@linfct = diag(n)
    logobj@misc$by = NULL
    
    con = contrast(logobj, "trt.vs.ctrlk", by = NULL)
    con = regrid(con, transform = "response")
    object@bhat = con@bhat * sgn
    object@V = con@V
    update(contrast(object, method, adjust = "none", ...), 
           estName = "effect.size", infer = c(TRUE, FALSE),
           initMesg = paste("sigma used for effect sizes:", signif(sigma, digits = 4)))
}
