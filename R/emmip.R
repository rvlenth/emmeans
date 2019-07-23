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

# emmip code - interaction plots

#' Interaction-style plots for estimated marginal means
#'
#' Creates an interaction plot of EMMs based on a fitted model and a simple
#' formula specification.
#' 
#' @export
emmip = function(object, formula, ...) {
    UseMethod("emmip")
}

# Our one method
#' @rdname emmip
#' @param object An object of class \code{emmGrid}, or a fitted model of a class
#'   supported by the \pkg{emmeans} package
#' @param formula Formula of the form 
#'   \code{trace.factors ~ x.factors | by.factors}. The EMMs are
#'   plotted against \code{x.factor} for each level of \code{trace.factors}.
#'   \code{by.factors} is optional, but if present, it determines separate
#'   panels. Each element of this formula may be a single factor in the model,
#'   or a combination of factors using the \code{*} operator.
#' @param type As in \code{\link{predict.emmGrid}}, this determines
#'   whether we want to inverse-transform the predictions
#'   (\code{type = "response"}) or not (any other choice). The default is
#'   \code{"link"}, unless the \code{"predict.type"} option is in force; see
#'   \code{\link{emm_options}}.
#' @param CIs Logical value. If \code{TRUE}, confidence intervals (or HPD intervals
#'   for Bayesian models) are added to the plot 
#'   (works only with \code{engine = "ggplot"}).
#' @param PIs Logical value. If \code{TRUE}, prediction intervals are added to the plot 
#'   (works only with \code{engine = "ggplot"}). If both \code{CIs} and
#'   \code{CIs} are \code{TRUE}, the prediction intervals will be somewhat
#'   longer, lighter, and thinner than the confidence intervals. Additional
#'   parameters to \code{\link{predict.emmGrid}} (e.g., \code{sigma}) may be passed via
#'   \code{...}. For Bayesian models, PIs require \code{frequentist = TRUE} and 
#'   a value for \code{sigma}.
#' @param engine Character value matching \code{"ggplot"} (default) or 
#'   \code{"lattice"}. The graphics engine to be used to produce the plot.
#'   These require, respectively, the \pkg{ggplot2} or \pkg{lattice} package to
#'   be installed.
#' @param pch The plotting characters to use for each group (i.e., levels of
#'   \code{trace.factors}). They are recycled as needed.
#' @param lty The line types to use for each group. Recycled as needed.
#' @param col The colors to use for each group, recycled as needed. If not
#'   specified, the default trellis colors are used.
#' @param plotit Logical value. If \code{TRUE}, the plot is displayed.
#'   Otherwise, one may use the \code{"lattice"} attribute of the returned
#'   object and print it, perhaps after additional manipulation.
#' @param ... Additional arguments passed to \code{\link{emmeans}} (when
#'   \code{object} is not already an \code{emmGrid} object),
#'   \code{predict.emmGrid},
#'   \code{\link[ggplot2]{ggplot}}, or \code{\link[lattice]{xyplot}}.
#'   
#' @section Details:
#' If \code{object} is a fitted model, \code{\link{emmeans}} is called with an
#' appropriate specification to obtain estimated marginal means for each
#' combination of the factors present in \code{formula} (in addition, any 
#' arguments in \code{\dots} that match \code{at}, \code{trend}, 
#' \code{cov.reduce}, or \code{fac.reduce} are passed to \code{emmeans}). 
#' Otherwise, if \code{object} is an \code{emmGrid} object, its first element is 
#' used, and it must contain one estimate for each combination of the factors
#' present in \code{formula}.
#'
#' @return If \code{plotit = FALSE}, a \code{data.frame} (actually, a
#'   \code{summary_emm} object) with the table of EMMs that would be plotted.
#'   The variables plotted are named \code{xvar} and \code{yvar}, and the trace
#'   factor is named \code{tvar}. This data frame has an added \code{"labs"}
#'   attribute containing the labels \code{xlab}, \code{ylab}, and \code{tlab}
#'   for these respective variables. The confidence limits are also
#'   included, renamed \code{LCL} and \code{UCL}.
#'   
#' @return If \code{plotit = TRUE}, the function
#'   returns an object of class \code{"ggplot"} or a \code{"trellis"}, depending
#'   on \code{engine}.
#'   
#' @note Conceptually, this function is equivalent to 
#'   \code{\link{interaction.plot}} where the summarization function is thought 
#'   to return the EMMs.
#' 
#' @seealso \code{\link{emmeans}}, \code{\link{interaction.plot}}
#' @export
#' @method emmip default
#'
#' @examples
#' #--- Three-factor example
#' noise.lm = lm(noise ~ size * type * side, data = auto.noise)
#'
#' # Separate interaction plots of size by type, for each side
#' emmip(noise.lm, type ~ size | side)
#'
#' # One interaction plot, using combinations of size and side as the x factor
#' # ... with added confidence intervals
#' emmip(noise.lm, type ~ side * size, CIs = TRUE)
#'
#' # One interaction plot using combinations of type and side as the trace factor
#' emmip(noise.lm, type * side ~ size)
#'
#' # Individual traces in panels
#' emmip(noise.lm, ~ size | type * side)
#'
emmip.default = function(object, formula, type, CIs = FALSE, PIs = FALSE,
                         engine = get_emm_option("graphics.engine"),
                         pch = c(1,2,6,7,9,10,15:20), 
                         lty = 1, col = NULL, plotit = TRUE, ...) {
    engine = match.arg(engine, c("ggplot", "lattice"))
    if (engine == "ggplot")
        .requireNS("ggplot2",
                   "The 'ggplot' engine requires the 'ggplot2' package be installed.")
    if (engine == "lattice")
        .requireNS("lattice", 
                "The 'lattice' engine requires the 'lattice' package be installed.")
    
    specs = .parse.by.formula(formula) # list of lhs, rhs, by
    
    # Glean the parts of ... to use in emmeans call
    # arguments allowed to be passed
    lsa.allowed = c("at","trend","cov.reduce","fac.reduce")
    xargs = list(...)
    emmopts = list(...)
    for (arg in names(xargs)) {
        idx = pmatch(arg, lsa.allowed)
        if (!is.na(idx)) {
            opt = lsa.allowed[idx]
            emmopts[[opt]] = xargs[[arg]]
            xargs[[arg]] = NULL
        }
    }
    
    emmopts$object = object
    emmopts$specs = .reformulate(unlist(specs))
    emmo = do.call("emmeans", emmopts)
    if(missing(type)) {
        type = get_emm_option("summary")$predict.type
        if (is.null(type))
            type = .get.predict.type(emmo@misc)
    }
    type = .validate.type(type)
    
    emms = summary(emmo, type = type, infer = c(CIs, F))
    if(PIs) {
        prd = predict(emmo, interval = "pred", ...)
        emms$LPL = prd$lower.PL
        emms$UPL = prd$upper.PL
    }
    # Ensure the estimate is named "yvar" and the conf limits are "LCL" and "UCL"
    nm = names(emms)
    tgts = c(attr(emms, "estName"), attr(emms, "clNames"))
    subs = c("yvar", "LCL", "UCL")
    for (i in 1:3)
        names(emms)[nm == tgts[i]] = subs[i] 
    attr(emms, "estName") = "yvar"
    
    
    # Set up trace vars and key
    tvars = specs$lhs
    if (one.trace <- (length(tvars) == 0)) {
        tlab = ""
        tvars = ".single."
        emms$.single. = 1
        my.key = function(tvars) list()
    }
    else {
        tlab = paste(tvars, collapse = ":")
        my.key = function(tvars) 
            list(space="right", 
                 title = paste(tvars, collapse=" * "), 
                 points = TRUE, 
                 lines=length(lty) > 1,
                 cex.title=1)
    }
    tv = do.call(paste, unname(emms[tvars]))
    emms$tvar = factor(tv, levels=unique(tv))
    
    xvars = specs$rhs
    xv = do.call(paste, unname(emms[xvars]))
    emms$xvar = factor(xv, levels = unique(xv))
    emms = emms[order(emms$xvar), ]
    plotform = yvar ~ xvar
    
    byvars = specs$by
    if (length(byvars) > 0) {
         plotform = as.formula(paste("yvar ~ xvar |", paste(byvars, collapse="*")))
    }
    xlab = ifelse(is.null(xargs$xlab),
                  paste("Levels of", paste(xvars, collapse=" * ")), xargs$xlab)
    rspLbl = paste("Predicted", 
                   ifelse(is.null(emmo@misc$inv.lbl), "response", emmo@misc$inv.lbl))
    ylab = ifelse(is.null(xargs$ylab),
                  ifelse(type == "response", rspLbl, "Linear prediction"),
                  xargs$ylab)
    
    # remove the unneeded stuff from xlabs
    xargs = xargs[setdiff(names(xargs), c("xlab","ylab"))]
    
    if (!plotit) {
        emms$.single. = NULL   # in case we have that trick column
        attr(emms, "labs") = list(xlab = xlab, ylab = ylab, tlab = tlab)
        return (emms)
    }
    
    if (engine == "lattice") {
        # The strips the way I want them
        my.strip = lattice::strip.custom(strip.names = c(TRUE,TRUE), strip.levels = c(TRUE,TRUE), sep = " = ")
        
        TP = TP.orig = lattice::trellis.par.get()
        TP$superpose.symbol$pch = pch
        TP$superpose.line$lty = lty
        if (!is.null(col)) TP$superpose.symbol$col = TP$superpose.line$col = col
        lattice::trellis.par.set(TP)
        
        plotspecs = list(x = plotform, data = emms, groups = ~ tvar, 
                         xlab = xlab, ylab = ylab,
                         strip = my.strip, auto.key = my.key(tvars), type=c("p","l"))
        grobj = do.call(lattice::xyplot, c(plotspecs, xargs))
        lattice::trellis.par.set(TP.orig)
    }  # --- end lattice method
    else {  # engine = "ggplot"
        pos = ggplot2::position_dodge(width = ifelse(CIs|PIs, .1, 0)) # use dodging if CIs
        if (!one.trace) {
            grobj = ggplot2::ggplot(emms, ggplot2::aes_(x = ~xvar, y = ~yvar, color = ~tvar)) +
                ggplot2::geom_point(position = pos) +
                ggplot2::geom_line(ggplot2::aes_(group = ~tvar), position = pos) +
                ggplot2::labs(x = xlab, y = ylab, color = tlab)
        }
        else { # just one trace per plot
            grobj = ggplot2::ggplot(emms, ggplot2::aes_(x = ~xvar, y = ~yvar)) +
                ggplot2::geom_point() +
                ggplot2::geom_line(ggplot2::aes_(group = ~tvar)) +
                ggplot2::labs(x = xlab, y = ylab)
            
        }
        if (PIs) # using linerange w/ extra width and semi-transparent
            grobj = grobj + ggplot2::geom_linerange(ggplot2::aes_(ymin = ~LPL, ymax = ~UPL), 
                                                    position = pos, lwd = 1.25, alpha = .33)
        if (CIs) # using linerange w/ extra width and semi-transparent
            grobj = grobj + ggplot2::geom_linerange(ggplot2::aes_(ymin = ~LCL, ymax = ~UCL), 
                                                    position = pos, lwd = 2, alpha = .5)
        if (length(byvars) > 0) {  # we have by variables 
            if (length(byvars) > 1) {
                byform = as.formula(paste(byvars[1], " ~ ", paste(byvars[-1], collapse="*")))
                grobj = grobj + ggplot2::facet_grid(byform, labeller = "label_both")
            }
            else
                grobj = grobj + ggplot2::facet_wrap(byvars, labeller = "label_both")
        }
    }
    
    grobj
}
