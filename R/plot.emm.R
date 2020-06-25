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

# S3 plot method for emmGrid objects
# ... are arguments sent to update()

# DELETED: #' @import ggplot2


#' @rdname plot
#' @importFrom graphics plot
#' @method plot emmGrid
#' @export
plot.emmGrid = function(x, y, type, CIs = TRUE, PIs = FALSE, comparisons = FALSE, 
                    colors = c("black", "blue", "blue", "red"),
                    alpha = .05, adjust = "tukey", int.adjust = "none", intervals, frequentist, ...) {
    if(!missing(intervals))
        CIs = intervals
    if(!missing(type))
        object = update(x, predict.type = type, ..., silent = TRUE)
    else
        object = update(x, ..., silent = TRUE)
    ptype = ifelse(is.null(object@misc$predict.type), "lp", object@misc$predict.type)
    # when we want comparisons, we have a transformation, and we want non-link scale, it's mandatory to regrid first:
    if(comparisons && !is.null(object@misc$tran) && 
            !(ptype %in% c("link", "lp", "linear.predictor")))
        object = regrid(object, transform = ptype)
    if (missing(int.adjust)) {
        int.adjust = object@misc$adjust
        if (is.null(int.adjust))
            int.adjust = "none"
    }
    
    summ = summary(object, infer = c(TRUE, FALSE), adjust = int.adjust, frequentist = frequentist, ...)
    if (is.null(attr(summ, "pri.vars"))) { ## new ref_grid - use all factors w/ > 1 level
        pv = names(x@levels)
        len = sapply(x@levels, length)
        if (max(len) > 1)
            pv = pv[len > 1]
        attr(summ, "pri.vars") = pv
    }
    if (PIs) {
        prd = predict(object, interval = "pred", frequentist = frequentist, ...)
        summ$lpl = prd$lower.PL
        summ$upl = prd$upper.PL
    }
    
    estName = attr(summ, "estName")
    extra = NULL
    bayes = 
    if(comparisons) {
        if ((!is.na(object@post.beta[1])) && (missing(frequentist) || !frequentist))
            stop("Comparison intervals are not implemented for Bayesian analyses")
        extra = object
        extra@misc$comp.alpha = alpha
        extra@misc$comp.adjust = adjust
    }
    .plot.srg(x = summ, CIs = CIs, PIs = PIs, colors = colors, extra = extra, ...)
}

# May use in place of plot.emmGrid but no control over level etc.
# extra is a placeholder for comparison-interval stuff



#' Plot an \code{emmGrid} or \code{summary_emm} object
#' 
#' Methods are provided to plot EMMs as side-by-side CIs, and optionally to display 
#'   \dQuote{comparison arrows} for displaying pairwise comparisons.
#'
#' @rdname plot
#' @param x Object of class \code{emmGrid} or \code{summary_emm}
#' @param y (Required but ignored)
#' @param horizontal Logical value specifying whether the intervals should be
#'   plotted horizontally or vertically
#' @param xlab Character label for horizontal axis
#' @param ylab Character label for vertical axis
#' @param layout Numeric value passed to \code{\link[lattice:xyplot]{dotplot}}
#' @param type Character value specifying the type of prediction desired
#'   (matching \code{"linear.predictor"}, \code{"link"}, or \code{"response"}).
#'   See details under \code{\link{summary.emmGrid}}.
#' @param CIs Logical value. If \code{TRUE}, confidence intervals are
#'   plotted for each estimate.
#' @param PIs Logical value. If \code{TRUE}, prediction intervals are
#'   plotted for each estimate. If \code{objecct} is a Bayesian model,
#'   this requires \code{frequentist = TRUE} and \code{sigma =} (some value).
#'   Note that the \code{PIs} option is \emph{not} available with
#'   \code{summary_emm} objects -- only for \code{emmGrid} objects.
#'   Also, prediction intervals are not available
#'   with \code{engine = "lattice"}.
#' @param comparisons Logical value. If \code{TRUE}, \dQuote{comparison arrows}
#'   are added to the plot, in such a way that the degree to which arrows
#'   overlap reflects as much as possible the significance of the comparison of
#'   the two estimates. (A warning is issued if this can't be done.)
#' @param colors Character vector of color names to use for estimates, CIs, PIs, 
#'   and comparison arrows, respectively. CIs and PIs are rendered with some
#'   transparency, and colors are recycled if the length is less than four;
#'   so all plot elements are visible even if a single color is specified. 
#' @param alpha The significance level to use in constructing comparison arrows
#' @param adjust Character value: Multiplicity adjustment method for comparison arrows \emph{only}.
#' @param int.adjust Character value: Multiplicity adjustment method for the plotted confidence intervals \emph{only}.
#' @param intervals If specified, it is used to set \code{CIs}. This is the previous
#'   name of \code{CIs} and is provided for backward compatibility.
#' @param frequentist Logical value. If there is a posterior MCMC sample and 
#'   \code{frequentist} is non-missing and TRUE, a frequentist summary is used for
#'   obtaining the plot data, rather than the posterior point estimate and HPD
#'   intervals. This argument is ignored when it is not a Bayesian model.
#' @param plotit Logical value. If \code{TRUE}, a graphical object is returned;
#'   if \code{FALSE}, a data.frame is returned containing all the values
#'   used to construct the plot.
#' @param ... Additional arguments passed to \code{\link{update.emmGrid}}, 
#'   \code{\link{predict.emmGrid}}, or
#'   \code{\link[lattice:xyplot]{dotplot}}
#'
#' @return If \code{plotit = TRUE}, a graphical object is returned.
#' 
#'   If \code{plotit = FALSE}, a \code{data.frame} with the table of
#'   EMMs that would be plotted. In the latter case, the estimate being plotted
#'   is named \code{the.emmean}, and any factors involved have the same names as
#'   in the object. Confidence limits are named \code{lower.CL} and
#'   \code{upper.CL}, prediction limits are named \code{lpl} and \code{upl}, and
#'   comparison-arrow limits are named \code{lcmpl} and \code{ucmpl}.
#'   There is also a variable named \code{pri.fac} which contains the factor 
#'   combinations that are \emph{not} among the \code{by} variables.

#' @section Details:
#' If any \code{by} variables are in force, the plot is divided into separate
#' panels. These functions use the \code{\link[lattice:xyplot]{dotplot}} function, and
#' thus require that the \pkg{lattice} package be installed. For
#' \code{"summary_emm"} objects, the \code{\dots} arguments in \code{plot}
#' are passed \emph{only} to \code{dotplot}, whereas for \code{"emmGrid"}
#' objects, the object is updated using \code{\dots} before summarizing and
#' plotting.
#' 
#' In plots with \code{comparisons = TRUE}, the resulting arrows are only
#' approximate, and in some cases may fail to accurately reflect the pairwise
#' comparisons of the estimates -- especially when estimates having large and
#' small standard errors are intermingled in just the wrong way. Note that the
#' maximum and minimum estimates have arrows only in one direction, since there
#' is no need to compare them with anything higher or lower, respectively. See
#' the \href{../doc/xplanations.html#arrows}{\code{vignette("xplanations",
#' "emmeans")}} for details on how these are derived.
#' 
#' If \code{adjust} or \code{int.adjust} are not supplied, they default to the 
#' internal \code{adjust} setting saved in \code{pairs(x)} and \code{x} 
#' respectively (see \code{\link{update.emmGrid}}).
#' 
#' @importFrom graphics plot
#' @method plot summary_emm
#' @export
#'
#' @examples
#' warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
#' warp.emm <- emmeans(warp.lm, ~ tension | wool)
#' plot(warp.emm)
#' plot(warp.emm, by = NULL, comparisons = TRUE, adjust = "mvt", 
#'      horizontal = FALSE, colors = "darkgreen")
plot.summary_emm = function(x, y, horizontal = TRUE, CIs = TRUE,
                            xlab, ylab, layout, 
                            colors = c("black", "blue", "blue", "red"), intervals, 
                            plotit = TRUE, ...) {
    if(!missing(intervals))
        CIs = intervals
    .plot.srg (x, y, horizontal, xlab, ylab, layout, CIs = CIs, colors = colors, plotit = plotit, ...)
}

# Workhorse for plot.summary_emm
.plot.srg = function(x, y, 
                     horizontal = TRUE, xlab, ylab, layout, colors,
                     engine = get_emm_option("graphics.engine"),
                     CIs = TRUE, PIs = FALSE, extra = NULL, 
                     plotit = TRUE, ...) {
    
    engine = match.arg(engine, c("ggplot", "lattice"))
    if (engine == "ggplot") 
        .requireNS("ggplot2", 
                   "The 'ggplot' engine requires the 'ggplot2' package be installed.")
    if (engine == "lattice")
        .requireNS("lattice", 
                   "The 'lattice' engine requires the 'lattice' package be installed.")
    
    summ = x # so I don't get confused
    estName = "the.emmean"
    names(summ)[which(names(summ) == attr(summ, "estName"))] = estName
    clNames = attr(summ, "clNames")
    if (is.null(clNames)) {
        warning("No information available to display confidence limits")
        lcl = ucl = summ[[estName]]
    }
    else {
        lcl = summ[[clNames[1]]]
        ucl = summ[[clNames[2]]]
    }
    
    if (engine == "lattice") { # ---------- lattice-specific stuff ----------
        
        # Panel functions...
        prepanel.ci = function(x, y, horizontal=TRUE, CIs=TRUE,
                               lcl, ucl, subscripts, ...) {
            x = as.numeric(x)
            lcl = as.numeric(lcl[subscripts])
            ucl = as.numeric(ucl[subscripts])
            if (!CIs) # no special scaling needed
                list()
            else if (horizontal)
                list(xlim = range(x, ucl, lcl, finite = TRUE)) 
            else
                list(ylim = range(y, ucl, lcl, finite = TRUE)) 
        }
        panel.ci <- function(x, y, horizontal=TRUE, CIs=TRUE,
                             lcl, ucl, lcmpl, rcmpl,                          subscripts, pch = 16, 
                             lty = dot.line$lty, lwd = dot.line$lwd, 
                             col = dot.symbol$col, col.line = dot.line$col, ...) {
            dot.line <- lattice::trellis.par.get("dot.line")
            dot.symbol <- lattice::trellis.par.get("dot.symbol")
            x = as.numeric(x)
            y = as.numeric(y)
            lcl = as.numeric(lcl[subscripts])
            ucl = as.numeric(ucl[subscripts])
            compare = !is.null(lcmpl)
            if(compare) {
                lcmpl = as.numeric(lcmpl[subscripts])
                rcmpl = as.numeric(rcmpl[subscripts])
            }
            if(horizontal) {
                lattice::panel.abline(h = unique(y), col = col.line, lty = lty, lwd = lwd)
                if(CIs) 
                    lattice::panel.arrows(lcl, y, ucl, y, col = col, length = .6, unit = "char", angle = 90, code = 3)
                if(compare) {
                    s = (x > min(x))
                    lattice::panel.arrows(lcmpl[s], y[s], x[s], y[s], length = .5, unit = "char", code = 1, col = "red", type = "closed", fill="red")
                    s = (x < max(x))
                    lattice::panel.arrows(rcmpl[s], y[s], x[s], y[s], length = .5, unit = "char", code = 1, col = "red", type = "closed", fill="red")
                }
            }
            else {
                lattice::panel.abline(v = unique(x), col = col.line, lty = lty, lwd = lwd)
                if(CIs)
                    lattice::panel.arrows(x, lcl, x, ucl, col=col, length = .6, unit = "char", angle = 90, code = 3)
                if(compare) {
                    s = (y > min(y))
                    lattice::panel.arrows(x[s], lcmpl[s], x[s], y[s], length = .5, unit = "char", code = 1, col = "red", type = "closed", fill="red")
                    s = (y < max(y))
                    lattice::panel.arrows(x[s], rcmpl[s], x[s], y[s], length = .5, unit = "char", code = 1, col = "red", type = "closed", fill="red")
                }
            }
            lattice::panel.xyplot(x, y, pch=16, ...)
        }
        my.strip = lattice::strip.custom(strip.names = c(TRUE,TRUE), strip.levels = c(TRUE,TRUE), sep = " = ")
        
    } # ---------- end lattice-specific -----------
    
    sep = get_emm_option("sep")
    priv = attr(summ, "pri.vars")
    pf = do.call(paste, c(unname(summ[priv]), sep = sep))
    if (length(pf) == 0)
        pf = "1"
    summ$pri.fac = factor(pf, levels=unique(pf))
    chform = ifelse(horizontal,
                    paste("pri.fac ~", estName),
                    paste(estName, "~ pri.fac"))
    
    byv = attr(summ, "by.vars")
    if (!is.null(byv) && length(byv) > 0) {
        chform = paste(chform, "|", paste(byv, collapse="*"))
        lbv = do.call("paste", c(unname(summ[byv]), sep = sep)) # strings for matching by variables
        ubv = unique(lbv)
    }
    else {
        lbv = rep(1, nrow(summ))
        ubv = 1
    }
    
    
    # Obtain comparison limits
    if (!is.null(extra)) {
        # we need to work on the linear predictor scale
        # typeid = 1 -> response, 2 -> other
        typeid = pmatch(extra@misc$predict.type, "response", nomatch = 2)
        if(length(typeid) < 1) typeid = 2        
        if (typeid == 1)
            est = predict(extra, type = "lp")
        else
            est = summ[[estName]]
        
        alpha = extra@misc$comp.alpha
        adjust = extra@misc$comp.adjust
        psumm = confint(pairs(extra), level = 1 - alpha, type = "lp", adjust = adjust)
        k = ncol(psumm)
        del = (psumm[[k]] - psumm[[k-1]]) / 4 # half the halfwidth, on lp scale
        diff = psumm[[attr(psumm, "estName")]]
        overlap = apply(psumm[ ,(k-1):k], 1, function(x) 2*min(-x[1],x[2])/(x[2]-x[1]))
        
        # figure out by variables and indexes (lbv, ubv already defined)
        if(is.null(byv))
            pbv = rep(1, nrow(psumm))
        else
            pbv = do.call("paste", c(unname(psumm[byv]), sep = sep))
        neach = length(lbv) / length(ubv)
        # indexes for pairs results -- est[id1] - est[id2]
        id1 = rep(seq_len(neach-1), rev(seq_len(neach-1)))
        id2 = unlist(sapply(seq_len(neach-1), function(x) x + seq_len(neach-x)))
        # list of psumm row numbers involved in each summ row
        involved = lapply(seq_len(neach), function(x) union(which(id2==x), which(id1==x)))
        
        # initialize arrays
        mind = numeric(length(lbv))   # for minima of del
        llen = rlen = numeric(neach)  # for left and right arrow lengths
        npairs = length(id1)
        iden = diag(rep(1, 2*neach))
        
        for (by in ubv) {
            d = del[pbv == by]
            rows = which(lbv == by)
            for(i in seq_len(neach)) 
                mind[rows[i]] = min(d[involved[[i]]])
            
            # Set up regression equations to match arrow overlaps with interval overlaps
            # We'll add rows later (with weights 1) to match with mind values
            lmat = rmat = matrix(0, nrow = npairs, ncol = neach)
            y = numeric(npairs)
            v1 = 1 - overlap[pbv == by]
            dif = diff[pbv == by]
            for (i in seq_len(npairs)) {
                #wgt = 6 * max(0, ifelse(v1[i] < 1, v1[i], 2-v1[i]))
                wgt = 3 + 20 * max(0, .5 - (1 - v1[i])^2)
                # really this is sqrt of weight
                if (dif[i] > 0)   # id2  <----->  id1
                    lmat[i, id1[i]] = rmat[i, id2[i]] = wgt*v1[i]
                else  # id1  <----->  id2
                    rmat[i, id1[i]] = lmat[i, id2[i]] = wgt*v1[i]
                y[i] = wgt * abs(dif[i])
            }
            X = rbind(cbind(lmat, rmat),iden)
            y = c(y, rep(mind[rows], 2))
            soln = qr.coef(qr(X), y)
            ll = llen[rows] = soln[seq_len(neach)]
            rl = rlen[rows] = soln[neach + seq_len(neach)]
            
            # Abort if negative lengths
            if (any(c(rl, ll) < 0)) {
                stop("Aborted -- Some comparison arrows have negative length!", 
                     call. = FALSE)
            }
            
            # Overlap check
            for (i in seq_len(npairs)) {
                v = 1 - v1[i]
                obsv = 1 - abs(dif[i]) / ifelse(dif[i] > 0, 
                                                ll[id1[i]] + rl[id2[i]], 
                                                rl[id1[i]] + ll[id2[i]])
                if (v*obsv < 0)
                    warning("Comparison discrepancy in group ", by, 
                            ", ", psumm[i, 1], 
                            ":\n    Target overlap = ", round(v, 4),
                            ", overlap on graph = ", round(obsv, 4),
                            call. = FALSE)
            }
        }
        # shorten arrows that go past the data range
        rng = range(est)
        diffr = diff(rng)
        ii = which(est - llen < rng[1])
        llen[ii] = est[ii] - rng[1] + .02 * diffr
        ii = which(est + rlen > rng[2])
        rlen[ii] = rng[2] - est[ii] + .02 * diffr
        
        # remove arrows completely from extremes
        llen[est < rng[1] + .0001 * diffr] = NA
        rlen[est > rng[2] - .0001 * diffr] = NA
        
        invtran = I
        if (typeid == 1) {
            tran = extra@misc$tran
            if(is.character(tran)) {
                link = try(make.link(tran), silent=TRUE)
                if (!inherits(link, "try-error"))
                    invtran = link$linkinv
            }
            else if (is.list(tran))
                invtran = tran$linkinv
        }
        
        lcmpl = summ$lcmpl = as.numeric(invtran(est - llen))
        rcmpl = summ$rcmpl = as.numeric(invtran(est + rlen))
    }
    else lcmpl = rcmpl = NULL
    
    if(!plotit) 
        return(as.data.frame(summ))
    
    
    facName = paste(priv, collapse=":")
    
    if(length(colors) < 4)
        colors = rep(colors, 4)
    dot.col = colors[1]
    CI.col = colors[2]
    PI.col = colors[3]
    comp.col = colors[4]
    
    if (engine == "lattice") {
        if (missing(layout)) {
            layout = c(1, length(ubv))
            if(!horizontal) 
                layout = rev(layout)
        }
        
        form = as.formula(chform)
        if (horizontal) {
            if (missing(xlab)) xlab = attr(summ, "estName")
            if (missing(ylab)) ylab = facName
            lattice::dotplot(form, prepanel=prepanel.ci, panel=panel.ci, 
                             strip = my.strip, horizontal = TRUE,
                             ylab = ylab, xlab = xlab,
                             data = summ, CIs = CIs, lcl=lcl, ucl=ucl, 
                             lcmpl=lcmpl, rcmpl=rcmpl, layout = layout, ...)
        }
        else {
            if (missing(xlab)) xlab = facName
            if (missing(ylab)) ylab = attr(summ, "estName")
            lattice::dotplot(form, prepanel=prepanel.ci, panel=panel.ci, 
                             strip = my.strip, horizontal = FALSE,
                             xlab = xlab, ylab = ylab,
                             data = summ, CIs = CIs, lcl=lcl, ucl=ucl, 
                             lcmpl=lcmpl, rcmpl=rcmpl, layout = layout, ...)
        }
    } # --- end lattice plot
    else {  ## ggplot method
        summ$lcl = lcl
        summ$ucl = ucl

        # construct horizontal plot - flip coords later if necessary
        grobj = ggplot2::ggplot(summ, ggplot2::aes_(x = ~the.emmean, y = ~pri.fac)) 
        if (PIs) 
            grobj = grobj + ggplot2::geom_segment(ggplot2::aes_(x = ~lpl, xend = ~upl, 
                                y = ~pri.fac, yend = ~pri.fac), 
                                color = PI.col, lwd = 2.5, alpha = .15)
        if (CIs) 
            grobj = grobj + ggplot2::geom_segment(ggplot2::aes_(x = ~lcl, xend = ~ucl, 
                                y = ~pri.fac, yend = ~pri.fac), 
                                color = CI.col, lwd = 4, alpha = .25)
        if (!is.null(extra)) {
            grobj = grobj + 
                ggplot2::geom_segment(ggplot2::aes_(x = ~the.emmean, 
                        xend = ~lcmpl, y = ~pri.fac, yend = ~pri.fac), 
                        arrow = ggplot2::arrow(length = ggplot2::unit(.07, "inches"), 
                            type = "closed"), color = comp.col, 
                        data = summ[!is.na(summ$lcmpl), ]) +
                ggplot2::geom_segment(ggplot2::aes_(x = ~the.emmean, 
                    xend = ~rcmpl, y = ~pri.fac, yend = ~pri.fac), 
                    arrow = ggplot2::arrow(length = ggplot2::unit(.07, "inches"), 
                            type = "closed"), color = comp.col,
                    data = summ[!is.na(summ$rcmpl), ])
        }
        if (length(byv) > 0)
            grobj = grobj + ggplot2::facet_grid(as.formula(paste(paste(byv, collapse = "+"), " ~ .")), 
                                       labeller = "label_both")
        grobj = grobj + ggplot2::geom_point(color = dot.col, size = 2)
        if (missing(xlab)) xlab = attr(summ, "estName")
        if (missing(ylab)) ylab = facName
            
        if(!horizontal)
            grobj = grobj + ggplot2::coord_flip()
        
        grobj + ggplot2::labs(x = xlab, y = ylab)
    }
}
