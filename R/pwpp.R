##############################################################################
#    Copyright (c) 2012-2019 Russell V. Lenth                                #
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

#' Pairwise P-value plot
#' 
#' Constructs a plot of P values associated with pairwise comparisons of 
#' estimated marginal means.
#' 
#' Factor levels (or combinations thereof) are plotted on the vertical scale, and P values
#' are plotted on the horizontal scale. Each P value is plotted twice -- at
#' vertical positions corresponding to the levels being compared -- and connected by
#' a line segment. Thus, it is easy to visualize which P values are small and large,
#' and which levels are compared. In addition, factor levels are color-coded, and the points
#' and half-line segments appear in the color of the other level.
#' The P-value scale is nonlinear, so as to stretch-out smaller P values and
#' compress larger ones.
#' 
#' If \code{xlab}, \code{ylab}, and \code{xsub} are not provided, reasonable labels
#' are created. \code{xsub} is used to note special features; e.g., equivalence
#' thresholds or one-sided tests.
#' 
#' @param emm An \code{emmGrid} object
#' @param method Character or list. Passed to \code{\link{contrast}}, and defines 
#'           the contrasts to be displayed. Any contrast method may be used,
#'           provided that each contrast includes one coefficient of \code{1},
#'           one coefficient of \code{-1}, and the rest \code{0}. That is, calling
#'           \code{contrast(object, method)} produces a set of comparisons, each with
#'           one estimate minus another estimate.
#' @param by Character vector of variable(s) in the grid to condition on. These will
#'           create different panels, one for each level or level-combination.
#'           Grid factors not in \code{by} are the \emph{primary} factors: 
#'           whose levels or level combinations are compared pairwise.
#' @param sort Logical value. If \code{TRUE}, levels of the factor combinations are
#'             ordered by their marginal means. If \code{FALSE}, they appear in
#'             order based on the existing ordering of the factor levels involved.
#'             Note that the levels are ordered the same way in all panels, and in
#'             many cases this implies that the means in any particular panel
#'             will \emph{not} be ordered even when \code{sort = TRUE}.
#' @param values Logical value. If \code{TRUE}, the values of the EMMs are included
#'               in the plot. When there are several side-by-side panels due
#'               to \code{by} variable(s), the labels showing values start
#'               stealing a lot of space from the plotting area; in those cases,
#'               it may be desirable to specify \code{FALSE} or use \code{rows}
#'               so that some panels are vertically stacked.
#' @param rows Character vector of which \code{by} variable(s) are used to define
#'           rows of the panel layout. Those variables in \code{by} not included in 
#'           \code{rows} define columns in the array of panels.
#'           A \code{"."} indicates that only one row
#'           is used, so all panels are stacked side-by-side.
#' @param xlab Character label to use in place of the default for the P-value axis.
#' @param ylab Character label to use in place of the default for the primary-factor axis.
#' @param xsub Character label used as caption at the lower right of the plot.
#' @param add.space Numeric value to adjust amount of space used for value labels. Positioning
#'                  of value labels is tricky, and depends on how many panels and the
#'                  physical size of the plotting region. This parameter allows the user to
#'                  adjust the position. Changing it by one unit should shift the position by
#'                  about one character width (right if positive, left if negative).
#' @param ... Additional arguments passed to \code{contrast} and \code{\link{summary.emmGrid}}
#' 
#' @note The \pkg{ggplot2} package must be installed in order for \code{pwpp} to work.
#'
#' @export
#' @examples
#' pigs.lm <- lm(log(conc) ~ source * factor(percent), data = pigs)
#' emm = emmeans(pigs.lm, ~ percent | source)
#' pwpp(emm)
#' pwpp(emm, method = "trt.vs.ctrl1", type = "response", side = ">")
pwpp = function(emm, method = "pairwise", by, sort = TRUE, values = TRUE, rows = ".",
                xlab, ylab, xsub = "", add.space = 0, ...) {
    if(missing(by)) 
        by = emm@misc$by.vars
    
    if(rows != "." && !(rows %in% by))
        stop("'rows' must be a subset of the 'by' variables")
    
    args = list(object = emm, method = method, by = by, ...)
    args$interaction = args$simple = args$offset = NULL
    con = do.call(contrast, args)
    
    args = list(object = emm, infer = c(FALSE, FALSE), by = by, ...)
    emm.summ = do.call(summary.emmGrid, args)
    
    args = list(object = con, infer = c(FALSE, TRUE), ...)
    args$null = NULL
    con.summ = do.call(summary.emmGrid, args)
    
    if(missing(xlab)) {
        adjust = .cap(attr(con.summ, "adjust"))
        delta = attr(con.summ, "delta")
        side = attr(con.summ, "side")
        xlab = "P value"
        if (adjust != "None")
            xlab = paste0(adjust, "-adjusted ", xlab)
        
        if (delta != 0)
            xsub = paste(c("Nonsuperiority", "Equivalence", "Noninferiority")[side + 2],
                         "test with threshold", delta)
        else
            xsub = c("Left-sided tests", "", "Right-sided tests")[side + 2]
    }
    if(missing(ylab))
        ylab = paste(attr(emm.summ, "pri.vars"), collapse = ":")
    
    # figure out levels being compared
    cf = coef(con)
    use = setdiff(names(cf), names(con@misc$orig.grid))
    idx = apply(as.matrix(cf[use]), 2, function(x) {
        if(!all(range(x) == c(-1,1)) || (sum(abs(x)) != 2))
            stop("Each contrast must be a comparison of two estimates")
        c(which(x == 1), which(x == -1))
    })
    primv = setdiff(attr(emm.summ, "pri.vars"), by)
    pf = do.call(paste, c(unname(emm.summ[primv]), sep = ":"))
    pemm = suppressMessages(emmeans(emm, primv))
    levs = do.call(paste, c(unname(pemm@grid[primv]), sep = ":"))
    if(sort) ord = order(predict(pemm))
    else ord = seq_along(pf)
    pf = emm.summ$pri.fac = factor(pf, levels = levs[ord])
    
    con.summ$plus = pf[idx[1, ]]
    con.summ$minus = pf[idx[2, ]]
    
    estName = attr(emm.summ, "estName")
    
    
########## The rest should probably be done in a separate function ################
    
    .requireNS("ggplot2", 
               "pwpp requires the 'ggplot2' package be installed.", call. = FALSE)
    
        # granulate values in each group so they won't overlap
        # do this on the transformed (plotted) scale
        byr = .find.by.rows(con.summ, by)
        for (r in byr) {
            pv = con.summ$p.value[r]
            con.summ$p.value[r] = gran(pv)
        }
        
        # form the reverse half & get midpoints
        con.summ$midpt = (as.numeric(con.summ$plus) + as.numeric(con.summ$minus)) / 2
        tmp = con.summ
        tmp$plus = con.summ$minus
        tmp$minus = con.summ$plus
        con.summ = rbind(con.summ, tmp)
        
        grobj = ggplot2::ggplot(data = con.summ, 
                                ggplot2::aes_(x = ~p.value, y = ~plus,
                                              color = ~minus, group = ~minus)) +
            ggplot2::geom_point(size = 2) +
            ggplot2::geom_segment(ggplot2::aes_(xend = ~p.value, yend = ~midpt))
        if (!is.null(by)) {
            cols = setdiff(by, rows)
            if (length(cols) > 0)
                ncols = length(unique(do.call(paste, unname(con.summ[cols]))))
            else {
                ncols = 1
                cols = "."
            }
            grobj = grobj + ggplot2::facet_grid(
                as.formula(paste( paste(rows, collapse = "+"), "~", 
                                  paste(cols, collapse = "+"))), 
                labeller = "label_both")
        }
        else
            ncols = 1
        if (values) {
            emm.summ$minus = emm.summ$pri.fac # for consistency in grouping/labeling
            dig = .opt.dig(emm.summ[, estName])
            emm.summ$fmtval = format(emm.summ[[estName]], digits = dig)
            tminp = .pval.tran(min(con.summ$p.value))
            pos = .pval.inv(tminp - .025)
            lpad = .012 * (add.space + max(nchar(emm.summ$fmtval))) * ncols # how much space needed for labels rel to (0,1)
            lpad = lpad * (1.1 - tminp) # scale closer to actual width of scales
            lpos = .pval.inv(tminp - lpad)   # pvalue at left end of label
            grobj = grobj + 
                ggplot2::geom_label(data = emm.summ, ggplot2::aes_(x = pos, y = ~minus,
                                                                   label = ~fmtval, hjust = "right"), size = 2.8) +
                ggplot2::geom_point(ggplot2::aes_(x = lpos, y = 1), alpha = 0) # invisible point to stake out space
        }
        else
            lpad = 0
        
        .requireNS("scales", 
                   "pwpp requires the 'scales' package be installed.", call. = FALSE)
        .pvtrans = scales::trans_new("Scaled P value", 
                                     transform = function(x) .pval.tran(x), 
                                     inverse = function(p) .pval.inv(p),
                                     format = function(x) format(x, drop0trailing = TRUE, scientific = FALSE),
                                     domain = c(0,1) )
        grobj = grobj + ggplot2::scale_x_continuous(trans = .pvtrans, 
                                                    breaks = .pvmaj.brk, minor_breaks = .pvmin.brk) +
            #### I plotted an extra point instead of expanding scale
            #### expand = ggplot2::expand_scale(add = c(.025 + lpad, .025))) +
            ggplot2::guides(color = "none")
        
        grobj + ggplot2::labs(x = xlab, y = ylab, caption = xsub)
}
    
# capitalize
.cap = function(s)
    paste0(toupper(substring(s, 1, 1)), substring(s, 2))

          
### Scale-transformation code: We stretch out small P values without stretching 
### extremely small ones too much -- via a combination of log and normal cdf functions

.tran.ctr = -2.5  # params of normal cdf transf of log(.value)
.tran.sd = 3
.tran.div = pnorm(0, .tran.ctr, .tran.sd)

.pvmaj.brk = c(.001, .01, .05, .1, .5, 1)
#.pvmin.brk = c(.0001, .0005, seq(.001, .009, by = .001), seq(.02 ,.09, by = .01), .2, .3, .4, .6, .7, .8, .9)
.pvmin.brk = c(.0001, .0005, .001, .005, seq(.01 ,.1, by = .01), .2, .3, .4, .5, 1)

# transforms x in (0, 1] to t in (0,1], while any x's < 0 or > 1 are preserved as is
# This allows me to position marginal labels etc. where I want
.pval.tran = function(x) {
    xx = sapply(x, function(.) min(max(., .00005), 1))
    rtn = pnorm(log(xx), .tran.ctr, .tran.sd) / .tran.div
    spec = which((x < .00005) | (x > 1))
    rtn[spec] = x[spec]
    rtn
}

.pval.inv = function(p) {
    pp = sapply(p, function(.) min(max(., .00005), .99995))
    rtn = exp(qnorm(.tran.div * pp, .tran.ctr, .tran.sd))
    spec = which((p < .00005) | (p > .99995))
    rtn[spec] = p[spec]
    rtn
}

# For scale_x_continuous(): (moved to body of fcn so I can check w/ requyireNamespace)
# .pvtrans = scales::trans_new("Scaled P value", 
#                              transform = function(x) .pval.tran(x), 
#                              inverse = function(p) .pval.inv(p),
#                              format = function(x) format(x, drop0trailing = TRUE, scientific = FALSE),
#                              domain = c(0,1) )

### quick & dirty algorithm to Stretch out P values so more distinguishable on transformed scale
gran = function(x, min_incr = .01) {
    if ((length(x) <= 3) || (diff(range(x)) == 0))
        return(x)
    
    kink = function(xx) { # linear spline basis; call with knot subtracted
        xx[xx < 0] = 0
        xx
    }
    
    x[x < .00009] = .00009   # forces granulation of extremely small P
    
    ord = order(x)
    tx = log(x[ord])
    df = diff(tx)
    incr = sapply(df, function(.) max(., min_incr))
    ttx = cumsum(c(tx[1], incr))
    if (length(x) > 5)
        tx = predict(lm(tx ~ ttx + kink(ttx + 1.5) + kink(ttx + 3) + kink(ttx + 4.5)))
    else
        tx = ttx
    tx[tx > 0] = 0
    x[ord] = exp(tx)
    x
}

