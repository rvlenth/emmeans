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
#' P values smaller than 0.0004 are altered and plotted in a way that makes 
#'   them more distinguishable from one another.
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
#' @param plim numeric vector of value(s) between 0 and 1. These are included
#'   among the observed p values so that the range of tick marks includes at
#'   least the range of \code{plim}. Choosing \code{plim = c(0,1)} will ensure
#'   the widest possible range.
#' @param add.space Numeric value to adjust amount of space used for value labels. Positioning
#'                  of value labels is tricky, and depends on how many panels and the
#'                  physical size of the plotting region. This parameter allows the user to
#'                  adjust the position. Changing it by one unit should shift the position by
#'                  about one character width (right if positive, left if negative).
#' @param ... Additional arguments passed to \code{contrast} and \code{\link{summary.emmGrid}}
#' 
#' 
#' @note The \pkg{ggplot2} and \pkg{scales} packages must be installed in order 
#'   for \code{pwpp} to work.
#'
#' @seealso A numerical display of essentially the same results is available
#'   from \code{\link{pwpm}}
#' @export
#' @examples
#' pigs.lm <- lm(log(conc) ~ source * factor(percent), data = pigs)
#' emm = emmeans(pigs.lm, ~ percent | source)
#' pwpp(emm)
#' pwpp(emm, method = "trt.vs.ctrl1", type = "response", side = ">")
pwpp = function(emm, method = "pairwise", by, sort = TRUE, values = TRUE, rows = ".",
                xlab, ylab, xsub = "", plim = numeric(0), add.space = 0, ...) {
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
    
    sep = get_emm_option("sep")
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
    pf = do.call(paste, c(unname(emm.summ[primv]), sep = sep))
    pemm = suppressMessages(emmeans(emm, primv))
    levs = do.call(paste, c(unname(pemm@grid[primv]), sep = sep))
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
        
        # find ranges to ensure we get tick marks:
        exmaj = c(0, .pvmaj.brk)
        pvtmp = c(plim, con.summ$p.value)
        pvtmp = pvtmp[!is.na(pvtmp)]
        tick.min = max(exmaj[exmaj <= min(pvtmp)])
        tick.max = min(exmaj[exmaj >= max(pvtmp)])
        
        grobj = ggplot2::ggplot(data = con.summ, 
                                ggplot2::aes_(x = ~p.value, y = ~plus,
                                              color = ~minus, group = ~minus)) +
            ggplot2::geom_point(size = 2) +
            ggplot2::geom_segment(ggplot2::aes_(xend = ~p.value, yend = ~midpt)) +
            ggplot2::geom_point(ggplot2::aes(x = tick.min, y = 1), alpha = 0) +
            ggplot2::geom_point(ggplot2::aes(x = tick.max, y = 1), alpha = 0)
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
            tminp = .pval.tran(min(c(tick.min, con.summ$p.value), na.rm=TRUE))
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

.pvmaj.brk = c(.001, .01, .05, .1, .2, .5, 1)
#.pvmin.brk = c(.0001, .0005, seq(.001, .009, by = .001), seq(.02 ,.09, by = .01), .2, .3, .4, .6, .7, .8, .9)
.pvmin.brk = c(.0005, .001, .005, seq(.01 ,.1, by = .01), .2, .3, .4, .5, 1)

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
    savex = x
    x = x[!is.na(x)]
    if ((length(x) <= 3) || (diff(range(x)) == 0))
        return(savex)
    
    kink = function(xx) { # linear spline basis; call with knot subtracted
        xx[xx < 0] = 0
        xx
    }
    
    ### x[x < .00009] = .00009   # forces granulation of extremely small P
    # spread-out the p values less than .0004
    rnk = rank(x)
    sm = which(x < .0004)
    if(length(sm) > 0)
        x[sm] = .0004 - .0003 * rnk[sm] / length(sm)
    
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
    savex[!is.na(savex)] = x
    savex
}



#' Pairwise P-value matrix (plus other statistics)
#'
#' This function presents results from \code{emmeans} and pairwise comparisons
#' thereof in a compact way. It displays a matrix (or matrices) of estimates,
#' pairwise differences, and P values. The user may opt to exclude any of these
#' via arguments \code{means}, \code{diffs}, and \code{pvals}, respectively.
#' To control the direction of the pairwise differences, use \code{reverse};
#' and to control what appears in the upper and lower triangle(s), use \code{flip}.
#' Optional arguments are passed to \code{contrast.emmGrid} and/or 
#' \code{summary.emmGrid}, making it possible to control what estimates 
#' and tests are displayed.
#'
#' @param emm An \code{emmGrid} object
#' @param by Character vector of variable(s) in the grid to condition on. These
#'   will create different matrices, one for each level or level-combination.
#'   If missing, \code{by} is set to \code{emm@misc$by.vars}.
#'   Grid factors not in \code{by} are the \emph{primary} factors:
#'   whose levels or level combinations are compared pairwise.
#' @param reverse Logical value passed to \code{\link{pairs.emmGrid}}.
#'   Thus, \code{FALSE} specifies \code{"pairwise"} comparisons 
#'   (earlier vs. later), and \code{TRUE} specifies \code{"revpairwise"}
#'   comparisons (later vs. earlier).
#' @param pvals Logical value. If \code{TRUE}, the pairwise differences 
#'   of the EMMs are included in each matrix according to \code{flip}.
#' @param means Logical value. If \code{TRUE}, the estimated marginal means
#'   (EMMs) from \code{emm} are included in the matrix diagonal(s).
#' @param diffs Logical value. If \code{TRUE}, the pairwise differences 
#'   of the EMMs are included in each matrix according to \code{flip}.
#' @param flip Logical value that determines where P values and differences 
#'   are placed. \code{FALSE} places the P values in the upper triangle
#'   and differences in the lower, and \code{TRUE} does just the opposite.
#' @param digits Integer. Number of digits to display. If missing,
#'   an optimal number of digits is determined.
#' @param ... Additional arguments passed to \code{\link{contrast.emmGrid}} and 
#'   \code{\link{summary.emmGrid}}. You should \emph{not} include \code{method}
#'   here, because pairwise comparisons are always used. 
#'
#' @return A matrix or `list` of matrices, one for each `by` level.
#' 
#' @seealso A graphical display of essentially the same results is available
#'   from \code{\link{pwpp}}
#' @export
#'
#' @examples
#' warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
#' warp.emm <- emmeans(warp.lm, ~ tension | wool)
#' 
#' pwpm(warp.emm)
#' 
#' # use dot options to specify noninferiority tests
#' pwpm(warp.emm, by = NULL, side = ">", delta = 5, adjust = "none")
pwpm = function(emm, by, reverse = FALSE,
                pvals = TRUE, means = TRUE, diffs = TRUE, 
                flip = FALSE, digits, ...) {
    if(missing(by)) 
        by = emm@misc$by.vars
    
    emm = update(emm, by = by)
    pri = paste(emm@misc$pri.vars, collapse = ":")
    mns = confint(emm, ...)
    mns$lbls = do.call(paste, c(unname(mns[attr(mns, "pri.vars")]), sep = get_emm_option("sep")))
    estName = attr(mns, "estName")
    prs = test(pairs(emm, reverse = reverse, ...), ...)
    diffName = attr(prs, "estName")
    null.hyp = "0"
    if (!reverse) 
        trifcn = lower.tri
    else {
        flip = !flip
        trifcn = upper.tri
    }
    
    if (!is.null(prs$null)) {
        null.hyp = as.character(signif(unique(prs$null), digits = 5))
        if (length(null.hyp) > 1)
            null.hyp = "(various values)"
    }
    mby = .find.by.rows(mns, by)
    pby = .find.by.rows(prs, by)

    
    if(opt.dig <- missing(digits)) {
        tmp = mns[[estName]] + mns[["SE"]] * cbind(rep(-2, nrow(mns)), 0, 2)
        digits = max(apply(tmp, 1, .opt.dig))
        opt.dig = TRUE
    }

    result = lapply(seq_along(mby), function(i) {
        if(opt.dig) {
            pv = prs$p.value[pby[[i]]]
            fpv = sprintf("%6.4f", pv) 
            fpv[pv < 0.0001] = "<.0001"
        }
        else
            fpv = format(prs$p.value[pby[[i]]], digits = digits)
        fmn = format(mns[mby[[i]], estName], digits = digits)
        fdiff = format(prs[pby[[i]], diffName], digits = digits) 
        
        lbls = mns$lbls[mby[[i]]]
        n = length(lbls)
        mat = matrix("", nrow = n, ncol = n, dimnames = list(lbls, lbls))
        if(pvals) {
            mat[trifcn(mat)] = fpv
            mat = t(mat)
        }
        if (diffs)
            mat[trifcn(mat)] = fdiff
        if (means)
            diag(mat) = fmn
        else { # trim off empty row and col
            idx = seq_len(n - 1)
            if (pvals && !diffs)
                mat = mat[idx, 1 + idx]
            if (!pvals && diffs)
                mat = mat[1 + idx, idx]
        }
        if (flip) t(mat)
        else mat
    })
    if (reverse) 
        flip = !flip
    names(result) = paste(paste(by, collapse = ", "), "=", names(mby))
    if (length(result) == 1)
        result = result[[1]]
    class(result) = c("pwpm", "list")
    attr(result, "parms") = c(pvals = pvals, diffs = diffs, means = means, pri = pri,
            estName = estName, diffName = diffName, reverse = reverse, flip = flip,
            type = attr(mns, "type"), adjust = attr(prs, "adjust"),
            side = attr(prs, "side"), delta = attr(prs, "delta"),
            null = null.hyp)
    result
}

#' @export
print.pwpm = function(x, ...) {
    parms = attr(x, "parms")
    attr(x, "class") = attr(x, "parms") = NULL
    if ((islist <- !is.matrix(x)))
        entries = seq_along(x)
    else {
        entries = 1
        m = x 
    }

    for (i in entries) {
        if (islist) {
            cat(paste0("\n", names(x)[i], "\n"))
            m = x[[i]]
        }
        if (parms["means"])
            diag(m) = paste0("[", diag(m), "]")
        print(m, quote = FALSE, right = TRUE, na.print = "nonEst")
    }
    
    # print a parm and its name if present unless it's in excl
    # optional subst is NAMED vector where each possibilitty MUST be present
    catparm = function(f, excl = "0", delim = "  ", quote = TRUE, subst) {
        if(!is.na(pf <- parms[f]) && !(pf %in% excl)) {
            if (!missing(subst)) pf = subst[pf]
            if (quote) pf = dQuote(pf)
            cat(paste0(delim, f, " = ", pf))
        }
    }
    cat(paste0("\nRow and column labels: ", parms["pri"], "\n"))
    if (parms["pvals"]) {
        cat(paste0(ifelse(parms["flip"], "Lower", "Upper"), " triangle: P values "))
        catparm("null", quote = FALSE)
        catparm("side", subst = c("-1" = "<", "1" = ">"))
        catparm("delta", quote = FALSE)
        catparm("adjust", "none")
        cat("\n")
    }
    if (parms["means"]) {
        cat(paste0("Diagonal: [Estimates] (", parms["estName"], ") "))
        catparm("type", "link")
        cat("\n")
    }
    if (parms["diffs"]) {
        cat(paste0(ifelse(parms["flip"], "Upper", "Lower"), " triangle: Comparisons (",
                   parms["diffName"], ")   "))
        if (parms["reverse"]) cat("later vs. earlier\n")
        else cat("earlier vs. later\n")
    }
    invisible(x)
}
    

