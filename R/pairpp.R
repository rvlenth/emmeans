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

#' Pairwise p-value plot
#' 
#' Constructs a plot of P values associated with pairwise comparisons of 
#' estimated marginal means. 
#' 
#' @param object An \code{emmGrid} object
#' @param sort Logical flag to determine whether means are sorted first
#' @param values Logical value. If \code{TRUE}, the values of the EMMs are included
#'               in the plot. When there are several side-by-side panels due
#'               to \code{by} variable(s), the labels showing values start
#'               stealing a lot of space from the plotting area; in those cases,
#'               it may be desiravle to specify \code{FALSE} or use \code{rows}
#'               so that some panels are vertically stacked.
#' @param by Character vector of variable(s) in the grid to condition on. These will
#'           create different panels, one for each level or level-combination.
#'           Grid factors not in \code{by} are the \emph{primary} factors: 
#'           whose levels or level combinations are compared pairwise.
#' @param rows Character vector of which \code{by} variable(s) are used to define
#'           rows of the panel layout. A \code{"."} indicates that only one row
#'           is used, so all panels are stacked side-by-side.
#' @param adjust Character \code{adjust} argument passed to \code{summary.emmGrid}.
#' @param xlab Character label to use in place of the default for the P-value axis.
#' @param ylab Character label to use in place of the default for the primary-factor axis.
#' @param ... Additional arguments passed to \code{\link{summary.emmGrid}}
#'
#' @export
pwpp = function(object, sort = TRUE, values = TRUE, by, rows = ".", 
                adjust = "Tukey", xlab, ylab, ...) {
    if(missing(by)) 
        by = object@misc$by.vars
    
    if(rows != "." && !(rows %in% by))
        stop("'rows' must be a subset of the 'by' variables")
    
    args = list(...)
    args$object = pairs(object, by = by)
    args$infer = c(FALSE, TRUE)
    args$adjust = adjust
    pvsumm = do.call(summary.emmGrid, args)
    pv.bys = .find.by.rows(pvsumm, by)
    
    emsumm = summary(object, by = by, infer = c(FALSE, FALSE))
    em.bys = .find.by.rows(emsumm, by)
    estName = attr(emsumm, "estName")
    dig = .opt.dig(emsumm[, estName])
    priv = attr(emsumm, "pri.vars")
    pf = do.call(paste, c(unname(emsumm[priv]), sep = ":"))
    emsumm$pri.fac = factor(pf, levels = unique(pf))
    if(missing(xlab))
        xlab = ifelse(tolower(adjust) == "none", "Unadjusted P value",
                      paste0(adjust, "-adjusted P value"))
    if(missing(ylab))
        ylab = paste(attr(emsumm, "pri.vars"), collapse = ":")

    data = list()
    for (i in seq_along(pv.bys)) {
        idx = pv.bys[[i]]
        emi = em.bys[[i]]
        p.value = pvsumm$p.value[idx]
        lvls = emsumm$pri.fac[emi]
        est = emsumm[emi, estName]
        pmat = diag(length(lvls)) * 1.025
        con = pairwise.emmc(lvls)
        for (j in seq_along(con)) {
            i = which(con[[j]] != 0)
            pmat[i[1],i[2]] = pmat[i[2], i[1]] = p.value[j]
        }
        df = expand.grid(level = lvls, versus = lvls)
        for(v in by)
            df[[v]] = pvsumm[idx, ][[v]][1]
        df$fmtval = rep(format(est, digits = dig), length(lvls))
        df$p.value = as.numeric(pmat)
        data = rbind(data, df)
    }
    key.dat = data[data$level == data$versus, ]
    # extra space needed for labels?
    if (!values)
        lpad = 0
    else {
        cols = setdiff(by, rows)
        if (length(cols) > 0)
            ncols = length(unique(do.call(paste, unname(data[cols]))))
        else
            ncols = 1
        lpad = .012 * max(nchar(key.dat$fmtval)) * ncols
    }
    grobj = ggplot2::ggplot(data = data, 
                ggplot2::aes_(x = ~p.value, y = ~level,
                    color = ~versus, group = ~versus)) +
        ggplot2::geom_point(size = 2) +
        ggplot2::geom_point(data = key.dat, size = 3, shape = "square") +
        ggplot2::geom_path(linetype = "dashed") +
        ggplot2::scale_x_continuous(trans = .pvtrans, 
            breaks = .pvmaj.brk, minor_breaks = .pvmin.brk,
            expand = expand_scale(add = c(.025 + lpad, .05))) +
        ggplot2::guides(color = "none")
    if (values) {
        # we have to remember that scale is reversed so labels are at pos > 1
        grobj = grobj + 
             ggplot2::geom_label(data = key.dat, aes_(x = 1.07, y = ~level,
                    label = ~fmtval, hjust = "right"), size = 2.8)
    }
    if (length(pv.bys) > 1) {
        cols = setdiff(by, rows)
        if (length(cols) == 0) cols = "."
        grobj = grobj + ggplot2::facet_grid(
            as.formula(paste( paste(rows, collapse = "+"), "~", 
                              paste(cols, collapse = "+"))), 
            labeller = "label_both")
    }
    grobj + ggplot2::labs(x = xlab, y = ylab)
}

### Scale-transformation code: We stretch out small P values without stretching 
### extremely small ones too much -- via a combination of log and normal cdf functions

.tran.ctr = -2.5  # params of normal cdf transf of log(.value)
.tran.sd = 3
.tran.div = pnorm(0, .tran.ctr, .tran.sd)

.pvmaj.brk = c(.001, .01, .05, .1, .5)
.pvmin.brk = c(seq(.001, .009, by = .001), seq(.02 ,.09, by = .01), .2, .3, .4, .6, .7, .8, .9)

# transforms x in (0, 1] to t in (0,1], while any x's < 0 or > 1 are preserved as is
# This allows me to position marginal labels etc. where I want
.pval.tran = function(x) {
    xx = sapply(x, function(.) min(max(., .0001), 1))
    rtn = pnorm(log(xx), .tran.ctr, .tran.sd) / .tran.div
    rtn[x < 0] = x[x < 0]
    rtn[x > 1] = x[x > 1]
    rtn
}

# Being careful with break points, we should never encounter values outside limits
.pval.inv = function(p) {
    pp = sapply(p, function(.) min(max(., .001), .999))
    rtn = exp(qnorm(.tran.div * pp, .tran.ctr, .tran.sd))
    rtn[p < .001] = 0
    rtn[p > .999] = 1
    rtn
}

# For scale_x_continuous():
.pvtrans = scales::trans_new("Scaled P value", 
                           transform = function(x) 1 - .pval.tran(x), 
                           inverse = function(p) .pval.inv(1 - p),
                           format = function(x) format(x, drop0trailing = TRUE),
                           domain = c(0,1) )

## 1st approx: 
## ggplot(data = dat, aes_(x=~tpv, y=~emm, color = ~minus, group=~minus)) + 
##   geom_point() +
##   geom_point(data = dat.diag, size = 3) +
##   geom_path(linetype = "dashed") + 
##   scale_x_continuous(trans = .trans, breaks = .maj.brk, minor_breaks = .min.brk) +
##   guides(color = "none") 



### function to structurally jitter values so they are all distinct
granulate = function(x, min_incr = .01) {
    ord = order(x)
    x = x[ord]
    df = diff(x)
    incr = sapply(df, function(.) max(., min_incr))
    act = which(df > min_incr)  # these are places we can make adjustments
    csd = cumsum(df)[act]
    sdf = rev(csd)[1]
    n = length(incr)
    first = 1
    while((sum(incr) > sdf) && (first <= n)) {
        csi = cumsum(incr)[act]
        j = c(which(csi > csd), n)
        j = min(j[j >= first])
        aj = act[j]
        ajj = act[max(1, j - 1)]
        del = max(0, aj - ajj - 1) * min_incr / 2
        incr[ajj] = max(min_incr, incr[ajj] - csi[j] + csd[j] + del)
        first = j + 1
    }
    while ((si <- sum(incr)) > sdf) {
        j = max(which(incr > min_incr))
        jdx = seq_len(j)
        incr[j] = max(min_incr, incr[j] - (si - sdf))
    }
    new = cumsum(c(x[1], incr))
    x[ord] = new
    x
}

