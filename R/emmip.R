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
#' @param object An object of class \code{emm}, or a fitted model of a class
#'   supported by the \pkg{emmeans} package
#' @param formula Formula of the form 
#'   \code{trace.factors ~ x.factors | by.factors}. The EMMs are
#'   plotted against \code{x.factor} for each level of \code{trace.factors}.
#'   \code{by.factors} is optional, but if present, it determines separate
#'   panels. Each element of this formula may be a single factor in the model,
#'   or a combination of factors using the \code{*} operator.
#' @param type As in \code{\link{predict.emm}}, this determines
#'   whether we want to inverse-transform the predictions
#'   (\code{type = "response"}) or not (any other choice). The default is
#'   \code{"link"}, unless the \code{"predict.type"} option is in force; see
#'   \code{\link{emm_options}}.
#' @param pch The plotting characters to use for each group (i.e., levels of
#'   \code{trace.factors}). They are recycled as needed.
#' @param lty The line types to use for each group. Recycled as needed.
#' @param col The colors to use for each group, recycled as needed. If not
#'   specified, the default trellis colors are used.
#' @param plotit Logical value. If \code{TRUE}, the plot is displayed.
#'   Otherwise, one may use the \code{"lattice"} attribute of the returned
#'   object and print it, perhaps after additional manipulation.
#' @param ... Additional arguments passed to \code{\link{emmeans}} (when
#'   \code{object} is not already an \code{emm} object), or to
#'   \code{\link[lattice]{xyplot}}.
#'   
#' @section Details:
#' If \code{object} is a fitted model, \code{\link{emmeans}} is called with an
#' appropriate specification to obtain estimated marginal means for each
#' combination of the factors present in \code{formula} (in addition, any 
#' arguments in \code{\dots} that match \code{at}, \code{trend}, 
#' \code{cov.reduce}, or \code{fac.reduce} are passed to \code{emmeans}). 
#' Otherwise, if \code{object} is an \code{emm} object, its first element is 
#' used, and it must contain one estimate for each combination of the factors
#' present in \code{formula}.
#'
#' @return (Invisibly), a \code{\link{data.frame}} with the table of
#'   least-squares means that were plotted, with an additional \code{"lattice"}
#'   attribute containing the \code{trellis} object for the plot.
#'   
#' @note This function uses the \code{\link[lattice]{xyplot}} function in the
#'   \code{lattice} package (an error is returned if \code{lattice} is not
#'   installed). Conceptually, it is equivalent to
#'   \code{\link{interaction.plot}} where the summarization function is thought
#'   to return the
#'   EMMs.
#' 
#' @seealso \code{\link{emmeans}}, \code{\link{interaction.plot}}
#' @export
#' @method emmip default
#'
#' @examples
#' #--- Two-factor example
#' warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
#'
#' # Following plot is the same as the usual interaction plot of the data
#' emmip(warp.lm, wool ~ tension)
#'
#' #--- Three-factor example
#' noise.lm = lm(noise ~ size * type * side, data = auto.noise)
#'
#' # Separate interaction plots of size by type, for each side
#' emmip(noise.lm, type ~ size | side)
#'
#' # One interaction plot, using combinations of size and side as the x factor
#' emmip(noise.lm, type ~ side * size)
#'
#' # One interaction plot using combinations of type and side as the trace factor
#' # customize the colors, line types, and symbols to suggest these combinations
#' emmip(noise.lm, type * side ~ size, lty=1:2, col=1:2, pch=c(1,1,2,2))
#'
#' # 3-way interaction is significant, but doesn't much visual difference...
#' noise.lm2 = update(noise.lm, . ~ . - size:type:side)
#' emmip(noise.lm2, type * side ~ size, lty=1:2, col=1:2, pch=c(1,1,2,2))
#'
emmip.default = function(object, formula, type,  
        pch=c(1,2,6,7,9,10,15:20), lty = 1, col = NULL, plotit = TRUE, ...) {
    if (!requireNamespace("lattice"))
        stop("This function requires the 'lattice' package be installed.")
    if (length(formula) < 3)
        formula = .reformulate(as.character(formula)[[2]], response = ".single.")

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
    
    allvars = setdiff(.all.vars(formula), ".single.")
    emmopts$object = object
    emmopts$specs = .reformulate(allvars)
    emmo = do.call("emmeans", emmopts)
    if(missing(type)) {
        type = get_emm_option("summary")$predict.type
        if (is.null(type))
            type = .get.predict.type(emmo@misc)
    }
    type = .validate.type(type)

    emm = predict(emmo, type = type)
    emms = cbind(emmo@grid, emmean = emm)

    # Set up trace vars and key
    tvars = .all.vars(update(formula, . ~ 1))
    if (all(tvars == ".single.")) {
        emms$.single. = 1
        my.key = function(tvars) list()
    }
    else {
        my.key = function(tvars) 
            list(space="right", 
                 title = paste(tvars, collapse=" * "), 
                 points = TRUE, 
                 lines=length(lty) > 1,
                 cex.title=1)
    }
    tv = do.call(paste, emms[tvars])
    emms$tvar = factor(tv, levels=unique(tv))
    
    # figure out 'x' and 'by' vars
    rhs = strsplit(as.character(formula[3]), "\\|")[[1]]
    xvars = .all.vars(stats::reformulate(rhs[[1]]))
    xv = do.call(paste, emms[xvars])
    emms$xvar = factor(xv, levels = unique(xv))
    emms = emms[order(emms$xvar), ]
    plotform = emmean ~ xvar
    
    # see if we have any 'by' vars
    if (length(rhs) > 1) {
        byvars = .all.vars(stats::reformulate(rhs[[2]]))
        plotform = as.formula(paste("emmean ~ xvar |", paste(byvars, collapse="*")))
    }

    # The strips the way I want them
    my.strip = lattice::strip.custom(strip.names = c(TRUE,TRUE), strip.levels = c(TRUE,TRUE), sep = " = ")
    
    TP = TP.orig = lattice::trellis.par.get()
    TP$superpose.symbol$pch = pch
    TP$superpose.line$lty = lty
    if (!is.null(col)) TP$superpose.symbol$col = TP$superpose.line$col = col
    lattice::trellis.par.set(TP)
    
    xlab = ifelse(is.null(xargs$xlab),
        paste("Levels of", paste(xvars, collapse=" * ")), xargs$xlab)
    rspLbl = paste("Predicted", 
        ifelse(is.null(emmo@misc$inv.lbl), "response", emmo@misc$inv.lbl))
    ylab = ifelse(is.null(xargs$ylab),
        ifelse(type == "response", rspLbl, "Linear prediction"),
        xargs$ylab)
    
    # remove the unneeded stuff from xlabs
    xargs = xargs[setdiff(names(xargs), c("xlab","ylab"))]
    plotspecs = list(x = plotform, data = emms, groups = ~ tvar, 
        xlab = xlab, ylab = ylab,
        strip = my.strip, auto.key = my.key(tvars), type=c("p","l"))
    grobj = do.call(lattice::xyplot, c(plotspecs, xargs))
    if (plotit)
        print(grobj)
    attr(emms, "lattice") = grobj
    lattice::trellis.par.set(TP.orig)
    invisible(emms)
}
