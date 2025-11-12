##############################################################################
#    Copyright (c) 2012-2025 Russell V. Lenth                                #
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
#'   In addition, the user may specify \code{type = "scale"} to create a
#'   transformed scale for the vertical axis based on \code{object}'s response 
#'   transformation or link function. 
#' @param CIs Logical value. If \code{TRUE}, confidence intervals (or HPD intervals
#'   for Bayesian models) are added to the plot 
#'   (works only with \code{engine = "ggplot"}).
#' @param PIs Logical value. If \code{TRUE}, prediction intervals are added to the plot 
#'   (works only with \code{engine = "ggplot"}). This is allowed only if the
#'   underlying model family is \code{"gaussian"}.
#'   If both \code{CIs} and
#'   \code{PIs} are \code{TRUE}, the prediction intervals will be somewhat
#'   longer, lighter, and thinner than the confidence intervals. Additional
#'   parameters to \code{\link{predict.emmGrid}} (e.g., \code{sigma}) may be passed via
#'   \code{...}. For Bayesian models, PIs require \code{frequentist = TRUE} and 
#'   a value for \code{sigma}.
#' @param style Optional character value. This has an effect only when the
#'   horizontal variable is a single numeric variable. If \code{style} is
#'   unspecified or \code{"numeric"}, the horizontal scale will be numeric and
#'   curves are plotted using lines (and no symbols). With \code{style =
#'   "factor"}, the horizontal variable is treated as the levels of a factor
#'   (equally spaced along the horizontal scale), and curves are plotted using
#'   lines and symbols. When the horizontal variable is character or factor, or
#'   a combination of more than one predictor, \code{"factor"} style is always used.
#' @param engine Character value matching \code{"ggplot"} (default), 
#'   \code{"lattice"}, or \code{"none"}. The graphics engine to be used to produce the plot.
#'   These require, respectively, the \pkg{ggplot2} or \pkg{lattice} package to
#'   be installed. Specifying \code{"none"} is equivalent to setting \code{plotit = FALSE}.
#' @param plotit Logical value. If \code{TRUE}, a graphical object is returned;
#'   if \code{FALSE}, a data.frame is returned containing all the values
#'   used to construct the plot.
#' @param nesting.order Logical value. If \code{TRUE}, factors that are nested
#'   are presented in order according to their nesting factors, even if those nesting
#'   factors are not present in \code{formula}. If \code{FALSE}, only the
#'   variables in \code{formula} are used to order the variables.
#' @param ... Additional arguments passed to \code{\link{emmeans}} (when
#'   \code{object} is not already an \code{emmGrid} object),
#'   \code{predict.emmGrid}, 
#'   \code{emmip_ggplot}, or \code{emmip_lattice}.
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
#' @seealso \code{\link{emmeans}}, \code{\link{interaction.plot}}.
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
#' # ... with added confidence intervals and some formatting changes
#' emmip(noise.lm, type ~ side * size, CIs = TRUE,
#'     CIarg = list(linewidth = 1.5, alpha = 1, color = "orange"),
#'     dotarg = list(size = 2, shape = "square", color = "black"))
#'     
#' # Same using legacy theme
#' with_emm_options(gg.theme = 1,
#'     emmip(noise.lm, type ~ side * size, CIs = TRUE,
#'         CIarg = list(linewidth = 1.5, alpha = 1, color = "orange"),
#'         dotarg = list(size = 2, shape = "square", color = "black")))
#' 
#' # Create a black-and-white version of above with different linetypes
#' # (Let the linetypes and symbols default to the palette)
#' emmip(noise.lm, type ~ side * size, CIs = TRUE, col = "black",
#'       linearg = list(), dotarg = list(size = 4), CIarg = list(alpha = 1)) +
#'     ggplot2::theme_bw()
#'
#' # One interaction plot using combinations of type and side as the trace factor
#' emmip(noise.lm, type * side ~ size)
#'
#' # Individual traces in panels
#' emmip(noise.lm, ~ size | type * side)
#' 
#' # Example for the 'style' argument
#' fib.lm = lm(strength ~ machine * sqrt(diameter), data = fiber)
#' fib.rg = ref_grid(fib.lm, at = list(diameter = c(12, 14, 15, 25, 36)))
#' emmip(fib.rg, machine ~ diameter)   # curves (because diameter is numeric)
#' emmip(fib.rg, machine ~ diameter, style = "factor")  # points and lines
#' 
#' # For an example using extra ggplot2 code, see 'vignette("messy-data")',
#' # in the section on nested models.
emmip.default = function(object, formula, type, CIs = FALSE, PIs = FALSE,
                         style,
                         engine = get_emm_option("graphics.engine"),
                         # pch = c(1,2,6,7,9,10,15:20), 
                         # lty = 1, col = NULL, 
                         plotit = TRUE, 
                         nesting.order = FALSE, ...) {
    object = .chk.list(object, ...)
    engine = match.arg(engine, c("ggplot", "lattice", "none"))
    if (engine == "ggplot")
        .requireNS("ggplot2",
                   "The 'ggplot' engine requires the 'ggplot2' package be installed.")
    else if (engine == "lattice")
        .requireNS("lattice", 
                "The 'lattice' engine requires the 'lattice' package be installed.")
    else
        plotit = FALSE
    
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
    
    # add possibility of type = "scale". If so, we use "response" and set a flag
    if(missing(type)) {
        type = get_emm_option("summary")$predict.type
        if (is.null(type))
            type = .get.predict.type(emmo@misc)
    }
    # If we say type = "scale", set it to "response" and set a flag
    if (nonlin.scale <- (type %.pin% "scale")) 
        type = "response"
            
    type = .validate.type(type)
    
    # get point.est & frequentist if specified (affects only Bayesian models)
    point.est = (\(point.est = "median", ...) point.est)(...)
    frequentist = (\(frequentist = FALSE, ...) frequentist)(...)
    emms = summary(emmo, type = type, infer = c(CIs, FALSE), 
                   point.est = point.est, frequentist = frequentist)
    if(PIs) {
        prd = predict(emmo, interval = "pred", type = type, ...)
        if(!is.null(prd)) {
            emms$LPL = prd$lower.PL
            emms$UPL = prd$upper.PL
        }
    }
    # Ensure the estimate is named "yvar" and the conf limits are "LCL" and "UCL"
    nm = names(emms)
    tgts = c(attr(emms, "estName"), attr(emms, "clNames"))
    subs = c("yvar", "LCL", "UCL")
    for (i in 1:3)
        names(emms)[nm == tgts[i]] = subs[i] 
    attr(emms, "estName") = "yvar"
    if(!CIs)
        emms$LCL = emms$UCL = NULL
    
    if(!nesting.order) { # re-order by factor levels actually in plot
        snm = intersect(nm, unlist(specs))
        ord = do.call(order, unname(emms[rev(snm)]))
        emms = emms[ord, ]
    }
    
    sep = get_emm_option("sep")
    
    # Set up trace vars and key
    tvars = specs$lhs
    if (one.trace <- (length(tvars) == 0)) {
        tlab = ""
        tvars = ".single."
        emms$.single. = 1
    }
    else 
        tlab = paste(tvars, collapse = sep)
    tv = do.call(paste, c(unname(emms[tvars]), sep = sep))
    emms$tvar = factor(tv, levels=unique(tv))
    
    xvars = specs$rhs
    xv = do.call(paste, c(unname(emms[xvars]), sep = sep))
    ltest = max(apply(table(xv,tv), 2, function(x) sum(x > 0))) # length of longest trace
    if (!missing(style))
        styl = match.arg(style, c("factor", "numeric"))
    if (missing(style) || styl == "numeric")
        styl = ifelse(length(xvars) == 1 && 
                           is.numeric(emms[[xvars]]) &&
                           ltest > 1,
                   "numeric", "factor")
    if (styl == "factor") {
        emms$xvar = factor(xv, levels = unique(xv))
        predicate = "Levels of "
        if (ltest <= 1)
            message("Suggestion: Add 'at = list(", xvars, " = ...)' ",
                    "to call to see > 1 value per group.")
    }
    else {
        emms$xvar = as.numeric(xv)
        predicate = ""
    }
    emms = emms[order(emms$xvar), ]

    byvars = specs$by
    xlab = ifelse(is.null(xargs$xlab),
                  paste0(predicate, paste(xvars, collapse = sep)), xargs$xlab)
    rspLbl = paste("Predicted", 
                   ifelse(is.null(emmo@misc$inv.lbl), "response", emmo@misc$inv.lbl))
    ylab = ifelse(is.null(xargs$ylab),
                  ifelse(type == "response", rspLbl, "Linear prediction"),
                  xargs$ylab)
    
    # remove the unneeded stuff from xlabs
    xargs = xargs[setdiff(names(xargs), c("xlab","ylab"))]
    
    emms$.single. = NULL   # in case we have that trick column
    attr(emms, "labs") = list(xlab = xlab, ylab = ylab, tlab = tlab)
    attr(emms, "vars") = list(byvars = byvars, tvars = setdiff(tvars, ".single."))
    if (!plotit && nonlin.scale)
        attr(emms, "scale") = .make.scale(emmo@misc)
    
    if (!plotit || engine == "none")
        return (emms)
    
    fcn = paste("emmip", engine, sep = "_")
    args = c(list(emms = emms, style = styl), xargs)
    if (nonlin.scale)
        args = c(args, list(scale = .make.scale(emmo@misc)))
    do.call(fcn, args)
}

### render emmip using ggplot
#' @rdname emmip
#' @param dodge Numerical amount passed to \code{ggplot2::position_dodge} 
#'   by which points and intervals are offset so they do not collide.
#' @param xlab,ylab,tlab Character labels for the horizontal axis, vertical
#'   axis, and traces (the different curves), respectively. The \code{emmip}
#'   function generates these automatically and provides therm via the \code{labs} 
#'   attribute, but the user may override these if desired.
#' @param facetlab Labeller for facets (when by variables are in play). Use
#'   \code{"label_value"} to show just the factor levels, or \code{"label_both"}
#'   to show both the factor names and factor levels; \code{"label_context"}
#'   decides which based on how many \code{by} factors there are. See the
#'   documentation for \code{ggplot2::labellers}.
#' @param scale If not missing, an object of class \code{scales::trans} specifying
#'   a (usually) nonlinear scaling for the vertical axis. For example, 
#'   \code{scales = scales::log_trans()} specifies a logarithmic scale. For
#'   fine-tuning purposes, additional
#'   arguments to \code{ggplot2::scale_y_continuous} may be included in \code{...} .
#' @param dotarg \code{list}
#'   of arguments passed to \code{geom_point} to customize appearance of points.
#'   If this includes \code{shape} and it is of length greater than 1, it is interpreted
#'   as a shape palette instead.
#' @param linearg \code{list}
#'   of arguments passed to \code{geom_line} to customize appearance of lines.
#'   If this includes \code{linetype} and it is of length greater than 1,
#'   it is interpreted as a linetype.
#' @param CIarg,PIarg \code{list}s
#'   of arguments passed to \code{geom_linerange} to customize appearance of intervals.
#'   (Note: the \code{linetype} aesthetic defaults to \code{"solid"} under the hood)
#'   
#' @section Rendering functions:
#' The functions \code{emmip_ggplot} and \code{emmip_lattice}
#' are called when \code{plotit == TRUE} to render the plots; 
#' but they may also be called later on an object saved via \code{plotit = FALSE}
#' (or \code{engine = "none"}). The functions require that \code{emms} contains variables
#' \code{xvar}, \code{yvar}, and \code{tvar}, and attributes \code{"labs"} and \code{"vars"}.
#' Confidence intervals are plotted if variables \code{LCL} and \code{UCL} exist;
#' and prediction intervals are plotted if \code{LPL} and \code{UPL} exist.
#' Finally, it must contain the variables named in \code{attr(emms, "vars")}.
#' 
#' In \code{emmip_ggplot}, colors, linetypes, and shapes are all assigned to
#' groups (according to \code{tvar}) unless overridden. So, for example, one may 
#' have different symbols for each group by simply specifying \code{dotarg = list()}.
#' 
#' 
#' @examples
#' 
#'### Options with transformations or link functions
#' neuralgia.glm <- glm(Pain ~ Treatment * Sex + Age, family = binomial(), 
#'                      data = neuralgia) 
#' 
#' # On link scale:
#' emmip(neuralgia.glm, Treatment ~ Sex)
#' 
#' # On response scale:
#' emmip(neuralgia.glm, Treatment ~ Sex, type = "response")
#' 
#' # With transformed axis scale and custom scale divisions
#' emmip(neuralgia.glm, Treatment ~ Sex, type = "scale",
#'     breaks = seq(0.10, 0.90, by = 0.10))
#' @export
#' @importFrom rlang .data
emmip_ggplot = function(emms, style = "factor", dodge = .2,
                        xlab = labs$xlab, ylab = labs$ylab, tlab = labs$tlab,
                        facetlab = "label_both",
                        scale, 
                        dotarg = list(shape = "circle"), 
                        linearg = list(linetype = "solid"),
                        CIarg = list(alpha = .40),
                        PIarg = list(alpha = .25),
                        col = NULL,
                        ...) {
    
    thm = theme_emm()    # Depends on gg.theme option
    
    labs = attr(emms, "labs")
    vars = attr(emms, "vars")
    ngrps = ifelse(is.null(emms$tvar), 1, length(unique(emms$tvar)))
    
    shape.pal = linetype.pal = NULL   # we use these to store shape and linetype specs in dotarg and linearg
    
    CIs = !is.null(emms$LCL)
    PIs = !is.null(emms$LPL)
    pos = ggplot2::position_dodge(width = ifelse(CIs|PIs, dodge, 0)) # use dodging if CIs
    
    if(missing(scale) && !is.null(attr(emms, "scale")))
        scale = attr(emms, "scale")

    dotarg$position = pos
    if(is.null(dotarg$size)) 
        dotarg$size = 3
    
    if (ngrps > 1) {
        if(!is.null(dotarg$shape) && length(dotarg$shape) > 1) {  # treat as a shape palette
            shape.pal = dotarg$shape
            dotarg$shape = NULL
        }
        
        linearg$mapping = ggplot2::aes(x = as.numeric(.data$xvar), linetype = .data[[tlab]])
        dotarg$mapping = ggplot2::aes(shape = .data[[tlab]])
        
        linearg$position = pos
        if(is.null(linearg$linewidth))
            linearg$linewidth = 0.8
        if(!is.null(linearg$linetype) && length(linearg$linetype) > 1) {  # treat as a linetype palette
            linetype.pal = linearg$linetype
            linearg$linetype = NULL
        }
        
        emms[[tlab]] = emms$tvar   # make the trace column have same name as its label
        
        grobj = ggplot2::ggplot(emms, ggplot2::aes(x = .data$xvar, y = .data$yvar, group = .data[[tlab]], 
                    color = .data[[tlab]]))
        
        if (style == "factor") {
            grobj = grobj + do.call(ggplot2::geom_point, dotarg) +
                ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = 0.1))
        }
        grobj = grobj +
                do.call(ggplot2::geom_line, linearg)
        
        # handle any custom shapes or linetypes
        if(!is.null(shape.pal))
            grobj = grobj + ggplot2::scale_shape_manual(values = shape.pal)
        if(!is.null(linetype.pal))
            grobj = grobj + ggplot2::scale_linetype_manual(values = linetype.pal)
    }
    else { # just one trace per plot
        if(is.null(col))
            col = .emm_palette[1]
        linearg$color = dotarg$color = CIarg$color = PIarg$color = col
        col = NULL
        grobj = ggplot2::ggplot(emms, ggplot2::aes(x = .data$xvar, y = .data$yvar))
        if (style == "factor")
            grobj = grobj + do.call(ggplot2::geom_point, dotarg) + 
                ggplot2::scale_x_discrete(expand = ggplot2::expansion(mult = 0.1))
        linearg$mapping = ggplot2::aes(x = as.numeric(.data$xvar))
        grobj = grobj +
            do.call(ggplot2::geom_line, linearg)
        
    }
    grobj = grobj + ggplot2::labs(x = xlab, y = ylab)
    if (PIs) {
        PIarg$mapping = ggplot2::aes(ymin = .data$LPL, ymax = .data$UPL)
        PIarg$position = pos
        if(is.null(PIarg$linetype))  PIarg$linetype = "solid"
        if(is.null(PIarg$linewidth))  PIarg$linewidth = 1.5
        grobj = grobj + do.call(ggplot2::geom_linerange, PIarg)
    }
    if (CIs) {
        CIarg$mapping = ggplot2::aes(ymin = .data$LCL, ymax = .data$UCL)
        CIarg$position = pos
        if(is.null(CIarg$linetype))  CIarg$linetype = "solid"
        if(is.null(CIarg$linewidth))  CIarg$linewidth = 2.5
        grobj = grobj + do.call(ggplot2::geom_linerange, CIarg)
    }
    if (length(byvars <- vars$byvars) > 0) {  # we have by variables 
        if (length(byvars) > 1) {
            byform = as.formula(paste(byvars[1], " ~ ", paste(byvars[-1], collapse="*")))
            grobj = grobj + ggplot2::facet_grid(byform, labeller = facetlab)
        }
        else
            grobj = grobj + ggplot2::facet_wrap(byvars, labeller = facetlab)
    }
    if (!missing(scale)) {
        args = list(...)
        pass = names(args) %.pin% names(as.list(args(ggplot2::scale_y_continuous)))
        args = c(list(trans = scale), args[pass])
        grobj = grobj + do.call(ggplot2::scale_y_continuous, args)
    }
    if (style == "factor")
        grobj = grobj + do.call(ggplot2::geom_point, dotarg)
    
    if(!is.null(col))
        grobj = grobj + ggplot2::scale_color_manual(values = rep(col, ngrps))
    
    
    grobj + thm
}

#' @importFrom grDevices palette.colors
##### Here is a palette of 45 colors. More than that, SOL
.emm_palette = c(
    grDevices::palette.colors(palette = "Okabe-Ito)")[c(6:9,1:5)], # Start w/ steel blue, yellow is last
    grDevices::palette.colors(palette = "Polychrome 36")[sample(1:36)]
)

theme_emm = function (base_size = 13, base_family = "sans", header_family = "sans", 
                      base_line_size = base_size/22, base_rect_size = base_size/22, 
                      ink ="#0e0033ff", paper = "white", accent = "#FF4000") 
{
    if(inherits(thm <- get_emm_option("gg.theme"), "theme"))
        return(thm)
    
    if(thm == 1) # legacy theme from version 1.x.x
        return(ggplot2::theme_grey())
    
    dark_color = "#4E4369"
    mid_color =  "#79718D"   ###   "#002789"
    lit_color =  "#d2d2d8ff"
  
    ggplot2::theme_light(base_size = base_size, base_family = base_family,
            header_family = header_family, base_line_size = base_line_size,
            base_rect_size = base_rect_size, ink = ink, paper = paper, accent = accent) +
      ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                    panel.border = ggplot2::element_rect(color = mid_color, linewidth = 0.85),
                    panel.grid.major = ggplot2::element_line(colour = lit_color),
                    strip.background = ggplot2::element_rect(fill = mid_color, color = mid_color),
                    strip.text = ggplot2::element_text(size = ggplot2::rel(0.8), family = "sans", face = "bold", lineheight = 1.1),
                    plot.title = ggplot2::element_text(colour = ink, family = "serif", size = ggplot2::rel(1.8), margin = ggplot2::margin(12, 0, 8, 0)),
                    plot.subtitle = ggplot2::element_text(size = ggplot2::rel(1.1), family = "serif", margin = ggplot2::margin(4, 0, 0, 0), color = dark_color),
                    legend.justification = "top",
                    palette.colour.discrete = .emm_palette # color-blind friendly, 45 total colors
    )
}


#' @rdname emmip
#' @param emms A \code{data.frame} created by calling \code{emmip} with
#'   \code{plotit = FALSE}. Certain variables and attributes are expected
#'   to exist in this data frame; see the section detailing the rendering functions.
#' @param pch (Lattice only) The plotting characters to use for each group (i.e., levels of
#'   \code{trace.factors}). They are recycled as needed.
#' @param lty (Lattice only) The line types to use for each group. Recycled as needed.
#' @param col With \code{emmip_ggplot}, this adds \code{color = col} (not
#'   \code{colour}) to all of the \code{*arg} lists. This is intended for setting a
#'   common color for everything, such as a black-and-white plot. 
#'   With \code{emmip_lattice}, \code{col} specifies the colors to use
#'   for each group, recycled as needed. If not specified, the default trellis
#'   colors are used.
#' @export
emmip_lattice = function(emms, style = "factor", 
                         xlab = labs$xlab, ylab = labs$ylab, tlab = labs$tlab, 
                         pch = c(1,2,6,7,9,10,15:20), 
                         lty = 1, col = NULL, ...) {
    labs = attr(emm, "labs")
    vars = attr(emms, "vars")
    
    # The strips the way I want them
    my.strip = lattice::strip.custom(strip.names = c(TRUE,TRUE), strip.levels = c(TRUE,TRUE), sep = " = ")
    
    if (length(vars$byvars) == 0)
        plotform = yvar ~ xvar
    else
        plotform = as.formula(paste("yvar ~ xvar |", paste(vars$byvars, collapse="*")))
    
    sep = get_emm_option("sep")
    my.key = function(tvars) 
        list(space="right", 
             title = paste(tvars, collapse = sep), 
             points = TRUE, 
             lines=length(lty) > 1,
             cex.title=1)
    
    TP = TP.orig = lattice::trellis.par.get()
    TP$superpose.symbol$pch = pch
    TP$superpose.line$lty = lty
    if (!is.null(col)) TP$superpose.symbol$col = TP$superpose.line$col = col
    lattice::trellis.par.set(TP)
    
    plty = if(style=="factor") c("p","l")   else "l"
    plotspecs = list(x = plotform, data = emms, groups = ~ tvar, 
                     xlab = xlab, ylab = ylab,
                     strip = my.strip, auto.key = my.key(vars$tvars), 
                     type = plty)
    if(length(vars$tvars) == 0)
        plotspecs$auto.key = NULL # no key when single trace
    grobj = do.call(lattice::xyplot, c(plotspecs, list(...)))
    lattice::trellis.par.set(TP.orig)
    
    grobj
}


