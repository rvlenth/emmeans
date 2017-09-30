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

# contrast() and related functions (previously in with emmeans code)

### 'contrast' S3 generic and method
#' Contrasts and linear functions of EMMs
#' 
#' These methods provide for follow-up analyses of \code{emm} objects:
#' Contrasts, pairwise comparisons, tests, and confidence intervals. They may
#' also be used to compute arbitrary linear functions of predictions or EMMs.
#'
#' @export
contrast = function(object, ...)
    UseMethod("contrast")

#' @rdname contrast 
#' @param object An object of class \code{emm}
#' @param method Character value giving the root name of a contrast method (e.g.
#'   \code{"pairwise"} -- see \link{emmc-functions}). Alternatively, a named 
#'   \code{list} of coefficients (for a contrast or linear function) that must
#'   each conform to the number of results in each \code{by} group.
#'   In a multi-factor situation, the factor levels are combined and treated
#'   like a single factor.
#' @param interaction Character vector or logical value. If this is specified,
#'   \code{method} is ignored. See the \dQuote{Interaction contrasts} section
#'   below for details.
#' @param by Character names of variable(s) to be used for ``by'' groups. The
#'   contrasts or joint tests will be evaluated separately for each combination
#'   of these variables. If \code{object} was created with by groups, those are
#'   used unless overridden. Use \code{by = NULL} to use no by groups at all.
#' @param offset Numeric vector of the same length as each \code{by} group.
#'   These values are added to their respective linear estimates. (It is ignored
#'   when \code{interaction} is specified.)
#' @param name Character name to use to override the default label for contrasts
#'   used in table headings or subsequent contrasts of the returned object.
#' @param options If non-\code{NULL}, a named \code{list} of arguments to pass
#'   to \code{\link{update.emm}}, just after the object is constructed.
#' @param type Character: prediction type (e.g., \code{"response"}) -- added to \code{options}
#' @param adjust Character: adjustment method (e.g., \code{"bonferroni"}) -- added to \code{options}
#' @param ... Additional arguments passed to other methods
#'
#' @return \code{contrast} and \code{pairs} return an object of class \code{emm}. Its grid will correspond to the levels of the contrasts and any \code{by} variables.
#' 
#' @section Pairs method:
#' The call \code{pairs(object)} is equivalent to \code{contrast(object, method = "pairwise")}; and \code{pairs(object, reverse = TRUE)} is the same as \code{contrast(object, method = "revpairwise")}.
#' 
#' @section Interaction contrasts:
#' When \code{interaction} is specified, interaction
#' contrasts are computed: Contrasts are generated for each factor separately,
#' one at a time; and these contrasts are applied to the object (the first time
#' around) or to the previous result (subsequently). (Any factors specified in
#' \code{by} are skipped.) The final result comprises contrasts of contrasts,
#' or, equivalently, products of contrasts for the factors involved. Processing
#' is done in the order of appearance in \code{object@levels}. With
#' \code{interaction = TRUE}, \code{method} (if specified as character) is used
#' for each contrast. If \code{interaction} is a character vector, the elements
#' specify the respective contrast method(s); they are recycled as needed.
#' 
#' @note When \code{object} has a nesting structure (this can be seen via
#'   \code{str(object)}), then any grouping factors involved are forced into
#'   service as \code{by} variables, and the contrasts are thus computed
#'   separately in each nest. This in turn may lead to an irregular grid in the
#'   returned \code{emm} object, which may not be valid for subsequent
#'   \code{emmeans} calls.
#' 
#' @method contrast emm
#' @export
#'
#' @examples
#' warp.lm <- lm(breaks ~ wool*tension, data = warpbreaks)
#' warp.emm <- emmeans(warp.lm, ~ tension | wool)
#' contrast(warp.emm, "poly")    # inherits 'by = "wool"' from warp.emm
#' pairs(warp.emm)               # ditto
#' contrast(warp.emm, "eff", by = NULL)  # contrasts of the 6 factor combs
#' 
#' # An interaction contrast for tension:wool
#' tw.emm <- contrast(warp.emm, interaction = c("poly", "consec"), by = NULL)
#' tw.emm          # see the estimates
#' coef(tw.emm)    # see the contrast coefficients
contrast.emm = function(object, method = "eff", interaction = FALSE, 
                        by, offset = NULL, name = "contrast", 
                        options = get_emm_option("contrast"), 
                        type, adjust, ...) 
{
    if(missing(by)) 
        by = object@misc$by.vars
    if(length(by) == 0) # character(0) --> NULL
        by = NULL
    
    nesting = object@model.info$nesting
    if (!is.null(nesting) || !is.null(object@misc$display))
        return (.nested_contrast(rgobj = object, method = method, by = by, adjust = adjust, ...))
    
    orig.grid = object@grid[, , drop = FALSE]
    orig.grid[[".wgt."]] = orig.grid[[".offset."]] = NULL
    
    if (is.logical(interaction) && interaction)
        interaction = method
    if (!is.logical(interaction)) { # i.e., interaction is not FALSE
        if (!is.character(interaction))
            stop("interaction requires named contrast function(s)")
        if(missing(adjust))
            adjust = "none"
        by = NULL
        vars = names(object@levels)
        k = length(vars)
        if(!is.null(by)) {
            vars = c(setdiff(vars, by), by)
            k = k - length(by)
        }
        interaction = rep(interaction, k)[1:k]
        tcm = NULL
        for (i in k:1) {
            nm = paste(vars[i], interaction[i], sep = "_")
            object = contrast.emm(object, interaction[i], by = vars[-i], name = nm)
            if(is.null(tcm))
                tcm = object@misc$con.coef
            else
                tcm = object@misc$con.coef %*% tcm
            vars[i] = nm
        }
        object = update(object, by = by, adjust = adjust, ...)
        object@misc$is.new.rg = NULL
        object@misc$orig.grid = orig.grid
        object@misc$con.coef = tcm
        if(!is.null(options)) {
            options$object = object
            object = do.call(update.emm, options)
        }
        return(object)
    }
    
    # else
    linfct = object@linfct[, , drop = FALSE]
    args = g = object@grid[, , drop = FALSE]
    args[[".offset."]] = NULL 
    args[[".wgt."]] = NULL # ignore auxiliary stuff in labels, etc.
    if (!is.null(by)) {
        by.rows = .find.by.rows(args, by)
        bylevs = args[, by, drop=FALSE]
        args = args[by.rows[[1]], , drop=FALSE]
        for (nm in by) args[[nm]] = NULL
    }
    args$sep = ","
    levs = do.call("paste", args)  # NOTE - these are levels for the first (or only) by-group
    
    
    if (is.list(method)) {
        cmat = as.data.frame(method, optional = TRUE)
        # I have no clue why they named that argument 'optional',
        # but setting it to TRUE keeps it from messing up the names
        method = function(levs) cmat
    }
    else if (is.character(method)) {
        fn = paste(method, "emmc", sep=".")
        method = if (exists(fn, mode="function")) 
            get(fn) 
        else 
            stop(paste("Contrast function '", fn, "' not found", sep=""))
    }
    # case like in old lsmeans, contr = list
    else if (!is.function(method))
        stop("'method' must be a function or the basename of an '.emmc' function")
    
    # Get the contrasts; this should be a data.frame
    cmat = method(levs, ...)
    if (!is.data.frame(cmat))
        stop("Contrast function must provide a data.frame")
    else if(ncol(cmat) == 0)
        cmat = data.frame(`(nothing)` = rep(NA, nrow(args)), check.names = FALSE)
    # warning("No contrasts were generated! Perhaps only one emmean is involved.\n",
    #      "  This can happen, for example, when your predictors are not factors.")
    else if (nrow(cmat) != nrow(args))
        stop("Nonconforming number of contrast coefficients")
    tcmat = t(cmat)
    
    if (is.null(by)) {
        linfct = tcmat %*% linfct
        grid = data.frame(.contrast.=names(cmat))
        if (hasName(object@grid, ".offset."))
            grid[[".offset."]] = t(cmat) %*% object@grid[[".offset."]]
        by.rows = list(seq_along(object@linfct[ , 1]))
    }
    
    # NOTE: The kronecker thing here depends on the grid being regular.
    # Irregular grids are handled by .neted_contrast
    else {
        tcmat = kronecker(.diag(rep(1,length(by.rows))), tcmat)
        linfct = tcmat %*% linfct[unlist(by.rows), , drop = FALSE]
        tmp = expand.grid(con = names(cmat), by = seq_len(length(by.rows)))###unique(by.id))
        grid = data.frame(.contrast. = tmp$con)
        n.each = ncol(cmat)
        row.1st = sapply(by.rows, function(x) x[1])
        xlevs = list()
        for (v in by)
            xlevs[[v]] = rep(bylevs[row.1st, v], each=n.each)
        grid = cbind(grid, as.data.frame(xlevs))
        if (hasName(object@grid, ".offset."))
            grid[[".offset."]] = tcmat %*% object@grid[unlist(by.rows), ".offset."]
    }
    
    # Rename the .contrast. column -- ordinarily to "contrast",
    # but otherwise a unique variation thereof
    con.pat = paste("^", name, "[0-p]?", sep = "")
    n.prev.con = length(grep(con.pat, names(grid)))
    con.col = grep("\\.contrast\\.", names(grid))
    con.name = paste(name, 
                     ifelse(n.prev.con == 0, "", n.prev.con), sep="")
    names(grid)[con.col] = con.name
    
    row.names(linfct) = NULL
    misc = object@misc
    misc$initMesg = NULL # initial annotation likely will no longer apply
    misc$estName = "estimate"
    if (!is.null(et <- attr(cmat, "type")))
        misc$estType = et
    else {
        is.con = all(abs(sapply(cmat, sum)) < .001)
        misc$estType = ifelse(is.con, "contrast", "prediction")
    }
    misc$methDesc = attr(cmat, "desc")
    misc$famSize = size = length(by.rows[[1]])
    misc$pri.vars = setdiff(names(grid), c(".offset.",".wgt."))
    if (missing(adjust)) adjust = attr(cmat, "adjust")
    if (is.null(adjust)) adjust = "none"
    if (!is.null(attr(cmat, "offset")))
        offset = attr(cmat, "offset")
    if (!is.null(offset)) {
        if(!hasName(grid, ".offset."))
            grid[[".offset."]] = 0
        grid[[".offset."]] = grid[[".offset."]] + rep(offset, length(by.rows))
    }
    misc$adjust = adjust
    misc$infer = c(FALSE, TRUE)
    misc$by.vars = by
    # save contrast coefs
    by.cols = seq_len(ncol(tcmat))
    if(!is.null(by.rows))
        by.cols[unlist(by.rows)] = by.cols # gives us inverse of by.rows order
    misc$orig.grid = orig.grid  # save original grid
    misc$con.coef = tcmat[ , by.cols, drop = FALSE] # save contrast coefs
    # zap the transformation info except in special cases
    if (!is.null(misc$tran)) {
        misc$orig.tran = misc$tran
        true.con = all(zapsmall(apply(cmat, 2, sum)) == 0) # each set of coefs sums to 0
        if (true.con && misc$tran %in% c("log", "genlog", "logit")) {
            misc$log.contrast = TRUE      # remember how we got here; used by summary
            misc$orig.inv.lbl = misc$inv.lbl
            if (misc$tran == "logit") {
                misc$inv.lbl = "odds.ratio"
                misc$tran = "log.o.r."
            }
            else {
                misc$inv.lbl = "ratio"
                misc$tran = "log"
            }
        }
        else
            misc$tran = misc$tran.mult = NULL
    }
    
    # ensure we don't inherit inappropriate settings
    misc$null = misc$delta = misc$side = NULL
    
    object@roles$predictors = "contrast"
    levels = list()
    for (nm in setdiff(names(grid), ".offset."))
        levels[[nm]] = unique(grid[[nm]])
    
    ### bypass new as we're not re-classing    result = new("emm", object, linfct = linfct, levels = levels, grid = grid, misc = misc)
    result = as(object, "emm")
    result@linfct = linfct
    result@levels = levels
    result@grid = grid
    result@misc = misc
    result@roles$predictors = setdiff(names(result@levels), result@roles$multresp)
    
    if (!missing(type))
        options = as.list(c(options, predict.type = type))
    if(!is.null(options)) {
        options$object = result
        result = do.call("update.emm", options)
    }
    result
}


# pairs method

#' @rdname contrast 
#' @param x An \code{emm} object
#' @param reverse Logical value - determines whether to use \code{"pairwise"} (if \code{TRUE}) or \code{"revpairwise"} (if \code{FALSE}).
#' @inheritParams contrast.emm 
#' @importFrom graphics pairs
#' @export
pairs.emm = function(x, reverse = FALSE, ...) {
    object = x # for my sanity
    if (reverse)
        contrast(object, method = "revpairwise", ...)
    else
        contrast(object, method = "pairwise", ...)
}


# coef method - returns contrast coefficients
#' @rdname contrast 
#' @return \code{coef} returns a \code{data.frame} containing the object's grid, along with columns named \code{c.1, c.2, ...} containing the contrast coefficients. If 
#' @export
#' @importFrom stats coef
#' @method coef emm
coef.emm = function(object, ...) {
    if (is.null(cc <- object@misc$con.coef)) {
        message("No contrast coefficients are available")
        return (NULL)
    }
    cc = as.data.frame(t(cc))
    names(cc) = paste("c", seq_len(ncol(cc)), sep = ".")
    cbind(object@misc$orig.grid, cc)
}

