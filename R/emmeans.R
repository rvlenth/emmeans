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

# emmeans and related functions


# emmeans utility for a list of specs
emmeans.list = function(object, specs, ...) {
    result = list()
    nms = names(specs)
    # Format a string describing the results
    .make.desc = function(meth, pri, by) {
        pri = paste(pri, collapse = ", ")
        desc = paste(meth, "of", pri)
        if (!is.null(by)) {
            by = paste(by, collapse = ", ")
            desc = paste(desc, "|", by)
        }
        desc
    }
    name.arg = match.call()$name
    
    for (i in seq_len(length(specs))) {
        nm = nms[i]
        # We'll rename the contrasts if spec is named
        if (is.null(name.arg)) {
            contr.name = ifelse(nm == "", "contrast", paste0(nm, ".contrast"))
            res = emmeans(object=object, specs = specs[[i]], name = contr.name, ...)
        }
        else
            res = emmeans(object=object, specs = specs[[i]], ...)
        if (is.data.frame(res)) { # happens e.g. when cld is used
            if (is.null(nm))
                nm = .make.desc("summary", attr(res, "pri.vars"), attr(res, "by.vars"))
            result[[nm]] = res
        }
        else if (is.list(res)) {
            for (j in seq_len(length(res))) {
                m = res[[j]]@misc
                if (is.null(nm))
                    names(res)[j] = .make.desc(m$methDesc, res[[1]]@misc$pri.vars, m$by.vars)
                else
                    names(res)[j] = paste(nm, m$methDesc)
            }
            result = c(result,res)
        }
        else{
            if (is.null(nm))
                nm = .make.desc(res@misc$methDesc, res@misc$pri.vars, res@misc$by.vars)
            result[[nm]] = res
        }
    }
    class(result) = c("emm_list", "list")
    result
}



# Here's our flagship function!
#' Estimated marginal means (Least-squares means)
#' 
#' Compute estimated marginal means (EMMs) for specified factors
#' or factor combinations in a linear model; and optionally, comparisons or
#' contrasts among them. EMMs are also known as least-squares means.
#' 
#' Users should also consult the documentation for \code{\link{ref_grid}}, 
#' because many important options for EMMs are implemented there, via the 
#' \code{...} argument.
#'
#' @param object An object of class \code{emmGrid}; or a fitted model object
#'   that is supported, such as the result of a call to \code{lm} or
#'   \code{lmer}. Many fitted-model objects are supported; see
#'   \href{../doc/models.html}{\code{vignette("models", "emmeans")}} for details.
#' @param specs A \code{character} vector specifying the names of the predictors
#'   over which EMMs are desired. \code{specs} may also be a \code{formula}
#'   or a \code{list} (optionally named) of valid \code{spec}s. Use of formulas
#'   is described in the Overview section below. 
#'   Specifying \code{.} as the only factor name
#'   creates a list of specifications for all model terms.
#'   
#'   \bold{Note:} We recommend \emph{against} using two-sided formulas; see the
#'   note below for \code{contr}.
#' @param by A character vector specifying the names of predictors to condition on.
#' @param fac.reduce A function that combines the rows of a matrix into a single
#'   vector. This implements the ``marginal averaging'' aspect of EMMs. 
#'   The default is the mean of the rows. Typically if it is overridden,
#'   it would be some kind of weighted mean of the rows. If \code{fac.reduce} is
#'   nonlinear, bizarre results are likely, and EMMs will not be
#'   interpretable. NOTE: If the \code{weights} argument is non-missing,
#'   \code{fac.reduce} is ignored.
#' @param contr A character value or \code{list} specifying contrasts to be
#'   added. See \code{\link{contrast}}. \bold{Note:} \code{contr} is ignored when
#'   \code{specs} is a formula. \bold{Note 2:}: We recommend \emph{against} using this
#'   argument; obtaining means and obtaining contrasts are two different things,
#'   and it is best to do them in separate steps, using the \code{\link{contrast}}
#'   function for the contrasts.
#' @param options If non-\code{NULL}, a named \code{list} of arguments to pass
#'   to \code{\link{update.emmGrid}}, just after the object is constructed.
#'   (Options may also be included in \code{...}; see the \sQuote{options}
#'   section below.)
#' @param weights Character value, numeric vector, or numeric matrix specifying
#'   weights to use in averaging predictions. See \dQuote{Weights} section below.
#'   Also, if \code{object} is not already a reference grid, \code{weights}
#'   (if it is character) is passed to \code{ref_grid} as \code{wt.nuis} in case 
#'   nuisance factors are specified. We can override this by specifying 
#'   \code{wt.nuis} explicitly.
#'   This more-or-less makes the weighting of nuisance factors consistent with
#'   that of primary factors.
#' @param offset Numeric vector or scalar. If specified, this adds an offset to
#'   the predictions, or overrides any offset in the model or its
#'   reference grid. If a vector of length differing from the number of rows in 
#'   the result, it is subsetted or cyclically recycled.
#' @param ... When \code{object} is not already a \code{"emmGrid"}
#'   object, these arguments are passed to \code{\link{ref_grid}}. Common
#'   examples are \code{at}, \code{cov.reduce}, \code{data}, \code{type}, 
#'   \code{regrid}, \code{df}, \code{nesting}, and \code{vcov.}.
#'   Model-type-specific options (see
#'   \href{../doc/models.html}{\code{vignette("models", "emmeans")}}), commonly
#'   \code{mode}, may be used here as well. In addition, if the model formula
#'   contains references to variables that are not predictors, you must provide
#'   a \code{params} argument with a list of their names.
#'   These arguments may also be used in lieu of \code{options}. See the 
#'   \sQuote{Options} section below.
#' @param tran Placeholder to prevent it from being included in \code{...}.
#'   If non-missing, it is added to `options`. See the \sQuote{Options}
#'   section.
#'   
#' @return   When \code{specs} is a \code{character} vector or one-sided formula,
#'   an object of class \code{"emmGrid"}. A number of methods
#'   are provided for further analysis, including
#'   \code{\link{summary.emmGrid}}, \code{\link{confint.emmGrid}},
#'   \code{\link{test.emmGrid}}, \code{\link{contrast.emmGrid}},
#'   and \code{\link{pairs.emmGrid}}.

#' When \code{specs} is a \code{list} or a \code{formula} having a left-hand
#' side, the return value is an \code{\link{emm_list}} object, which is simply a
#' \code{list} of \code{emmGrid} objects. 
#' 
#' @section Overview:
#' Estimated marginal means or EMMs (sometimes called least-squares means) are
#' predictions from a linear model over a \emph{reference grid}; or marginal
#' averages thereof. The \code{\link{ref_grid}} function identifies/creates the
#' reference grid upon which \code{emmeans} is based.
#' 
#' For those who prefer the terms \dQuote{least-squares means} or
#' \dQuote{predicted marginal means}, functions \code{lsmeans} and
#' \code{pmmeans} are provided as wrappers. See \code{\link{wrappers}}.
#' 
#' If \code{specs} is a \code{formula}, it should be of the form \code{~ specs},
#' \code{~ specs | by}, \code{contr ~ specs}, or \code{contr ~ specs | by}. The
#' formula is parsed and the variables therein are used as the arguments
#' \code{specs}, \code{by}, and \code{contr} as indicated. The left-hand side is
#' optional (and we don't recommend it), but if specified it should be the 
#' name of a contrast family (e.g., \code{pairwise}). Operators like
#' \code{*} or \code{:} are needed in the formula to delineate names, but
#' otherwise are ignored.
#' 
#' We now also allow using \code{.} in \code{specs}. If this is done, we run
#' \code{\link{joint_tests}} on the side to determine all relevant model terms,
#' then replace \code{specs} with a corresponding list of specifications. This is
#' a convenience, but it can create a sizeable \code{emm_list} object and it is
#' coded rather inefficiently. While it is permissible to include contrasts
#' via a \code{contr} argument or formula left-hand-side, we recommend instead
#' doing this in a follow-up test with \code{\link{contrast.emm_list}}. \emph{Caution}:
#' In models with nested fixed effects, using \code{.} creates results where
#' nested factors interact with nesting factors; in those cases, any contrasts 
#' you specify will go across nests, which is likely not what is desired.
#' 
#' In the special case where the mean (or weighted mean) of all the predictions
#' is desired, specify \code{specs} as \code{~ 1} or \code{"1"}.
#' 
#' A number of standard contrast families are provided. They can be identified 
#' as functions having names ending in \code{.emmc} -- see the documentation
#' for \code{\link{emmc-functions}} for details -- including how to write your
#' own \code{.emmc} function for custom contrasts.
#' 
#' @section Weights:
#' If \code{weights} is a vector, its length must equal
#'   the number of predictions to be averaged to obtain each EMM.
#'   If a matrix, each row of the matrix is used in turn, wrapping back to the
#'   first row as needed.  When in doubt about what is being averaged (or how
#'   many), first call \code{emmeans} with \code{weights = "show.levels"}.
#'   
#' If \code{weights} is a string, it should partially match one of the following:
#' \describe{
#' \item{\code{"equal"}}{Use an equally weighted average.}
#' \item{\code{"proportional"}}{Weight in proportion to the frequencies (in the
#'   original data) of the factor combinations that are averaged over.}
#' \item{\code{"outer"}}{Weight in proportion to each individual factor's
#'   marginal frequencies. Thus, the weights for a combination of factors are the
#'   outer product of the one-factor margins}
#' \item{\code{"cells"}}{Weight according to the frequencies of the cells being
#'   averaged.}
#' \item{\code{"flat"}}{Give equal weight to all cells with data, and ignore
#'   empty cells.}
#' \item{\code{"show.levels"}}{This is a convenience feature for understanding
#'   what is being averaged over. Instead of a table of EMMs, this causes the
#'   function to return a table showing the levels that are averaged over, in the
#'   order that they appear.}
#' }
#' Outer weights are like the 'expected' counts in a chi-square test of
#' independence, and will yield the same results as those obtained by
#' proportional averaging with one factor at a time. All except \code{"cells"}
#' uses the same set of weights for each mean. In a model where the predicted
#' values are the cell means, cell weights will yield the raw averages of the
#' data for the factors involved. Using \code{"flat"} is similar to
#' \code{"cells"}, except nonempty cells are weighted equally and empty cells
#' are ignored.
#' 
#' @section Counterfactuals:
#' Counterfactual reference grids (see the documentation for \code{\link{ref_grid}})
#' contain pairs of imputed and actual factor levels, and are handled in a special way.
#' For starters, the \code{weights} argument is ignored and we always use 
#' \code{"cells"} weights.
#' Our understanding is that if factors \code{A, B} are specified as counterfactuals, 
#' the marginal means for \code{A} should still be the same as if \code{A} were the only 
#' counterfactual. Accordingly, in computing these marginal means, we
#' exclude all cases where \code{B != actual_B}, because if \code{A} were the only
#' counterfactual, \code{B} will stay at its actual level.
#' We also take special pains to "remember" information about actual and 
#' imputed levels of counterfactuals so that appropriate results are obtained when
#' \code{emmeans} is applied to a previous \code{emmeans} result.
#' 
#' 
#' @section Offsets:
#' Unlike in \code{ref_grid}, an offset need not be scalar. If not enough values
#' are supplied, they are cyclically recycled. For a vector of offsets, it is 
#' important to understand that the ordering of results goes with the first 
#' name in \code{specs} varying fastest. If there are any \code{by} factors,
#' those vary slower than all the primary ones, but the first \code{by} variable
#' varies the fastest within that hierarchy. See the examples.
#' 
#' @section Options and \code{...}:
#'   Arguments that could go in \code{options} may instead be included in \code{...},
#'   typically, arguments such as \code{type}, \code{infer}, etc. that in essence
#'   are passed to \code{\link{summary.emmGrid}}. Arguments in both places are 
#'   overridden by the ones in \code{...}.
#'   
#'   There is a danger that \code{...} arguments could partially match those used
#'   by both \code{ref_grid} and \code{update.emmGrid}, creating a conflict.
#'   If these occur, usually they can be resolved by providing complete (or at least 
#'   longer) argument names; or by isolating non-\code{ref_grid} arguments in
#'   \code{options}; or by calling \code{ref_grid} separately and passing the
#'   result as \code{object}. See a not-run example below.
#'   
#'   Also, when \code{specs} is a two-sided formula, or \code{contr} is specified,
#'   there is potential confusion concerning which \code{...} arguments
#'   apply to the means, and which to the contrasts. When such confusion is possible,
#'   we suggest doing things separately 
#'   (a call to \code{emmeans} with no contrasts, followed by a call to 
#'   \code{\link{contrast}}). We treat
#'   \code{adjust} as a special case: it is applied to the \code{emmeans} results 
#'   \emph{only} if there are
#'   no contrasts specified, otherwise it is passed only to \code{contrast}.

#'
#' @export
#' 
#' @seealso \code{\link{ref_grid}}, \code{\link{contrast}}, 
#' \href{../doc/models.html}{vignette("models", "emmeans")}
#'
#' @examples
#' warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
#' emmeans (warp.lm,  ~ wool | tension)
#' # or equivalently emmeans(warp.lm, "wool", by = "tension")
#' 
#' # 'adjust' argument ignored in emmeans, passed to contrast part...
#' emmeans (warp.lm, poly ~ tension | wool, adjust = "sidak")
#' 
#' # 'adjust' argument NOT ignored ...
#' emmeans (warp.lm, ~ tension | wool, adjust = "sidak")
#' 
#' # Get all sets of EMMs for this model
#' ( allsets <- emmeans(warp.lm, ".") )
#' contrast(allsets, "eff")    # all effects
#' 
#' 
#' \dontrun{
#'   ### Offsets: Consider a silly example:
#'   emmeans(warp.lm, ~ tension | wool, offset = c(17, 23, 47)) @ grid
#'   # note that offsets are recycled so that each level of tension receives
#'   # the same offset for each wool.
#'   # But using the same offsets with ~ wool | tension will probably not
#'   # be what you want because the ordering of combinations is different.
#' }
emmeans = function(object, specs, by = NULL, 
                   fac.reduce = function(coefs) apply(coefs, 2, mean), 
                   contr, options = get_emm_option("emmeans"), 
                   weights, offset, ..., tran) {
    
    object = .chk.list(object, ...)
    if(!is(object, "emmGrid")) {
        args = .zap.args(object = object, ..., omit = "submodel")
        if (is.null(args$wt.nuis)) # pass weights as wt.nuis
            args$wt.nuis = ifelse(!missing(weights) && is.character(weights), weights, "equal")
        object = do.call(ref_grid, args)
    }
    
    # Check to see if we want all sets of means
    specs = .parse.specs.for.all(object, specs, by)
    if(!is.null(rtn <- attr(specs, "form.rtn"))) {
        if(length(rtn$by) > 0)
            by = rtn$by
        if(length(rtn$lhs) > 0) 
            contr = rtn$lhs
    }
    
    if (is.list(specs)) {
        return (emmeans.list(object, specs, by = by, 
                             contr = contr, weights = weights, 
                             offset = offset, ...))
    }
    if (inherits(specs, "formula")) {
        spc = .parse.by.formula(specs)
        specs = spc$rhs
        if (length(spc$by) > 0) 
            by = setdiff(union(spc$by, by), spc$rhs)
        if (length(spc$lhs) > 0) 
            contr = spc$lhs
    }
    
    if (!missing(tran)) {
        options $tran = tran
    }
    
    if(is.null(nesting <- object@model.info$nesting)) 
        {
        RG = object
        facs = union(specs, by)
        
        # Check that grid is complete
        # This isn't a 100% reliable check, but...
        if(nrow(RG@grid) != prod(sapply(RG@levels, length)))
            stop("Irregular reference grid: Marginal means cannot be determined.\n",
                 "You can possibly fix this with the 'force_regular' function.")
        
        if (!is.null(RG@misc$display)) {
            RG@misc$display = NULL
            warning("emmeans() results may be corrupted by removal of a nesting structure")
        }
        
        ## Handle counterfactuals...
        # cf.grid is either NULL, logical, or a "parent" grid that should replace RG
        if(!is.null(cf.grid <- RG@misc$cf.grid)) {
            weights = "cells"  # always use cells weights with counterfactuals
            if (inherits(cf.grid, "emmGrid"))
                RG = cf.grid
            # which counterfactuals do we average over?
            cf.ao = intersect(setdiff(names(RG@levels), facs), RG@roles$counterfactuals)
            # zero-out weights for any cases where cf.ao levels differ from actual levels
            for(f in cf.ao) {
                fc = paste0("actual_", f)
                excl = (RG@grid[ , f] != RG@grid[, fc])
                RG@grid[excl, ".wgt."] = 0
            }
            # Fix up the returned grid to play along
            cf.surv = intersect(facs, RG@roles$counterfactuals)
            if(length(cf.surv) == 0) # no cfs remain
                RG@misc$cf.grid = RG@roles$counterfactuals = NULL
            else {
                if(length(intersect(facs, paste0("actual_", cf.surv))) == length(cf.surv))
                    RG@misc$cf.grid = TRUE
                else
                    RG@misc$cf.grid = RG   # save this as "parent" grid
            }
        }
        
        # Ensure object is in standard order
        ord = .std.order(RG@grid, RG@levels) ###do.call(order, unname(RG@grid[rev(names(RG@levels))]))
        if(any(ord != seq_along(ord)))
            RG = RG[ord]
        
        # xxx if ((length(facs) == 1) && (facs == "1")) {  ### just want grand mean
        if("1" %in% facs) {
            RG@levels[["1"]] = "overall"
            RG@grid[ ,"1"] = 1
        }

        
        # Figure out the structure of the grid
        wgt = RG@grid[[".wgt."]]
        if(!is.null(wgt) && all(zapsmall(wgt) == 0)) wgt = wgt + 1 ### repl all zero wgts with 1
        dims = sapply(RG@levels, length)
        row.idx = array(seq_len(nrow(RG@linfct)), dims)
        use.mars = match(facs, names(RG@levels)) # which margins to use
        avgd.mars = setdiff(seq_along(dims)[dims>1], use.mars) # margins that we average over
        
        # Reconcile weights, if there are any margins left
        if ((length(avgd.mars) > 0) && !missing(weights)) {
            if (is.character(weights)) {
                if (is.null(wgt))
                    warning("'weights' requested but no weighting information is available")
                else {
                    wopts = c("equal","proportional","outer","cells","flat","show.levels","invalid")
                    weights = switch(wopts[pmatch(weights, wopts, 7)],
                                     equal = rep(1, prod(dims[avgd.mars])),
                                     proportional = as.numeric(apply(row.idx, avgd.mars,
                                                                           function(idx) sum(wgt[idx]))),
                                     outer = {
                                         ftbl = apply(row.idx, avgd.mars,
                                                            function(idx) sum(wgt[idx]))
                                         # Fix up the dimensions
                                         ftbl = array(ftbl, dim(row.idx)[avgd.mars])
                                         w = N = sum(ftbl)
                                         for (d in seq_along(dim(ftbl)))
                                             w = outer(w, apply(ftbl, d, sum) / N)
                                         as.numeric(w)
                                     },
                                     cells = "fq",
                                     flat = "fl",
                                     show.levels = {
                                         cat("emmeans are obtained by averaging over these factor combinations\n")
                                         return(do.call(expand.grid, RG@levels[avgd.mars]))
                                     },
                                     invalid = stop("Invalid 'weights' option: '", weights, "'")
                    )
                }
            }
            if (is.matrix(weights)) {
                wtrow = 0
                fac.reduce = function(coefs) {
                    wtvec = weights[wtrow+1, ] / sum(weights[wtrow+1, ])
                    ans = apply(sweep(coefs, 1, wtvec, "*"), 2, sum)
                    wtrow <<- (1 + wtrow) %% nrow(weights)
                    ans
                }
            }
            else if (is.numeric(weights)) {
                wtsum = sum(weights)
                if (wtsum <= 1e-8) wtsum = NA
                fac.reduce = function(coefs) {
                    if (nrow(coefs) != length(weights))
                        stop("Nonconforming number of weights -- need ", nrow(coefs))
                    apply(sweep(coefs, 1, weights, "*"), 2, sum) / wtsum
                }
            }
        }
        
        # Get the required factor combs
        levs = list()
        for (f in facs) {
            levs[[f]] = RG@levels[[f]]
            if (!hasName(levs, f))
                stop(paste("No variable named", f, "in the reference grid"))
        }
        combs = do.call("expand.grid", levs)
        if (!missing(weights) && is.character(weights) && (weights %in% c("fq", "fl")))
            K = apply(row.idx, use.mars, function(idx) {
                fq = RG@grid[[".wgt."]][idx]
                if (weights == "fl")
                    fq = 0 + (fq > 0)  # fq = 1 if > 0, else 0
                apply(sweep(RG@linfct[idx, , drop=FALSE], 1, fq, "*"), 2, sum) / sum(fq)
            })
        else
            K = apply(row.idx, use.mars, function(idx) {
                fac.reduce(RG@linfct[idx, , drop=FALSE])
            })
        
        linfct = t(matrix(K, nrow = ncol(RG@linfct),
                          dimnames = list(colnames(RG@linfct), NULL)))

        if(.some.term.contains(union(facs, RG@roles$trend), RG@model.info$terms))
            if(get_emm_option("msg.interaction"))
                message("NOTE: Results may be misleading due ", 
                        "to involvement in interactions")
        
        # Figure offset, if any
        if (hasName(RG@grid, ".offset.")) {
            combs[[".offset."]] = as.numeric(apply(row.idx, use.mars, function(idx)
                fac.reduce(as.matrix(RG@grid[idx, ".offset.", drop = FALSE]))))
        }
        
        avgd.over = names(RG@levels[avgd.mars])
        
        # add/override .offset. column if requested
        if(!missing(offset)) {
            combs[[".offset."]] = rep(offset, nrow(combs))[seq_len(nrow(combs))]
            
        }
        
        # Update .wgt column of grid, if it exists
        if (!is.null(wgt)) {
            combs[[".wgt."]] = as.numeric(apply(row.idx, use.mars, 
                                                      function(idx) sum(wgt[idx])))
        }
        
        RG@roles$responses = character()
        RG@misc$is.new.rg = NULL
        RG@misc$famSize = nrow(linfct)
        if(RG@misc$estName == "prediction") 
            RG@misc$estName = "emmean"
        RG@misc$adjust = "none"
        RG@misc$infer = c(TRUE,FALSE)
        RG@misc$pri.vars = setdiff(facs, by)
        RG@misc$by.vars = by
        RG@misc$avgd.over = union(RG@misc$avgd.over, avgd.over)
        RG@misc$methDesc = "emmeans"
        RG@roles$predictors = setdiff(names(levs), RG@roles$multresp)
        # if ((length(RG@roles$multresp) > 0) && !(RG@roles$multresp %in% names(levs))) 
        #     RG@roles$multresp = character(0)
        # REPLACED BY:
        RG@roles$multresp = intersect(RG@roles$multresp, names(levs))
        
        result = as.emmGrid(RG)
        result@linfct = linfct
        result@levels = levs
        result@grid = combs
        
        result = if (missing(contr))
            .update.options(result, options, ...)
        else
            .update.options(result, options, ..., exclude.opts = "adjust")
    }
    
    else {  # handle a nested structure
        object@model.info$nesting = NULL
        result = .nested_emm(object, specs, by = by, fac.reduce = fac.reduce, 
                       options = options, weights = weights, offset = offset, nesting = nesting)
        if(!is.null(type <- list(...)$type))
            result = update(result, type = type)
    }
    

    if(!missing(contr)) { # return a list with emmeans and contrasts
        # # skip doing contrasts when we have no primary factor
        # if(length(result@misc$pri.vars) == 0) 
        #     return(result)
        warn.contr = interactive() && (sys.parent() == 0) && is.character(contr) && 
            (length(result@misc$pri.vars) > 1)
        # NULL-out a bunch of arguments to not pass. 
        dontpass = c("data", "avgd.over", "by.vars", "df", "initMesg", "estName", "estType",
                     "famSize", "inv.lbl", "methDesc", "nesting", "pri.vars", 
                     "tran", "tran.mult", "tran.offset", "tran2", "is.new.rg")
        args = .zap.args(object = result, method = contr, by = by, ..., omit = dontpass)
        ctrs = do.call(contrast, args)
        result = .cls.list("emm_list", emmeans = result, contrasts = ctrs)
        if(!is.null(lbl <- object@misc$methDesc))
            names(result)[1] = lbl
        
        if (warn.contr && (nrow(result$contr@grid) > 21))
            warning("You may have generated more contrasts than you really wanted. In the future,\n",
                    "we suggest you avoid things like 'pairwise ~ fac1*fac2' when you have\n",
                    "more than one factor. Instead, call emmeans() with just '~ fac1*fac2' and do the\n",
                    "contrasts you need in a later step. See vignette(\"QuickStart\", \"emmeans\").\n",
                     call. = FALSE)
    }

    result
}



# Construct a new emmGrid object with given arguments

#' Construct an \code{emmGrid} object from scratch
#' 
#' This allows the user to incorporate results obtained by some analysis
#' into an \code{emmGrid} object, enabling the use of \code{emmGrid} methods
#' to perform related follow-up analyses.
#' 
#'  The arguments must be conformable. This includes that the length of
#'  \code{bhat}, the number of columns of \code{linfct}, and the number of
#'  columns of \code{post.beta} must all be equal. And that the product of
#'  lengths in \code{levels} must be equal to the number of rows of
#'  \code{linfct}. The \code{grid} slot of the returned object is generated 
#'  by \code{\link{expand.grid}} using \code{levels} as its arguments. So the
#'  rows of \code{linfct} should be in corresponding order.
#'  
#' The functions \code{qdrg} and \code{\link{emmobj}} are close cousins, in that
#' they both produce \code{emmGrid} objects. When starting with summary
#' statistics for an existing grid, \code{emmobj} is more useful, while
#' \code{qdrg} is more useful when starting from an unsupported fitted model.
#'
#'
#' @param bhat Numeric. Vector of regression coefficients
#' @param V Square matrix. Covariance matrix of \code{bhat}
#' @param levels Named list or vector. Levels of factor(s) that define the
#'   estimates defined by \code{linfct}. If not a list, we assume one factor
#'   named \code{"level"}
#' @param linfct Matrix. Linear functions of \code{bhat} for each combination 
#'   of \code{levels}. 
#' @param df Numeric value or function with arguments \code{(x, dfargs)}. If a
#'   number, that is used for the degrees of freedom. If a function, it should
#'   return the degrees of freedom for \code{sum(x*bhat)}, with any additional
#'   parameters in \code{dfargs}.
#' @param dffun Overrides \code{df} if specified. This is a convenience
#'   to match the slot names of the returned object.
#' @param dfargs List containing arguments for \code{df}.
#'   This is ignored if df is numeric.
#' @param post.beta Matrix whose columns comprise a sample from the posterior
#'   distribution of the regression coefficients (so that typically, the column
#'   averages will be \code{bhat}). A 1 x 1 matrix of \code{NA} indicates that
#'   such a sample is unavailable.
#' @param nesting Nesting specification as in \code{\link{ref_grid}}. This is
#'   ignored if \code{model.info} is supplied.
#' @param se.bhat,se.diff Alternative way of specifying \code{V}. If \code{se.bhat}
#'   is specified, \code{V} is constructed using \code{se.bhat}, the standard errors of \code{bhat},
#'   and \code{se.diffs}, the standard errors of its pairwise differences. \code{se.diff}
#'   should be a vector of length \code{k*(k-1)/2} where \code{k} is the length
#'   of \code{se.bhat}, and its elements should be in the order \code{12,13,...,1k,23,...2k,...}.
#'   If \code{se.diff} is missing, \code{V} is computed as if the \code{bhat} are independent.
#' @param ... Arguments passed to \code{\link{update.emmGrid}}
#' 
#' @seealso \code{\link{qdrg}}, an alternative that is useful when starting 
#'   with a fitted model not supported in \pkg{emmeans}.
#'
#' @return An \code{emmGrid} object
#' @export
#'
#' @examples
#' # Given summary statistics for 4 cells in a 2 x 2 layout, obtain 
#' # marginal means and comparisons thereof. Assume heteroscedasticity
#' # and use the Satterthwaite method
#' levels <- list(trt = c("A", "B"), dose = c("high", "low"))
#' ybar <- c(57.6, 43.2, 88.9, 69.8)
#' s <-    c(12.1, 19.5, 22.8, 43.2)
#' n <-    c(44,   11,   37,   24)
#' se2 = s^2 / n
#' Satt.df <- function(x, dfargs)
#'     sum(x * dfargs$v)^2 / sum((x * dfargs$v)^2 / (dfargs$n - 1))
#'     
#' expt.emm <- emmobj(bhat = ybar, V = diag(se2),
#'     levels = levels, linfct = diag(c(1, 1, 1, 1)),
#'     df = Satt.df, dfargs = list(v = se2, n = n), estName = "mean")
#' 
#' ( trt.emm <- emmeans(expt.emm, "trt") )
#' ( dose.emm <- emmeans(expt.emm, "dose") )
#' 
#' rbind(pairs(trt.emm), pairs(dose.emm), adjust = "mvt")
#' 
#' ### Create an emmobj from means and SEs
#' ### (This illustration reproduces the MOats example for Variety = "Victory")
#' means = c(71.50000,  89.66667, 110.83333, 118.50000)
#' semeans = c(5.540591, 6.602048, 8.695358, 7.303221)
#' sediffs = c(7.310571,  9.894724,  7.463615, 10.248306,  4.935698,  8.694507)
#' foo = emmobj(bhat = means, se.bhat = semeans, se.diff = sediffs, 
#'              levels = list(nitro = seq(0, .6, by = .2)), df = 10)
#' plot(foo, comparisons = TRUE)
#' 
emmobj = function(bhat, V, levels, linfct = diag(length(bhat)), df = NA, dffun, dfargs = list(), 
                  post.beta = matrix(NA), nesting = NULL, 
                  se.bhat, se.diff, ...) {
    if (!missing(se.bhat)) {
        V = D = diag(vb <- se.bhat^2)
        if (!missing(se.diff)) {
            V[lower.tri(V)] = se.diff^2
            V = (-V - t(V) + sweep(D, 1, vb, "+") + + sweep(D, 2, vb, "+"))/2
            diag(V) = vb
        }
    }
    if ((nrow(V) != ncol(V)) || (nrow(V) != ncol(linfct)) || (length(bhat) != ncol(linfct)))
        stop("bhat, V, and linfct are incompatible")
    if (!is.list(levels))
        levels = list(level = levels)
    grid = do.call(expand.grid, levels)
    if (nrow(grid) != nrow(linfct))
        stop("linfct should have ", nrow(grid), "rows")
    pri.vars = names(grid)
    dotargs = list(...)
    model.info = dotargs$model.info
    if(is.null(model.info))
        model.info = list(call = str2lang("emmobj"), xlev = levels, 
                          nesting = .parse_nest(nesting))
    roles = list(predictors= names(grid), responses=character(0), 
                 multresp=character(0))
    for (nm in names(dotargs$extras))
        grid[[nm]] = dotargs$extras[[nm]]
    if (!missing(dffun))
        df = dffun
    if (is.function(df)) {
        dffun = df
    } 
    else {
        dffun = function(x, dfargs) dfargs$df
        dfargs = list(df = df)
    }
    ### let's not be so brute-force about what's needed in misc
    # misc = list(estName = "estimate", estType = "prediction", infer = c(TRUE,FALSE), level = .95,
    #             adjust = "none", famSize = nrow(linfct), 
    #             avgd.over = character(0), pri.vars = pri.vars,
    #             methDesc = "emmobj", display = dotargs$display, .pairby = dotargs$.pairby)
    dflt = list(estName = "estimate", estType = "prediction", infer = c(TRUE,FALSE), level = .95,
                adjust = "none", famSize = nrow(linfct),
                avgd.over = character(0), pri.vars = pri.vars,
                methDesc = "emmobj")
    misc = dflt[setdiff(names(dflt), names(dotargs))]  # use defaults when not specified
    nm = setdiff(names(dotargs), c(names(misc), "extras", "model.info")) # otherwise use dots
    misc = c(misc, dotargs[nm])
    
    result = new("emmGrid", model.info=model.info, roles=roles, grid=grid,
                 levels = levels, matlevs=list(),
                 linfct=linfct, bhat=bhat, nbasis=all.estble, V=V,
                 dffun=dffun, dfargs=dfargs, misc=misc, post.beta=post.beta)
    dotargs$model.info = dotargs$extras = dotargs$display = NULL
    do.call(update, c(object = result, dotargs, silent = TRUE))
}

#' Convert to and from \code{emmGrid} objects
#' 
#' These are useful utility functions for creating a compact version of an
#' \code{emmGrid} object that may be saved and later reconstructed, or for
#' converting old \code{ref.grid} or \code{lsmobj} objects into \code{emmGrid}
#' objects.
#' 
#' An \code{emmGrid} object is an S4 object, and as such cannot be saved in a
#' text format or saved without a lot of overhead. By using \code{as.list},
#' the essential parts of the object are converted to a list format that can be
#' easily and compactly saved for use, say, in another session or by another user.
#' Providing this list as the arguments for \code{\link{emmobj}} allows the user 
#' to restore a working \code{emmGrid} object.
#' 
#' @param object Object to be converted to class \code{emmGrid}. It may
#'   be a \code{list} returned by \code{as.list.emmGrid}, or a \code{ref.grid}
#'   or \code{lsmobj} object created by \pkg{emmeans}'s predecessor, the 
#'   \pkg{lsmeans} package. An error is thrown if \code{object} cannot
#'   be converted.
#' @param ... In \code{as.emmGrid}, additional arguments passed to 
#'   \code{\link{update.emmGrid}} before returning the object. This
#'   argument is ignored in \code{as.list.emmGrid}
#'   
#' @return \code{as.emmGrid} returns an object of class \code{emmGrid}. 
#'     However, in fact, both \code{as.emmGrid} and \code{as.emm_list} check for an
#'     attribute in \code{object} to decide whether to return an \code{emmGrid} 
#'     or \code{emm_list)} object.
#' 
#' @seealso \code{\link{emmobj}}
#' @export
#' 
#' @examples
#' pigs.lm <- lm(log(conc) ~ source + factor(percent), data = pigs)
#' pigs.sav <- as.list(ref_grid(pigs.lm))
#' 
#' pigs.anew <- as.emmGrid(pigs.sav)
#' emmeans(pigs.anew, "source")
as.emmGrid = function(object, ...) {
    if (cls <- class(object)[1] %in% c("ref.grid", "lsmobj")) {
        object = as.list.emmGrid(object)
        if (is.null(object$misc$is.new.rg))
            object$misc$is.new.rg = (cls == "ref.grid")
    }
    # above keeps us from having to define these classes in emmeans
    if (is.list(object)) {
        if (!is.null(attr(object, "emm_list")))
            return(as.emm_list(object))
        else
            result = do.call(emmobj, object)
    }
    else {
        result = try(as(object, "emmGrid", strict = FALSE), silent = TRUE)
        if (inherits(result, "try-error"))
            stop("Object cannot be coerced to class 'emmGrid'")
    }
    update(result, ...)
}

#' @rdname as.emmGrid
#' @order 2
#' @param x An \code{emmGrid} object
#' @param model.info.slot Logical value: Include the \code{model.info} slot?
#'   Set this to \code{TRUE} if you want to preserve the original call and 
#'   information needed by the \code{submodel} option.
#'   If \code{FALSE}, only the nesting information (if any) is saved
#' @return \code{as.list.emmGrid} returns an object of class \code{list}.
#' @method as.list emmGrid
#' @export
as.list.emmGrid = function(x, model.info.slot = FALSE, ...) {
    slots = c("bhat", "V", "levels", "linfct", "dffun", "dfargs", "post.beta")
    result = lapply(slots,function(nm) slot(x, nm))
    names(result) = slots
    result = c(result, x@misc)
    if(model.info.slot)
        result$model.info = x@model.info
    else
        result$nesting = x@model.info$nesting
    nm = intersect(names(x@grid), c(".wgt.", ".offset."))
    if (length(nm) > 0)
        result$extras = x@grid[nm]
    result$pri.vars = NULL
    result
}



#### --- internal stuff used only by emmeans -------------

# Check if model contains a term containing all elts of facs
# terms is terms() component of model
.some.term.contains = function(facs, terms) {
    for (trm in attr(terms, "term.labels")) {
        if(all(sapply(facs, function(f) length(grep(f,trm))>0)))
            if (length(.all.vars(as.formula(paste("~",trm)))) > length(facs)) 
                return(TRUE)
    }
    return(FALSE)
}

### Sort grid in standard order according to ordering of entries in levels. 
### Thus .std.order(do.call(expand.grid, levels), levels) --> 1,2,...,nrow
.std.order = function(grid, levels) {
    tmp = lapply(rev(names(levels)), function(nm) {
        x = grid[[nm]]
        if (inherits(x, "factor")) as.integer(x)
        else as.integer(factor(x, levels = as.character(levels[[nm]])))
        # Note: need as.character(levels) here so we handle such as Date vectors correctly
    })
    do.call(order, tmp)
}

