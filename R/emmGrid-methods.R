##############################################################################
#    Copyright (c) 2012-2022 Russell V. Lenth                                #
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


### =========== Various methods for emmGrid class =============================
# (note: some major ones have their own file)



### S4 show method
## use S3 for this setMethod("summary", "emmGrid", summary.emmGrid)
setMethod("show", "emmGrid", function(object) {
    isnewrg = object@misc$is.new.rg
    if (is.null(isnewrg)) 
        isnewrg = FALSE
    
    if (isnewrg)
        str.emmGrid(object)
    else
        print(summary.emmGrid(object))
})


### Others are all S3 methods

#' @rdname emmGrid-methods
#' @method str emmGrid
#' @export
str.emmGrid <- function(object, ...) {
    showlevs = function(x) { # internal convenience function
        if (is.null(x)) cat("(predicted by other variables)")
        else {
            if (length(x) > 9) x = c(x[1:9], "...", x[length(x)])
            cat(paste(format(x, digits = 5, justify = "none"), collapse=", "))
        }
    }
    showtran = function(misc, label) { # internal convenience fcn
        cat(paste(label, dQuote(.fmt.tran(misc)), "\n"))
    }
    levs = object@levels
    cat(paste("'", class(object)[1], "' object with variables:\n", sep=""))
    for (nm in union(object@roles$predictors, union(object@roles$multresp, object@roles$responses))) {
        cat(paste("    ", nm, " = ", sep = ""))
        if (hasName(object@matlevs, nm)) {
            if (nm %in% object@roles$responses)
                cat("multivariate response with means: ")
            else
                cat("matrix with column reference values: ")
            cat("\n        ")
            showlevs(object@matlevs[[nm]])
        }
        else if (nm %in% object@roles$multresp) {
            cat("multivariate response levels: ")
            showlevs(levs[[nm]])
        }
        else if (nm %in% object@roles$responses) {
            cat("response variable with mean ")
            showlevs(levs[[nm]])
        }
        else
            showlevs(levs[[nm]])
        cat("\n")
    }
    if (length(nuis <- object@roles$nuisance) > 0) {
        cat("Nuisance factors that have been collapsed by averaging:\n    ")
        tmp = paste0(names(nuis), "(", sapply(nuis, length), ")")
        cat(paste(tmp, collapse = ", "))
        cat("\n")
    }
    if (length(cf <- object@roles$counterfactuals) > 0) {
        cat("Counterfactuals\n    ")
        cat(paste(cf, collapse = ", "))
        cat("\n")
    }
    if(!is.null(object@model.info$nesting)) {
        cat("Nesting structure:  ")
        cat(.fmt.nest(object@model.info$nesting))
        cat("\n")
    }
    if(!is.null(tran <- object@misc$tran)) {
        showtran(object@misc, "Transformation:")
        if (!is.null(tran2 <- object@misc$tran2))
            showtran(list(tran = tran2), "Additional response transformation:")
    }
    if(!is.na(object@nbasis[1])) {
        cat(paste0("Some things are non-estimable (null space dim = ",
                  ncol(object@nbasis), ")\n"))
    }
}


#' @rdname emmGrid-methods
#' @param export Logical value. If \code{FALSE}, the object is printed.
#'   If \code{TRUE}, a list is invisibly returned, which contains character
#'   elements named \code{summary} and \code{annotations} that may be saved
#'   or displayed as the user sees fit. \code{summary} is a character matrix
#'   (or list of such matrices, if a \code{by} variable is in effect).
#'   \code{annotations} is a character vector of the annotations that would
#'   have been printed below the summary or summaries.
#' @method print emmGrid
#' @param x An \code{emmGrid} object
#' @export
print.emmGrid = function(x, ..., export = FALSE) {
    print(summary.emmGrid(x, ...), export = export)
    # if(export)
    #     print(summary.emmGrid(x, ...), export = export)
    # else
    #     show(x)
}

### Former print method which I changed to work around a bug in pkgdown
# print.emmGrid = function(x, ..., export = FALSE)
#     print(summary.emmGrid(x, ...), export = export)

# vcov method
#' Miscellaneous methods for \code{emmGrid} objects
#' @rdname emmGrid-methods
#' 
#' @param object An \code{emmGrid} object
#' @param ... (required but not used)
#' @param sep separator for pasting levels in creating row and column
#'        names for \code{vcov()} results
#' 
#' @return The \code{vcov} method returns a symmetric matrix of variances and
#'   covariances for \code{predict.emmGrid(object, type = "lp")}
#'
#' @method vcov emmGrid
#' @export
#' @examples
#' warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
#' warp.emm <- emmeans(warp.lm, ~ tension | wool)
#' vcov(warp.emm) |> zapsmall()
#' 
#' vcov(pairs(warp.emm), sep = "|") |> zapsmall()
vcov.emmGrid = function(object, ..., sep = get_emm_option("sep")) {
    tol = get_emm_option("estble.tol")
    if (!is.null(hook <- object@misc$vcovHook)) {
        if (is.character(hook)) 
            hook = get(hook)
        hook(object, tol = tol, ...)
    }
    else {
        if(is.null(disp <- object@misc$display))
            disp = seq_len(nrow(object@linfct))
        X = object@linfct[disp, , drop = FALSE]
        estble = estimability::is.estble(X, object@nbasis, tol) 
        X[!estble, ] = NA
        X = X[, !is.na(object@bhat), drop = FALSE]
        rtn = X %*% tcrossprod(object@V, X)
        largs = as.list(object@grid[disp, seq_along(object@levels), drop = FALSE])
        largs$sep = sep
        rownames(rtn) = colnames(rtn) = do.call(paste, largs)
        rtn
    }
}


#' @rdname emmGrid-methods
#'
#' @returns The \code{linfct} function and method returns the \code{linfct} slot of \code{object}.
#' @export
linfct = function(object, ...)
    UseMethod("linfct")

#' @rdname emmGrid-methods
#' @export
linfct.default = function(object, ...) {
    attr(object, "linfct")
}


# Method to alter contents of misc slot

#' Update an \code{emmGrid} object
#' 
#' Objects of class \code{emmGrid} contain several settings that affect such things as
#' what arguments to pass to \code{\link{summary.emmGrid}}. 
#' The \code{update} method allows safer management of these settings than
#' by direct modification of its slots.
#'
#' @param object An \code{emmGrid} object
#' @param ... Options to be set. These must match a list of known options (see
#'   Details)
#' @param silent Logical value. If \code{FALSE} (the default), a message is
#'   displayed if any options are not matched. If \code{TRUE}, no messages are
#'   shown.
#'
#' @return an updated \code{emmGrid} object.
#' 
#' @method update emmGrid
#' @order 1
#' @export
#' 
#' @section Details:
#' The names in \code{\dots} are partially matched against those that are valid, and if a match is found, it adds or replaces the current setting. The valid names are
#' 
#' \describe{
#' \item{\code{tran}, \code{tran2}}{(\code{list} or \code{character}) specifies
#' the transformation which, when inverted, determines the results displayed by
#' \code{\link{summary.emmGrid}}, \code{\link{predict.emmGrid}}, or \code{\link{emmip}} when
#' \code{type="response"}. The value may be the name of a standard
#' transformation from \code{\link{make.link}} or additional ones supported by
#' name, such as \code{"log2"}; or, for a custom transformation, a \code{list}
#' containing at least the functions \code{linkinv} (the inverse of the
#' transformation) and \code{mu.eta} (the derivative thereof). The
#' \code{\link{make.tran}} function returns such lists for a number of popular
#' transformations. See the help page of \code{\link{make.tran}} for details as
#' well as information on the additional named transformations that are
#' supported. \code{tran2} is just like \code{tran} except it is a second
#' transformation (i.e., a response transformation in a generalized linear
#' model).}
#' 
#' \item{\code{tran.mult}}{Multiple for \code{tran}. For example, for the
#' response transformation \samp{2*sqrt(y)} (or \samp{sqrt(y) + sqrt(y + 1)},
#' for that matter), we should have \code{tran = "sqrt"} and \code{tran.mult =
#' 2}. If absent, a multiple of 1 is assumed.}
#' 
#' \item{\code{tran.offset}}{Additive constant before a transformation is applied.
#' For example, a response transformation of \code{log(y + pi)} has
#' \code{tran.offset  = pi}. If no value is present, an offset of 0 is assumed.}
#' 
#' \item{\code{estName}}{(\code{character}) is the column label used for
#' displaying predictions or EMMs.}
#' 
#' \item{\code{inv.lbl}}{(\code{character)}) is the column label to use for
#' predictions or EMMs when \code{type="response"}.}
#' 
#' \item{\code{by.vars}}{(\code{character)} vector or \code{NULL}) the variables
#' used for grouping in the summary, and also for defining subfamilies in a call
#' to \code{\link{contrast}}.}
#' 
#' \item{\code{pri.vars}}{(\code{character} vector) are the names of the grid
#' variables that are not in \code{by.vars}. Thus, the combinations of their
#' levels are used as columns in each table produced by \code{\link{summary.emmGrid}}.}
#' 
#' \item{\code{alpha}}{(numeric) is the default significance level for tests, in
#' \code{\link{summary.emmGrid}} as well as \code{\link{plot.emmGrid}}
#' when \samp{CIs = TRUE}. Be cautious that methods that depend on
#' specifying \code{alpha} are prone to abuse. See the
#' discussion in \href{../doc/basics.html#pvalues}{\code{vignette("basics", "emmeans")}}.}
#' 
#' \item{\code{adjust}}{(\code{character)}) is the default for the \code{adjust}
#' argument in \code{\link{summary.emmGrid}}.}
#' 
#' \item{\code{cross.adjust}}{(\code{character)}) is the default for the \code{cross.adjust}
#' argument in \code{\link{summary.emmGrid}} (used for adjusting between groups).}
#' 
#' \item{\code{famSize}}{(integer) is the number of means involved in a family of
#' inferences; used in Tukey adjustment}
#' 
#' \item{\code{infer}}{(\code{logical} vector of length 2) is the default value
#' of \code{infer} in \code{\link{summary.emmGrid}}.}
#' 
#' \item{\code{level}}{(numeric) is the default confidence level, \code{level},
#' in \code{\link{summary.emmGrid}}. \emph{Note:} You must specify all five letters 
#' of \sQuote{level} to distinguish it from the slot name \sQuote{levels}.}
#' 
#' \item{\code{df}}{(numeric) overrides the default degrees of freedom with a
#' specified single value.}
#' 
#' \item{\code{calc}}{(list) additional calculated columns. See \code{\link{summary.emmGrid}}.}
#' 
#' \item{\code{null}}{(numeric) null hypothesis for \code{summary} or
#' \code{test} (taken to be zero if missing).}
#' 
#' \item{\code{side}}{(numeric or character) \code{side} specification for for
#' \code{summary} or \code{test} (taken to be zero if missing).}
#' 
#' \item{\code{sigma}}{(numeric) Error SD to use in predictions and for bias-adjusted
#' back-transformations}
#' 
#' \item{\code{delta}}{(numeric) \code{delta} specification for \code{summary}
#' or \code{test} (taken to be zero if missing).}
#' 
#' \item{\code{predict.type} or \code{type}}{(character) sets the default method
#' of displaying predictions in \code{\link{summary.emmGrid}},
#' \code{\link{predict.emmGrid}}, and \code{\link{emmip}}. Valid values are
#' \code{"link"} (with synonyms \code{"lp"} and \code{"linear"}), or
#' \code{"response"}.}
#' 
#' \item{\code{bias.adjust}, \code{frequentist}}{(logical) These
#' are used by \code{summary} if the value of these arguments are not specified.}
#' 
#' \item{\code{estType}}{(\code{character}) is used internally to determine 
#' what \code{adjust} methods are appropriate. It should match one of 
#' \samp{c("prediction", "contrast", "pairs")}. As an example of why this is needed,
#' the Tukey adjustment should only be used for pairwise comparisons 
#' (\code{estType = "pairs"}); if \code{estType} is some other string, Tukey
#' adjustments are not allowed.}
#' 
#' \item{\code{avgd.over}}{(\code{character)} vector) are the names of the 
#' variables whose levels are averaged over in obtaining marginal averages of 
#' predictions, i.e., estimated marginal means. Changing this might produce a 
#' misleading printout, but setting it to \code{character(0)} will suppress the 
#' \dQuote{averaged over} message in the summary.}
#' 
#' \item{\code{initMesg}}{(\code{character}) is a string that is added to the
#' beginning of any annotations that appear below the \code{\link{summary.emmGrid}}
#' display.}
#' 
#' \item{\code{methDesc}}{(\code{character}) is a string that may be used for
#' creating names for a list of \code{emmGrid} objects. }
#' 
#' \item{\code{nesting}}{(Character or named \code{list}) specifies the nesting
#' structure. See \dQuote{Recovering or overriding model information} in the
#' documentation for \code{\link{ref_grid}}. The current nesting structure is
#' displayed by \code{\link{str.emmGrid}}.}
#'
#' \item{\code{levels}}{named \code{list} of new levels for the elements of the
#' current \code{emmGrid}. The list name(s) are used as new variable names, and
#' if needed, the list is expanded using \code{expand.grid}. These results replace
#' current variable names and levels. This specification changes the \code{levels},
#' \code{grid}, \code{roles}, and \code{misc} slots in the updated \code{emmGrid},
#' and resets \code{pri.vars}, \code{by.vars}, \code{adjust}, \code{famSize},
#' and \code{avgd.over}. In addition, if there is nesting of factors, that may be 
#' altered; a warning is issued if it involves something other than mere name changes.
#' \emph{Note:} All six letters of \code{levels} is needed in order to distinguish
#' it from \code{level}.}
#' 
#' \item{\code{submodel}}{\code{formula} or \code{character} value specifying a 
#' submodel (requires this feature being supported by underlying methods 
#' for the model class). When specified, the \code{linfct} slot is replaced by 
#' its aliases for the specified sub-model. Any factors in the sub-model that 
#' do not appear in the model matrix are ignored, as are any interactions that 
#' are not in the main model, and any factors associate with multivariate responses. 
#' The estimates displayed are then computed as if 
#' the sub-model had been fitted. (However, the standard errors will be based on the
#' error variance(s) of the full model.) 
#' \emph{Note:} The formula should refer only to predictor names, \emph{excluding} any
#' function calls (such as \code{factor} or \code{poly}) that appear in the 
#' original model formula. See the example.
#' 
#' The character values allowed should partially 
#' match \code{"minimal"} or \code{"type2"}. With \code{"minimal"}, the sub-model
#' is taken to be the one only involving the surviving factors in \code{object}
#' (the ones averaged over being omitted). Specifying \code{"type2"} is the same as
#' \code{"minimal"} except only the highest-order term in the submodel is retained,
#' and all effects not containing it are orthogonalized-out. Thus, in a purely linear
#' situation such as an \code{lm} model, the joint test
#' of the modified object is in essence a type-2 test as in \code{car::Anova}.
#' 
#' Please note that it is possible (or even likely) that there will be disparity
#' between the \code{grid} and \code{linfct} slots when a submodel is used. This is
#' because \code{grid} contains the \emph{claimed} values of the predictors and
#' \code{linfct} contains \emph{aliases} of them computed from the submodel.
#' 
#' For some objects such as generalized linear models, specifying \code{submodel}
#' will typically not produce the same estimates or type-2 tests as would be
#' obtained by actually fitting a separate model with those specifications.
#' The reason is that those models are fitted by iterative-reweighting methods,
#' whereas the \code{submodel} calculations preserve the final weights used in
#' fitting the full model.}
#' 
#' \item{(any other slot name)}{If the name matches an element of
#' \code{slotNames(object)} other than \code{levels}, that slot is replaced by 
#' the supplied value, if it is of the required class (otherwise an error occurs). 
#' 
#' The user must be very careful in
#' replacing slots because they are interrelated; for example, the lengths
#' and dimensions of \code{grid}, \code{linfct}, \code{bhat}, and \code{V} must
#' conform.}
#' } %%%%%%% end \describe 
#' 
#' @note
#' When it makes sense, an option set by \code{update} will persist into 
#' future results based on that object. But some options are disabled as well.
#' For example, a \code{calc} option will be nulled-out if \code{contrast}
#' is called, because it probably will not make sense to do the same 
#' calculations on the contrast results, and in fact the variable(s) needed
#' may not even still exist.
#' \code{factor(percent)}.
#'
#' @seealso \code{\link{emm_options}}
#' @examples
#' # Using an already-transformed response:
#' pigs.lm <- lm(log(conc) ~ source * factor(percent), data = pigs)
#' 
#' # Reference grid that knows about the transformation
#' # and asks to include the sample size in any summaries:
#' pigs.rg <- update(ref_grid(pigs.lm), tran = "log", 
#'                     predict.type = "response",
#'                     calc = c(n = ~.wgt.))
#' emmeans(pigs.rg, "source")
#' 
#' # Obtain estimates for the additive model
#' # [Note that the submodel refers to 'percent', not 'factor(percent)']
#' emmeans(pigs.rg, "source", submodel = ~ source + percent)
#' 
#' # Type II ANOVA
#' joint_tests(pigs.rg, submodel = "type2")
update.emmGrid = function(object, ..., silent = FALSE) {
    args = list(...)
    # see .valid.misc below this function for list of legal options
    valid.slots = slotNames(object)
    valid.choices = union(.valid.misc, valid.slots)
    misc = object@misc
    for (nm in names(args)) {
        fullname = try(match.arg(nm, valid.choices), silent=TRUE)
        if(inherits(fullname, "try-error")) {
            if (!silent)
                message("Argument ", sQuote(nm), " was ignored. Valid choices are:\n",
                        paste(valid.choices, collapse=", "))
        }
        else {
            if (fullname == "type") fullname = "predict.type"
            if (fullname == "levels") {
                lvls = args[[nm]]
                if (!is.list(lvls))
                    stop("'levels' must be a named list.")
                nm = names(lvls)
                if (is.null(nm) || any(nm == ""))
                     stop("'levels' must be a named list.")
                grd = do.call(expand.grid, lvls)
                if (nrow(object@grid) != nrow(grd))
                    stop("Length of replacement levels does not match the number of rows in the grid")
                oldlvls = object@levels
                if(!is.null(object@model.info$nesting) &&
                   ((length(oldlvls) != length(lvls)) ||
                    any(sapply(oldlvls, length) != sapply(lvls, length))))
                    warning("Changes to levels may have altered nesting structure.\n",
                            "You likely need to also run 'update(..., nesting = ...)'")
                object@levels = lvls
                for (nm in c(".wgt.", ".offset."))
                        grd[[nm]] = object@grid[[nm]]
                object@grid = grd
                object@roles$predictors = misc$pri.vars = names(lvls)
                misc$by.vars = misc$avgd.over = NULL
                if (!is.null(object@model.info$nesting))
                    object@model.info$nesting = .find_nests(grd, NULL, FALSE, lvls)
                misc$adjust = "none"
                misc$famSize = nrow(grd)
            }
            else if((fullname == "estName") && 
                    !is.null(misc$predict.type) && startsWith(misc$predict.type, "r"))
                misc$inv.lbl = args[[nm]]
            else if (fullname %in% valid.slots) # all slots but "levels"
                slot(object, fullname) = args[[nm]]
            else {
                if (fullname == "by.vars") {
                    allvars = union(misc$pri.vars, misc$by.vars)
                    misc$pri.vars = setdiff(allvars, args[[nm]])
                }
                if (fullname == "pri.vars") {
                    allvars = union(misc$pri.vars, misc$by.vars)
                    misc$by.vars = setdiff(allvars, args[[nm]])
                }
                # special case - I keep nesting in model.info. Plus add'l checks
                if (fullname == "nesting") {
                    object@model.info$nesting = lst = .parse_nest(args[[nm]])
                    if(!is.null(lst)) {
                        nms = union(names(lst), unlist(lst))
                        if(!all(nms %in% names(object@grid)))
                            stop("Nonexistent variables specified in 'nesting'")
                        object@misc$display = .find.nonempty.nests(object, nms)
                    }
                }
                if (fullname == "submodel") {
                    if(!is.null(A <- .alias.matrix(object, args[[nm]]))) {
                        rcols = attr(A, "rcols")
                        L = object@linfct
                        k = ncol(L) / ncol(A)
                        if(abs(k - (k<-as.integer(k))) > .01)
                            stop("Incompatible columns in alias setup")
                        ixmat = matrix(seq_along(L[1,]), ncol = k)
                        for (j in seq_len(k)) {
                            rc = ixmat[rcols, j]
                            fc =ixmat[, j]
                            object@linfct[, fc] = L[, rc, drop = FALSE] %*% A
                        }
                        object@model.info$model.matrix = "" # silent message
                        misc$initMesg = c(misc$initMesg, paste("submodel: ~", attr(A, "submodstr")))
                    }
                }
                    
                else
                    misc[[fullname]] = args[[nm]]
            }
        }
    }
    object@misc = misc
    object
}

### List of valid strings to match in update() ###
.valid.misc = c("adjust","alpha","avgd.over","bias.adjust","by.vars","calc",
                "cross.adjust", "delta","df",
               "initMesg","estName","estType","famSize","frequentist","infer","inv.lbl",
               "level","methDesc","nesting","null","predict.type","pri.vars",
               "side","sigma","tran","tran.mult","tran.offset","tran2","type","is.new.rg",
               "submodel")


#' Set or change emmeans options
#'
#' Use \code{emm_options} to set or change various options that are used in
#' the \pkg{emmeans} package. These options are set separately for different contexts in
#' which \code{emmGrid} objects are created, in a named list of option lists.
#' 
#' \pkg{emmeans}'s options are stored as a list in the system option \code{"emmeans"}. 
#' Thus, \code{emm_options(foo = bar)} is the same as 
#' \code{options(emmeans = list(..., foo = bar))} where \code{...} represents any
#' previously existing options. The list \code{emm_defaults} contains the default
#' values in case the corresponding element of system option \code{emmeans} is \code{NULL}.
#' 
#' Currently, the following main list entries are supported:
#' \describe{
#' \item{\code{ref_grid}}{A named \code{list} of defaults for objects created by
#' \code{\link{ref_grid}}. This could affect other objects as well. For example,
#' if \code{emmeans} is called with a fitted model object, it calls
#' \code{ref_grid} and this option will affect the resulting \code{emmGrid}
#' object.}
#' \item{\code{emmeans}}{A named \code{list} of defaults for objects created by
#'   \code{\link{emmeans}} or \code{\link{emtrends}}.}
#' \item{\code{contrast}}{A named \code{list} of defaults for objects created by
#'   \code{\link{contrast.emmGrid}} or \code{\link{pairs.emmGrid}}.}
#' \item{\code{summary}}{A named \code{list} of defaults used by the methods
#'   \code{\link{summary.emmGrid}}, \code{\link{predict.emmGrid}}, \code{\link{test.emmGrid}},
#'   \code{\link{confint.emmGrid}}, and \code{\link{emmip}}. The only option that can
#'   affect the latter four is \code{"predict.type"}.}
#' \item{\code{allow.na.levs}}{A logical value that if \code{TRUE} (the default), allows
#'   \code{NA} to be among the levels of a factor. Older versions of \pkg{emmeans} did not
#'   allow this. So if problems come up (say in a messy dataset that includes incomplete cases), 
#'   try setting this to \code{FALSE}.}
#' \item{\code{sep}}{A character value to use as a separator in labeling factor combinations.
#'   Such labels are potentially used in several places such as \code{\link{contrast}} and 
#'   \code{\link{plot.emmGrid}} when combinations of factors are compared or plotted.
#'   The default is \code{" "}.}
#' \item{\code{parens}}{Character vector that determines which labels are parenthesized
#'   when they are contrasted. The first element is a regular expression, and the second and
#'   third elements are used as left and right parentheses. 
#'   See details for the \code{parens} argument in \code{\link{contrast}}. The default
#'   will parenthesize labels containing the four arithmetic operators, 
#'   using round parentheses.}
#' \item{\code{cov.keep}}{The default value of \code{cov.keep} in \code{\link{ref_grid}}.
#'   Defaults to \code{"2"}, i.e., two-level covariates are treated like factors.}
#' \item{\code{graphics.engine}}{A character value matching 
#'   \code{c("ggplot", "lattice")}, setting the default engine to use in
#'   \code{\link{emmip}} and \code{\link{plot.emmGrid}}.  Defaults to \code{"ggplot"}.}
#' \item{\code{gg.theme}}{An integer. Set it to 1 or 2 if you want the appearance
#'   of plots used in \pkg{emmeans} version 1.x.x or 2.x.x. Or set it to the actual theme
#'   that you want to use, e.g., \code{emm_options(gg.theme = ggplot2::theme_dark())}}
#' \item{\code{msg.interaction}}{A logical value controlling whether or not
#'   a message is displayed when \code{emmeans} averages over a factor involved
#'   in an interaction. It is probably not appropriate to do this, unless
#'   the interaction is weak. Defaults to \code{TRUE}.}
#' \item{\code{msg.nesting}}{A logical value controlling whether or not to
#'   display a message when a nesting structure is auto-detected. The existence
#'   of such a structure affects computations of EMMs. Sometimes, a nesting
#'   structure is falsely detected -- namely when a user has omitted some
#'   main effects but included them in interactions. This does not change the
#'   model fit, but it produces a different parameterization that is picked
#'   up when the reference grid is constructed. Defaults to \code{TRUE}.}
#' \item{\code{rg.limit}}{An integer value setting a limit on the number of rows
#'   in a newly constructed reference grid. This is checked based on the number of
#'   levels of the factors involved; but it excludes the levels of any multivariate
#'   responses because those are not yet known. The reference grid consists of all
#'   possible combinations of the predictors, and this can become huge if there are
#'   several factors. An error is thrown if this limit is exceeded. One can use the 
#'   \code{nuisance} argument of \code{\link{ref_grid}} to collapse on nuisance
#'   factors, thus making the grid smaller. Defaults to 10,000.}
#' \item{\code{simplify.names}}{A logical value controlling whether to
#'   simplify (when possible) names in the model formula that refer to datasets --
#'   for example, should we simplify a predictor name like \dQuote{\code{data$trt}}
#'   to just \dQuote{\code{trt}}? Defaults to \code{TRUE}.}
#' \item{\code{opt.digits}}{A logical value controlling the precision with which
#'   summaries are printed. If \code{TRUE} (default), the number of digits
#'   displayed is just enough to reasonably distinguish estimates from the ends
#'   of their confidence intervals; but always at least 3 digits. If
#'   \code{FALSE}, the system value \code{getOption("digits")} is used.}
#' \item{\code{pval.digits}}{An integer indicating the precision with which
#'   p-values are printed. Integers between \code{2} and \code{6} are acceptable; the default is 
#' \code{4}. Floating point numbers are truncated. If a number outside the accepted range is
#' used, the nearest acceptable integer will be used instead. Non-numeric arguments are ignored.}
#' \item{\code{back.bias.adj}}{A logical value controlling whether we 
#'   try to adjust bias when back-transforming. If \code{FALSE}, we use naive
#'   back transformation. If \code{TRUE} \emph{and \code{sigma} is available and valid}, a
#'   second-order adjustment is applied to estimate the mean on the response
#'   scale. A warning is issued if no valid \code{sigma} is available}
#' \item{\code{enable.submodel}}{A logical value. If \code{TRUE}, enables support 
#'   for selected model classes to implement the \code{submodel} option. If
#'   \code{FALSE}, this support is disabled. Setting this option to \code{FALSE}
#'   could save excess memory consumption.}
#'   
#' }%%% end describe{}
#' Some other options have more specific purposes:
#' \describe{
#' \item{\code{estble.tol}}{Tolerance for determining estimability in
#' rank-deficient cases. If absent, the value in \code{emm_defaults$estble.tol)}
#' is used.}
#' \item{\code{save.ref_grid}}{Logical value of \code{TRUE} if you wish the 
#' latest reference grid created to be saved in \code{.Last.ref_grid}.
#' The default is \code{FALSE}.}
#' \item{Options for \code{lme4::lmerMod} models}{Options \code{lmer.df},
#' \code{disable.pbkrtest}, \code{pbkrtest.limit}, \code{disable.lmerTest},
#' and \code{lmerTest.limit}
#' options affect how degrees of freedom are computed for \code{lmerMod} objects
#' produced by the \pkg{lme4} package). See that section of the "models" vignette
#' for details.}
#' \item{\code{post.ci.method}}{method for estimating confidence intervals for Bayesian models. 
#' Options are \code{"HPD"} (default) or \code{"quantile"}.} 
#' } %%%%%% end \describe
#'
#' @param ... Option names and values (see Details)
#' @param disable If non-missing, this will reset all options to their defaults 
#'   if \code{disable} tests \code{TRUE} (but first save them for possible later 
#'   restoration). Otherwise, all previously saved options
#'   are restored. This is important for bug reporting; please see the section below
#'   on reproducible bugs. When \code{disable} is specified, the other arguments are ignored.
#' 
#' @return \code{emm_options} returns the current options (same as the result 
#'   of \samp{getOption("emmeans")}) -- invisibly, unless called with no arguments.
#' 
#' @section Reproducible bugs:
#' Most options set display attributes and such that are not likely to be associated
#' with bugs in the code. However, some other options (e.g., \code{cov.keep})
#' are essentially configuration settings that may affect how/whether the code
#' runs, and the settings for these options may cause subtle effects that may be
#' hard to reproduce. Therefore, when sending a bug report, please create a reproducible
#' example and make sure the bug occurs with all options set at their defaults.
#' This is done by preceding it with  \code{emm_options(disable = TRUE)}. 
#' 
#' By the way, \code{disable} works like a stack (LIFO buffer), in that \code{disable = TRUE}
#' is equivalent to \code{emm_options(saved.opts = emm_options())} and 
#' \code{emm_options(disable = FALSE)} is equivalent to 
#' \code{options(emmeans = get_emm_option("saved.opts"))}. To completely erase
#' all options, use \code{options(emmeans = NULL)}
#' 
#' @section Startup options:
#' \pkg{emmeans}'s options are kept in the system options in a list names \code{emmeans}.
#' Thus, if the user wants certain options set on startup, this can be arranged by
#' setting the \code{emmeans} option in one's startup script. For example,
#' the user's \code{.Rprofile} file could contain:
#' \preformatted{    options(emmeans = list(summary = list(predict.type = "response"), 
#'                            lmer.df = "satterthwaite"))}
#' 
#' 
#' @seealso \code{\link{update.emmGrid}}
#' @export
#' @examples
#' \dontrun{
#' emm_options(ref_grid = list(level = .90),
#'             contrast = list(infer = c(TRUE,FALSE)),
#'             estble.tol = 1e-6)
#' # Sets default confidence level to .90 for objects created by ref.grid
#' # AS WELL AS emmeans called with a model object (since it creates a 
#' # reference grid). In addition, when we call 'contrast', 'pairs', etc.,
#' # confidence intervals rather than tests are displayed by default.
#' }
#' 
#' \dontrun{
#' emm_options(disable.pbkrtest = TRUE)
#' # This forces use of asymptotic methods for lmerMod objects.
#' # Set to FALSE or NULL to re-enable using pbkrtest.
#' }
#' 
#' # See tolerance being used for determining estimability
#' get_emm_option("estble.tol")
#' 
#' \dontrun{
#' # Set all options to their defaults
#' emm_options(disable = TRUE)
#' # ... and perhaps follow with code for a minimal reproducible bug,
#' #     which may include emm_options() clls if they are pertinent ...
#' 
#' # restore options that had existed previously
#' emm_options(disable = FALSE)
#' }
#'
emm_options = function(..., disable) {
    opts = getOption("emmeans", list())
    newopts = list(...)
    display = TRUE  # flag to display all options if ... is empty
    if (!missing(disable)) {
        if (disable)
            opts = list(saved.opts = opts)
        else if (!is.null(saved <- opts$saved.opts))
            opts = saved
        newopts = list()
        display = FALSE
    }
    for (nm in names(newopts))
        opts[[nm]] = newopts[[nm]]
    options(emmeans = opts)
    if (length(newopts) > 0)
        invisible(opts)
    else if (display) {
        opts = c(opts, emm_defaults)
        opts[sort(names(opts))]
    }
}

# equivalent of getOption()
#' @rdname emm_options
#' @param x Character value - the name of an option to be queried
#' @param default Value to return if \code{x} is not found
#' @return \code{get_emm_option} returns the currently stored option for \code{x}, 
#'   or its default value if not found.
#' @export
get_emm_option = function(x, default = emm_defaults[[x]]) {
    opts = getOption("emmeans", list())
    if(is.null(default) || hasName(opts, x))
        opts[[x]]
    else 
        default
}

#' @rdname emm_options
#' @param expr Expression to evaluate. If missing, the last element of \code{...} is used.
#' @return \code{with_emm_options()} temporarily sets the options in \code{...}, then 
#'   evaluates \code{try(expr)} and returns the result.
#' @export
#' @examples
#' 
#' # Temporarily alter the formatting of labels in an interaction contrast:
#' pigs.lm <- lm(inverse(conc) ~ source * factor(percent), data = pigs)
#' with_emm_options(parens = c("\\-", "<", ">"), sep = ": ",
#'     contrast(contrast(ref_grid(pigs.lm), "consec", by = "percent"),
#'         "consec", by = NULL))
#' 
with_emm_options = function(..., expr) {
    cl = match.call()
    if(missing(expr)) {
        expr = cl[[length(cl)]]
        cl = cl[-length(cl)]
    }
    cl[[1]] = as.name("emm_options")
    oldopts = getOption("emmeans")
    eval(cl, parent.frame())
    result = try(eval(expr, parent.frame()), silent = TRUE)
    options(emmeans = oldopts)
    result
}



### Exported defaults for certain options
#' @rdname emm_options
#' @export
emm_defaults = list (
    ref_grid = list(is.new.rg = TRUE, infer = c(FALSE, FALSE)),
    emmeans = list(infer = c(TRUE, FALSE)),
    contrast = list(infer = c(FALSE, TRUE)),
    save.ref_grid = FALSE,    # save new ref_grid in .Last.ref_grid
    allow.na.levs = TRUE,     # Do we allow factors to have NA as a level?
    cov.keep = "2",           # default for cov.keep arg in ref_grid
    sep = " ",                # separator for combining factor levels
    parens = c("-|\\+|\\/|\\*", "(", ")"), # patterns for what/how to parenthesize in contrast
    graphics.engine = "ggplot",  # default for emmip and plot.emmGrid
    gg.theme = 2,             # Use the ggplot theme used in this emmeans version 
###    msg.data.call = TRUE,     # message when there's a call in data or subset
    msg.interaction = TRUE,   # message about averaging w/ interactions
    msg.nesting = TRUE,       # message when nesting is detected
    estble.tol = 1e-8,        # tolerance for estimability checks
    simplify.names = TRUE,    # simplify names like data$x to just "x"
    back.bias.adj = FALSE,    # Try to bias-adjust back-transformations?
    opt.digits = TRUE,        # optimize displayed digits for estimates and CI's?
    pval.digits = 4L,          # how many digits to list for p-values
    enable.submodel = TRUE,   # enable saving extra info for submodel
    rg.limit = 10000,         # limit on number of rows in a reference grid
    lmer.df = "kenward-roger",  # Use Kenward-Roger for df
    disable.pbkrtest = FALSE, # whether to bypass pbkrtest routines for lmerMod
    pbkrtest.limit = 3000,    # limit on N for enabling K-R
    disable.lmerTest = FALSE, # whether to bypass lmerTest routines for lmerMod
    lmerTest.limit = 3000,     # limit on N for enabling Satterthwaite
    post.ci.method = "HPD"    #  set CI method for Bayesian methods

)


# override levels<- method
#' @rdname update.emmGrid
#' @order 5
#' @export
#' @param x an \code{emmGrid} object
#' @param value \code{list} or replacement levels. See the documentation for
#' \code{update.emmGrid} with the \code{levels} argument, 
#' as well as the section below on \dQuote{Replaciong levels}
#'   
#' @return \code{levels<-} replaces the levels of the object in-place.
#'   See the section on replacing levels for details.
#' @section Replacing levels:
#' The \code{levels<-} method uses \code{update.emmGrid} to replace the
#' levels of one or more factors. This method allows selectively replacing
#' the levels of just one factor (via subsetting operators), whereas 
#' \code{update(x, levels = list(...))} requires a list of \emph{all} factors
#' and their levels. If any factors are to be renamed, we must replace all
#' levels and include the new names in the replacements. See the examples.
#' 
#' @examples
#' 
#' ## Changing levels of one factor
#' newrg <- pigs.rg
#' levels(newrg)$source <- 1:3
#' newrg
#' 
#' ## Unraveling a previously standardized covariate
#' zd = scale(fiber$diameter)
#' fibz.lm <- lm(strength ~ machine * zd, data = fiber)
#' (fibz.rg <- ref_grid(fibz.lm, at = list(zd = -2:2)))   ### 2*SD range
#' lev <- levels(fibz.rg)
#' levels(fibz.rg) <- list (
#'     machine = lev$machine,
#'     diameter = with(attributes(zd), 
#'                     `scaled:center` + `scaled:scale` * lev$zd) )
#' fibz.rg
#' 
"levels<-.emmGrid" = function(x, value) {
    update.emmGrid(x, levels = value)
}


### Utility to change the internal structure of an emmGrid object
### Returned emmGrid object has linfct = I and bhat = estimates
### Primary reason to do this is with transform = TRUE, then can 
### work with linear functions of the transformed predictions

#' Reconstruct a reference grid with a new transformation or simulations
#' 
#' The typical use of this function is to cause EMMs to be computed on
#' a different scale, e.g., the back-transformed scale rather than the 
#' linear-predictor scale. In other words, if you want back-transformed 
#' results, do you want to average and then back-transform, or 
#' back-transform and then average?
#' 
#' The \code{regrid} function reparameterizes an existing \code{ref.grid} so
#' that its \code{linfct} slot is the identity matrix and its \code{bhat} slot
#' consists of the estimates at the grid points. If \code{transform} is
#' \code{TRUE}, the inverse transform is applied to the estimates. Outwardly,
#' when \code{transform = "response"}, the result of \code{\link{summary.emmGrid}}
#' after applying \code{regrid} is identical to the summary of the original
#' object using \samp{type="response"}. But subsequent EMMs or
#' contrasts will be conducted on the new scale -- which is
#' the reason this function exists. 
#' 
#' This function may also be used to simulate a sample of regression
#' coefficients for a frequentist model for subsequent use as though it were a
#' Bayesian model. To do so, specify a value for \code{N.sim} and a sample is
#' simulated using the function \code{sim}. The grid may be further processed in
#' accordance with the other arguments; or if \code{transform = "pass"}, it is
#' simply returned with the only change being the addition of the simulated
#' sample.
#' 
#' @param object An object of class \code{emmGrid}
#' @param transform Character, list, or logical value. If \code{"response"},
#'   \code{"mu"}, or \code{TRUE}, the inverse transformation is applied to the
#'   estimates in the grid (but if there is both a link function and a response
#'   transformation, \code{"mu"} back-transforms only the link part); if
#'   \code{"none"} or \code{FALSE}, \code{object} is re-gridded so that its
#'   \code{bhat} slot contains \code{predict(object)} and its \code{linfct} slot
#'   is the identity. Any internal transformation information is preserved. If
#'   \code{transform = "pass"}, the object is not re-gridded in any way (this
#'   may be useful in conjunction with \code{N.sim}).
#'   
#'   If \code{transform} is a character value in \code{links} (which is the set
#'   of valid arguments for the \code{\link{make.link}} function, excepting
#'   \code{"identity"}), or if \code{transform} is a list of the same form as
#'   returned by \code{make.links} or \code{\link{make.tran}}, the results are
#'   formulated as if the response had been transformed with that link function.
#'   
#' @param inv.link.lbl Character value. This applies only when \code{transform} 
#'   is in \code{links}, and is used to label the predictions if subsequently summarized
#'   with \code{type = "response"}.
#' @param predict.type Character value. If provided, the returned object is
#'   updated with the given type to use by default by \code{summary.emmGrid}
#'   (see \code{\link{update.emmGrid}}).  This may be useful if, for example,
#'   when one specifies \code{transform = "log"} but desires summaries to be
#'   produced by default on the response scale.
#' @param bias.adjust Logical value for whether to adjust for bias in
#'   back-transforming (\code{transform = "response"}). This requires a valid value of 
#'   \code{sigma} to exist in the object or be specified.
#' @param sigma Error SD assumed for bias correction (when 
#'   \code{transform = "response"} and a transformation
#'   is in effect). If not specified,
#'   \code{object@misc$sigma} is used, and a warning is issued if it is not found.
#' @param N.sim Integer value. If specified and \code{object} is based on a 
#'   frequentist model (i.e., does not have a posterior sample), then a fake 
#'   posterior sample is generated using the function \code{sim}.
#' @param sim A function of three arguments (no names are assumed).
#'   If \code{N.sim} is supplied with a frequentist model, this function is called
#'   with respective arguments \code{N.sim}, \code{object@bhat}, and \code{object@V}.
#'   The default is the multivariate normal distribution.
#' @param ... Ignored.
#' 
#' @section Degrees of freedom:  
#' In cases where the
#' degrees of freedom depended on the linear function being estimated (e.g.,
#' Satterthwaite method), the d.f.
#' from the reference grid are saved, and a kind of \dQuote{containment} method
#' is substituted in the returned object, whereby the calculated d.f. for a new
#' linear function will be the minimum d.f. among those having nonzero
#' coefficients. This is kind of an \emph{ad hoc} method, and it can
#' over-estimate the degrees of freedom in some cases. An annotation is
#' displayed below any subsequent summary results stating that the 
#' degrees-of-freedom method is inherited from the previous method at
#' the time of re-gridding.
#'
#' @note Another way to use \code{regrid} is to supply a \code{regrid} 
#'   argument to \code{\link{ref_grid}} (either directly of indirectly via
#'   \code{\link{emmeans}}), in which case its value is passed to \code{regrid} as
#'   \code{transform}. This is often a simpler approach if the reference
#'   grid has not already been constructed.
#'
#' @return An \code{emmGrid} object with the requested changes
#' @export
#'
#' @examples
#' pigs.lm <- lm(log(conc) ~ source + factor(percent), data = pigs)
#' rg <- ref_grid(pigs.lm)
#' 
#' # This will yield EMMs as GEOMETRIC means of concentrations:
#' (emm1 <- emmeans(rg, "source", type = "response"))
#' pairs(emm1) ## We obtain RATIOS
#' 
#' # This will yield EMMs as ARITHMETIC means of concentrations:
#' (emm2 <- emmeans(regrid(rg, transform = "response"), "source"))
#' pairs(emm2)  ## We obtain DIFFERENCES
#' # Same result, useful if we hadn't already created 'rg'
#' # emm2 <- emmeans(pigs.lm, "source", regrid = "response")
#' 
#' # Simulate a sample of regression coefficients
#' set.seed(2.71828)
#' rgb <- regrid(rg, N.sim = 200, transform = "pass")
#' emmeans(rgb, "source", type = "response")  ## similar to emm1
regrid = function(object, transform = c("response", "mu", "unlink", "none", "pass", links), 
                  inv.link.lbl = "response", predict.type, 
                  bias.adjust = get_emm_option("back.bias.adj"), sigma, 
                  N.sim, sim = mvtnorm::rmvnorm, ...) 
{
    object = .chk.list(object, ...)
    links = c("logit", "probit", "cauchit", "cloglog", "log", "log10", "log2", "sqrt", "1/mu^2", "inverse")
    if (is.logical(transform))   # for backward-compatibility
        transform = ifelse(transform, "response", "none")
    else if (is.list(transform)) {
        userlink = transform
        transform = "user"
    }
    else
        transform = match.arg(transform)
    
    if (is.na(object@post.beta[1]) && !missing(N.sim)) {
        message("Simulating a sample of size ", N.sim, " of regression coefficients.")
        bh = object@bhat[!is.na(object@bhat)]
        object@post.beta = sim(N.sim, bh, object@V)
    }

    if (transform == "pass")
        return(object)

    # if we have two transformations to undo, do the first one recursively
    if ((transform == "response") && (!is.null(object@misc$tran2)))
             object = regrid(object, transform = "mu")

    # Save post.beta stuff
    PB = object@post.beta
    NC = attr(PB, "n.chains")
    
    if (!is.na(PB[1])) { # fix up post.beta BEFORE we overwrite parameters
        PB = PB %*% t(object@linfct[ , !is.na(object@bhat), drop = FALSE])
        if (".offset." %in% names(object@grid))
            PB = t(apply(PB, 1, function(.) . + object@grid[[".offset."]]))
    }
    
    est = .est.se.df(object, do.se = TRUE) ###FALSE)
    
    estble = !(is.na(est[[1]]))
    object@V = vcov(object)[estble, estble, drop = FALSE]
    object@bhat = est[[1]]
    object@linfct = diag(1, length(estble))
    object@misc$regrid.flag = TRUE
    pargs = object@grid[names(object@levels)]
    lbls = do.call(paste, c(pargs, sep = "."))
    if (!is.null(disp <- object@misc$display)) {  # fix up for the bookkeeping in nested models
        object@V = object@V[estble, estble, drop = FALSE]
        object@linfct = matrix(0, nrow = length(disp), ncol = length(estble))
        object@linfct[disp, ] = diag(1, length(estble))
        lbls = lbls[disp]
    }
    colnames(object@linfct) = lbls
    if(all(estble))
        object@nbasis = estimability::all.estble
    else
        object@nbasis = object@linfct[, !estble, drop = FALSE]
    
    # override the df function
    df = est$df
    edf = df[estble]
    if (length(edf) == 0) edf = NA
    # note both NA/NA and Inf/Inf test is.na() = TRUE
    prev.df.msg = attr(object@dffun, "mesg")
    if (any(is.na(edf/edf)) || (diff(range(edf)) < .01)) { # use common value
        object@dfargs = list(df = mean(edf, na.rm = TRUE))
        object@dffun = function(k, dfargs) dfargs$df
    }
    else { # use containment df
        object@dfargs = list(df = df)
        object@dffun = function(k, dfargs) {
            idx = which(zapsmall(k) != 0)
            ifelse(length(idx) == 0, NA, min(dfargs$df[idx], na.rm = TRUE))
        }
    }
    if(!is.null(prev.df.msg)) 
        attr(object@dffun, "mesg") = ifelse(
            startsWith(prev.df.msg, "inherited"), prev.df.msg,
            paste("inherited from", prev.df.msg, "when re-gridding"))
    
    if(transform %in% c("response", "mu", "unlink", links, "user") && !is.null(object@misc$tran)) {
        flink = link = attr(est, "link")
        if (bias.adjust) {
            if(missing(sigma))
                sigma = object@misc$sigma
            link = .make.bias.adj.link(link, sigma)
            #### (not needed)   bias.adjust = attr(link, "bias.adjust")
            if (!is.na(PB[1])) # special frequentist version when sigma is MCMC sample
                flink = .make.bias.adj.link(flink, mean(sigma))
            else
                flink = link
        }
        if (!is.na(PB[1]))
            PB = matrix(flink$linkinv(PB), ncol = ncol(PB))
        D = flink$mu.eta(object@bhat[estble])
        object@bhat = flink$linkinv(object@bhat)
        # efficient repl for D'VD with D <- diag(D)
        object@V = sweep(sweep(object@V, 1, D, "*"), 2, D, "*")

        inm = object@misc$inv.lbl
        if (!is.null(inm)) {
            object@misc$estName = inm
            if (!is.null(object@misc$log.contrast) && object@misc$log.contrast) # relabel ratios
                for (v in setdiff(object@misc$pri.vars, object@misc$by.vars))
                    object@grid[[v]] = gsub(" - ", "/", object@grid[[v]])
        }
        if((transform %in% c("mu", "unlink")) && !is.null(object@misc$tran2)) {
            object@misc$tran = object@misc$tran2
            object@misc$tran2 = object@misc$tran.mult = object@misc$tran.offset = object@misc$inv.lbl = NULL
        }
        else
            object@misc$tran = object@misc$tran.mult = object@misc$tran.offset = object@misc$inv.lbl = NULL
        sigma = object@misc$sigma = NULL
    }
    if (transform %in% c(links, "user")) { # fake a transformation
        link = if (transform == "user") userlink
               else                     .make.link(transform)
        bounds = range(link$linkinv(c(-1e6, -100, -1, 0, 1, 100, 1e6)))
        nas = which(is.na(object@bhat)) # already NA
        incl = vincl = which((object@bhat > bounds[1]) & (object@bhat < bounds[2]))
        if (length(nas) > 0) 
            vincl = which((object@bhat[-nas] > bounds[1]) & (object@bhat[-nas] < bounds[2]))
        negs = setdiff(seq_along(object@bhat), incl)
        if (length(negs) > 0) {
            message("Invalid response predictions are flagged as non-estimable")
            object@bhat[negs] = NA
            tmp = seq_along(object@bhat)
            object@nbasis = sapply(c(nas, negs), function(ii) 0 + (tmp == ii))
        }
        object@bhat = link$linkfun(object@bhat)
        Vee = object@V
        if(length(incl) > 0) {
            D = 1 / link$mu.eta(object@bhat[incl])
            object@V = sweep(sweep(Vee[vincl, vincl, drop = FALSE], 1, D, "*"), 2, D, "*")
        }
        if (!is.na(PB[1])) {
            PB[PB <= 0] = NA
            PB = link$linkfun(PB)
            PB[1] = ifelse(is.na(PB[1]), 0, PB[1]) # make sure 1st elt isn't NA
        }
        object@misc$tran = if (transform == "user") link   
                           else                     transform
        object@misc$inv.lbl = inv.link.lbl
    }
    
    if(!is.na(PB[1])) {
        PB = PB[ , !is.na(object@bhat), drop = FALSE]
        attr(PB, "n.chains") = NC
        object@post.beta = PB
    }
    
    # Nix out things that are no longer needed or valid
    object@grid$.offset. = object@misc$offset.mult =
        object@misc$estHook = object@misc$vcovHook = NULL
    object@model.info$model.matrix = "Submodels are not available with regridded objects"
    if(!missing(predict.type))
        object = update(object, predict.type = predict.type)

    object
}

