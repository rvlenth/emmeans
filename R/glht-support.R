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

### Code for an enhancement of 'glht' in 'multcomp' package
### Provides for using 'emm' in similar way to 'mcp'
### This is implemented via the class "emmlf" -- linear functions for emmeans

## NOTE: Registration of S3 methods for glht is done dynamically in zzz.R

# emm(specs) will be used as 'linfct' argument in glht
# all we need to do is class it and save the arguments

#' Support for \code{multcomp::glht}
#' 
#' These functions and methods provide an interface between \pkg{emmeans} and
#' the \code{multcomp::glht} function for simultaneous inference provided
#' by the \pkg{multcomp} package.
#' 
#' \code{emm} is meant to be called only \emph{from} \code{"glht"} as its second
#' (\code{linfct}) argument. It works similarly to \code{multcomp::mcp},
#' except with \code{specs} (and optionally \code{by} and \code{contr}
#' arguments) provided as in a call to \code{\link{emmeans}}.
#' 
#' @rdname glht-support
#' @aliases glht-support glht.emmGrid glht.emmlf modelparm.emmwrap
#' @param ... In \code{emm}, the \code{specs}, \code{by}, and \code{contr}
#'   arguments you would normally supply to \code{\link{emmeans}}. Only
#'   \code{specs} is required. Otherwise, arguments that are passed to other
#'   methods.
#'
#' @return \code{emm} returns an object of an intermediate class for which
#'   there is a \code{multcomp::glht} method.
#' @export
emm <- function(...) {
    result <- list(...)
    class(result) <- "emmlf"
    result
}

# New S3 method for emmlf objects
glht.emmlf <- function(model, linfct, ...) {
    # Pass the arguments we should pass to ref_grid:
    args = linfct
    args[[1]] = model
    names(args)[1] = "object"
    # Now pass the ref_grid to emmeans:
    linfct$object <- do.call("ref_grid", args)
    emmo <- do.call("emmeans", linfct)
    if (is.list(emmo)) 
        emmo = emmo[[length(emmo)]]
    # Then call the method for emmo objject
    glht.emmGrid(model, emmo, ...)
}


# S3 method for an emmGrid object
# Note: model is redundant, really, so can be omitted
# See related roxygen stuff just before glht.emmlf
glht.emmGrid <- function(model, linfct, by, ...) {
    .requireNS("multcomp", sQuote("glht")," requires ", dQuote("multcomp"), 
               " to be installed", call. = FALSE)
    object = linfct # so I don't get confused
    if (missing(model)) 
        model = .cls.list("emmwrap", object = object)
    args = list(model = model, ...)
    # add a df value if not supplied
    if (is.null(args$df)) {
        df = summary(linfct)$df
        df[is.infinite(df)] = NA
        if(any(!is.na(df))) {
            args$df = max(1, as.integer(mean(df, na.rm=TRUE) + .25))
            if (any(args$df != df))
                message("Note: df set to ", args$df)
        }
    }
    if (missing(by)) by = object@misc$by.vars
    
    nms = setdiff(names(object@grid), c(by, ".offset.", ".freq.", ".wgt."))
    if (is.null(object@misc$estHook))
        lf = object@linfct
    else # custom estimation setup - use the grid itself as the parameterization
        lf = diag(1, nrow(object@linfct))
    dimnames(lf)[[1]] = as.character(interaction(object@grid[, nms], sep=", "))
    
    if (is.null(by)) {
        args$linfct = lf
        return(do.call(multcomp::glht, args))
    }
    
    # (else...)
    by.rows = .find.by.rows(object@grid, by)
    result = lapply(by.rows, function(r) {
        args$linfct = lf[r, , drop=FALSE]
        do.call(multcomp::glht, args)
    })
    bylevs = lapply(by, function(byv) unique(object@grid[[byv]]))
    names(bylevs) = by
    bygrid = do.call("expand.grid", bylevs)
    levlbls = unname(lapply(by, function(byv) paste(byv, "=", bygrid[[byv]])))
    levlbls$sep = ", "
    names(result) = do.call("paste", levlbls)
    class(result) = c("glht_list", "list")
    result
}

### as. glht -- convert my object to glht object
#' @rdname glht-support
#' @param object An object of class \code{emmGrid} or \code{emm_list}
#' 
#' @return \code{as.glht} returns an object of class \code{glht} or \code{glht_list}
#'   according to whether \code{object} is of class \code{emmGrid} or \code{emm_list}. 
#'   See Details below for more on \code{glht_list}s.
#'   
#' @section Details:
#' A \code{glht_list} object is simply a \code{list} of \code{glht} objects. 
#' It is created as needed -- for example, when there is a \code{by} variable. 
#' Appropriate convenience methods \code{coef},
#' \code{confint}, \code{plot}, \code{summary}, and \code{vcov} are provided,
#' which simply apply the corresponding \code{glht} methods to each member.
#' 
#' @note The multivariate-\eqn{t} routines used by \code{glht} require that all
#'   estimates in the family have the same integer degrees of freedom. In cases
#'   where that is not true, a message is displayed that shows what df is used.
#'   The user may override this via the \code{df} argument.
#' 
#' @examples
#' if(require(multcomp)) { # --- multcomp must be installed
#' 
#' warp.lm <- lm(breaks ~ wool*tension, data = warpbreaks)
#' 
#' # Using 'emm'
#' summary(glht(warp.lm, emm(pairwise ~ tension | wool)))
#' 
#' # Same, but using an existing 'emmeans' result
#' warp.emm <- emmeans(warp.lm, ~ tension | wool)
#' summary(as.glht(pairs(warp.emm)))
#' 
#' # Same contrasts, but treat as one family
#' summary(as.glht(pairs(warp.emm), by = NULL))
#' 
#' } # --- was tested only if multcomp is installed
#' @export
as.glht <- function(object, ...) {
    UseMethod("as.glht")
}

#' @method as.glht default
#' @export
as.glht.default <- function(object, ...)
    stop("Cannot convert an object of class ", sQuote(class(object)[1]),
         " to a ", sQuote("glht"), " object")

#' @rdname glht-support
#' @method as.glht emmGrid
#' @export
as.glht.emmGrid <- function(object, ...)
    glht.emmGrid( , object, ...)   # 1st arg not necessary

#' @method as.glht emm_list
#' @export
as.glht.emm_list <- function(object, ..., which = 1)
    as.glht(object[[which]], ...)


# S3 modelparm method for emmwrap (S3 wrapper for an emmGrid obj - see glht.emmGrid)
#--- dynamically registered in zzz.R --- #' @export
modelparm.emmwrap <- function(model, coef., vcov., df, ...) {
    object = model$object
    if (is.null(object@misc$estHook)) {
        bhat = object@bhat
        V = object@V
    }
    else { # Have custom vcov and est methods. Use the grid itself as parameterization
        bhat = predict(object)
        V = vcov(object)
    }
    if(missing(df) || is.na(df) || is.infinite(df))
        df = 0
    .cls.list("modelparm", coef = bhat, vcov = V,
                df = df, estimable = !is.na(bhat))
    # This is NOT what we mean by 'estimable', but it is what glht wants...
}


# S3 methods for glht_list

### Doesn't work so excluded...
# cld.glht_list = function(object, ...)
#     lapply(object, cld, ...)

#' @method coef glht_list
#' @export
coef.glht_list = function(object, ...)
    lapply(object, coef, ...)

#' @method confint glht_list
#' @export
confint.glht_list = function(object, ...)
    lapply(object, confint, ...)

#' @method plot glht_list
#' @export
plot.glht_list = function(x, ...)
    lapply(x, plot, ...)

#' @method summary glht_list
#' @export
summary.glht_list = function(object, ...)
    lapply(object, summary, ...)

#' @method vcov glht_list
#' @export
vcov.glht_list = function(object, ...)
    lapply(object, vcov, ...)






