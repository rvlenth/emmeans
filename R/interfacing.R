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

#' Support functions for model extensions
#' 
#' This documents the methods that \code{\link{ref_grid}} calls. A user
#' or package developer may add \pkg{emmeans} support for a model
#' class by writing \code{recover_data} and \code{emm_basis} methods
#' for that class.
#' 
## #' @rdname extending-emmeans
#' @name extending-emmeans
#' @param object An object of the same class as is supported by a new method.
#' @param ... Additional parameters that may be supported by the method.
#' 
#' @section Details:
#' To create a reference grid, the \code{ref_grid} function needs to reconstruct
#' the data used in fitting the model, and then obtain a matrix of linear
#' functions of the regression coefficients for a given grid of predictor
#' values. These tasks are performed by calls to \code{recover_data} and
#' \code{emm_basis} respectively. A vignette giving details and examples
#' is available via \href{../doc/xtending.html}{vignette("xtending", "emmeans")}
#' 
#' To extend \pkg{emmeans}'s support to additional model types, one need only
#' write S3 methods for these two functions. The existing methods serve as
#' helpful guidance for writing new ones.  Most of the work for
#' \code{recover_data} can be done by its method for class \code{"call"},
#' providing the \code{terms} component and \code{na.action} data as additional
#' arguments. Writing an \code{emm_basis} method is more involved, but the
#' existing methods (e.g., \code{emmeans:::emm_basis.lm}) can serve as models.
#' Certain \code{recover_data} and \code{emm_basis} methods are exported from
#' \pkg{emmeans}. (To find out, do \code{methods("recover_data")}.) If your
#' object is based on another model-fitting object, it
#' may be that all that is needed is to call one of these exported methods and
#' perhaps make modifications to the results. Contact the developer if you need
#' others of these exported.
#' 
#' If the model has a multivariate response, \code{bhat} needs to be
#' \dQuote{flattened} into a single vector, and \code{X} and \code{V} must be
#' constructed consistently.
#' 
#' In models where a non-full-rank result is possible (often, you can tell by
#' seeing if there is a \code{singular.ok} argument in the model-fitting
#' function), \code{\link{summary.emmGrid}} and its relatives check the
#' estimability of each
#' prediction, using the \code{\link[estimability]{nonest.basis}} function in
#' the \pkg{estimability} package.
#' 
#' The models already supported are detailed in \href{../doc/models.html}{the
#' "models" vignette}. Some packages may provide additional \pkg{emmeans}
#' support for its object classes.
#'
#'
#' @return The \code{recover_data} method must return a \code{\link{data.frame}}
#'   containing all the variables that appear as predictors in the model,
#'   and attributes \code{"call"}, \code{"terms"}, \code{"predictors"},
#'   and \code{"responses"}. (\code{recover_data.call} will 
#'   provide these attributes.)
#' 
#' @note Without an explicit \code{data} argument, \code{recover_data} returns
#'    the \emph{current version} of the dataset. If the dataset has changed
#'    since the model was fitted, then this will not be the data used to fit
#'    the model. It is especially important to know this in simulation studies
#'    where the data are randomly generated or permuted, and in cases where
#'    several datasets are processed in one step (e.g., using \code{dplyr}).
#'    In those cases, users should be careful to provide the actual data
#'    used to fit the model in the \code{data} argument.
#'   
#' @seealso \href{../doc/xtending.html}{Vignette on extending emmeans}
#' 
#' @export
recover_data = function(object, ...) {
    # look for outside methods first
    for (cl in .chk.cls(object)) {
        rd <- try(getS3method("recover_data", cl, envir = .GlobalEnv), silent = TRUE)
        if(!inherits(rd, "try-error"))
            return(rd(object, ...))
    }
    UseMethod("recover_data")
}

# get classes that are OK for external code to modify
# We don't allow overriding certain anchor classes, 
# nor ones in 3rd place or later in inheritance
.chk.cls = function(object) {
    sacred = c("call", "lm", "glm", "mlm", "aovlist", "lme", "qdrg")
    setdiff(class(object)[1:2], sacred)
}







#--------------------------------------------------------------
### call' objects
# This recover_data method serves as the workhorse for the others
# For model objects, call this with the object's call and its terms component
# Late addition: if data is non-null, use it in place of recovered data
# Later addition: na.action arg req'd - vector of row indexes removed due to NAs
#    na.action is ignored when data is non-NULL

#' @rdname extending-emmeans
#' @param trms The \code{\link{terms}} component of \code{object} (typically with
#'   the response deleted, e.g. via \code{\link{delete.response}})
#' @param na.action Integer vector of indices of observations to ignore; or
#'   \code{NULL} if none
#' @param data Data frame. Usually, this is \code{NULL}. However, if non-null,
#'   this is used in place of the reconstructed dataset. It must have all of the
#'   predictors used in the model, and any factor levels must match those used
#'   in fitting the model.
#' @param params Character vector giving the names of any variables in the model
#'   formula that are \emph{not} predictors. For example, a spline model may involve
#'   a local variable \code{knots} that is not a predictor, but its value is
#'   needed to fit the model. Names of parameters not actually used are harmless,
#'   and the default value \code{"pi"} (the only numeric constant in base R)
#'   is provided in case the model involves it.
#' 
#' @method recover_data call
#' @export
recover_data.call = function(object, trms, na.action, data = NULL, params = "pi", ...) {
    fcall = object # because I'm easily confused
    vars = setdiff(.all.vars(trms), params)
    tbl = data
    if (length(vars) == 0 || vars[1] == "1") {
        tbl = data.frame(c(1,1))
        vars = names(tbl) = 1
    }
    if (is.null(tbl)) {
        possibly.random = FALSE
        m = match(c("formula", "data", "subset", "weights"), names(fcall), 0L)
        fcall = fcall[c(1L, m)]
        
        # check to see if there are any function calls to worry about
        # [e.g., subset = sample(1:n, 50) will give us a 
        #    different subset than model used]
        mm = match(c("data", "subset"), names(fcall), 0L)
        if(any(mm > 0)) {
            # Flag cases where there is a function call in data or subset
            # May indicate a situation where data are randomized
            fcns = unlist(lapply(fcall[mm], 
                     function(x) setdiff(all.names(x), 
                                         c("::",":::","[[","]]",all.vars(x)))))
            possibly.random = (max(nchar(c("", fcns))) > 1)
        }
        
        fcall$drop.unused.levels = TRUE
        fcall[[1L]] = as.name("model.frame")
        fcall$xlev = NULL # we'll ignore xlev
        
        if(!is.numeric(na.action))   ### In case na.action is not a vector of indices
            na.action = NULL
        
        # If we have an explicit list of cases to exclude, let everything through now
        if (!is.null(na.action))
            fcall$na.action = na.pass
        else  # exclude incomplete cases
            fcall$na.action = na.omit
        form = .reformulate(vars)
        fcall$formula = update(trms, form)
        env = environment(trms)
        if (is.null(env)) 
            env = parent.frame()
        tbl = try(eval(fcall, env, parent.frame()), silent = TRUE)
        if(inherits(tbl, "try-error"))
            return(.rd.error(vars, fcall))
        if (possibly.random) {
            chk = eval(fcall, env, parent.frame())
            if (!all(chk == tbl))
                stop("Data appear to be randomized -- ", 
                     "cannot consistently recover the data\n",
                     "Move the randomization ",
                     "outside of the model-fitting call.")
        }
        
        # Now we can drop na.action's rows
        if (!is.null(na.action))
            tbl = tbl[-(na.action),  , drop=FALSE]
    }
    
    else {
        tbl = tbl[, vars, drop = FALSE] # consider only the variables actually needed
        tbl = tbl[complete.cases(tbl), , drop=FALSE]
    }
    
    attr(tbl, "call") = object # the original call
    attr(tbl, "terms") = trms
    attr(tbl, "predictors") = setdiff(.all.vars(delete.response(trms)), params)
    attr(tbl, "responses") = setdiff(vars, union(attr(tbl, "predictors"), params))
    tbl
}

# error message for recover_data.call
.rd.error = function(vars, fcall) {
    if ("pi" %in% vars) 
        return("\nTry re-running with 'params = c\"pi\", ...)'")
    if (is.list(fcall$data)) fcall$data = "(raw data structure)"
    dataname = as.character(fcall$data)[1]
    if ((!is.na(dataname)) && (nchar(dataname) > 50))
        dataname = paste(substring(dataname, 1, 50), "...")
    mesg = "We are unable to reconstruct the data.\n"
    mesg = paste0(mesg, "The variables needed are:\n\t",
           paste(vars, collapse = " "), "\n",
           "Are any of these actually constants? (specify via 'params = ')\n")
    if (is.na(dataname))
        mesg = paste(mesg, "Try re-running with 'data = \"<name of dataset>\"'\n")
    else 
        mesg = paste0(mesg, "The dataset name is:\n\t", dataname, "\n",
           "Does the data still exist? Or you can specify a dataset via 'data = '\n")
    mesg
}

#----------------------------------------------------------
### emm_basis methods create a basis for the reference grid
#
# Required args:
#     object - the model object
#     trms   - terms component of object
#     xlev   - named list of factor levels (but not the coerced ones)
#     grid   - reference grid
# All methods must return a list with these elements:
#     X      - basis for linear fcns over reference grid
#     bhat   - regression coefficients for fixed effects (INCLUDING any NAs)
#     nbasis - matrix whose columns for a basis for non-estimable functions of beta; all.estble if none
#     V      - estimated covariance matrix of bhat
#     dffun  - function(k, dfargs) to find df for k'bhat having std error se
#     dfargs - additional arguments, if any, for dffun
#     misc   - Extra info ...
#              -- if extra levels need to be added (e.g. mlm, polr),
#                 put them in misc$ylevs
#              -- For transformations or link fcns, use misc$tran
#                 for name (see 'make.link'), and use misc$inv.lbl
#                 for label to use in 'summary' when tran is inverted
#                 (ref_grid looks at lhs of model for tran if none found)
# Note: if no df exists, set dffun = function(...) NA and dfargs = list()
#--------------------------------------------------------------

#' @rdname extending-emmeans
#' @param xlev Named list of factor levels (\emph{excluding} ones coerced to 
#'   factors in the model formula)
#' @param grid A \code{data.frame} (provided by \code{ref_grid}) containing 
#'   the predictor settings needed in the reference grid
#'
#' @return The \code{emm_basis} method should return a \code{list} with the
#'   following elements:
#' \describe{
#' \item{X}{The matrix of linear functions over \code{grid}, having the same
#'   number of rows as \code{grid} and the number of columns equal to the length
#'   of \code{bhat}.}
#' \item{bhat}{The vector of regression coefficients for fixed effects. This
#'   should \emph{include} any \code{NA}s that result from rank deficiencies.}
#' \item{nbasis}{A matrix whose columns form a basis for non-estimable functions
#'   of beta, or a 1x1 matrix of \code{NA} if there is no rank deficiency.}
#' \item{V}{The estimated covariance matrix of \code{bhat}.}
#' \item{dffun}{A function of \code{(k, dfargs)} that returns the degrees of
#'   freedom associated with \code{sum(k * bhat)}.}
#' \item{dfargs}{A \code{list} containing additional arguments needed for
#'   \code{dffun}}.
#' } %%% end of describe
#' @export
#' 
#' @section Communication between methods:
#' If the \code{recover_data} method generates information needed by \code{emm_basis},
#' that information may be incorporated by creating a \code{"misc"} attribute in the
#' returned recovered data. That information is then passed as the \code{misc} 
#' argument when \code{ref_grid} calls \code{emm_basis}.
#' 
#' @section Optional hooks:
#' Some models may need something other than standard linear estimates and
#' standard errors. If so, custom functions may be pointed to via the items
#' \code{misc$estHook}, \code{misc$vcovHook} and \code{misc$postGridHook}. If
#' just the name of the hook function is provided as a character string, then it
#' is retrieved using \code{\link{get}}.
#' 
#' The \code{estHook} function should have arguments \samp{(object, do.se, tol,
#' ...)} where \code{object} is the \code{emmGrid} object,
#' \code{do.se} is a logical flag for whether to return the standard error, and
#' \code{tol} is the tolerance for assessing estimability. It should return a
#' matrix with 3 columns: the estimates, standard errors (\code{NA} when
#' \code{do.se==FALSE}), and degrees of freedom (\code{NA} for asymptotic). The
#' number of rows should be the same as \samp{object@linfct}. The
#' \code{vcovHook} function should have arguments \samp{(object, tol, ...)} as
#' described. It should return the covariance matrix for the estimates. Finally,
#' \code{postGridHook}, if present, is called at the very end of
#' \code{ref_grid}; it takes one argument, the constructed \code{object}, and
#' should return a suitably modified \code{emmGrid} object.
emm_basis = function(object, trms, xlev, grid, ...) {
    # look for outside methods first
    for (cl in .chk.cls(object)) {
        emb <- try(getS3method("emm_basis", cl, envir = .GlobalEnv), silent = TRUE)
        if(!inherits(emb, "try-error"))
            return(emb(object, trms, xlev, grid, ...))
    }
    UseMethod("emm_basis")
}

# Hidden courtesy function that provides access to all recover_data methods
#' @rdname extending-emmeans
#' @export
.recover_data = function(object, ...)
    recover_data(object, ...)

# Hidden courtesy function that provides access to all emm_basis methods
#' @rdname extending-emmeans
#' @return \code{.recover_data} and \code{.emm_basis} are hidden exported versions of 
#'   \code{recover_data} and \code{emm_basis}, respectively. They run in \pkg{emmeans}'s
#'   namespace, thus providing access to all existing methods.
#' @export
.emm_basis = function(object, trms, xlev, grid, ...)
    emm_basis(object, trms, xlev, grid, ...)




#--------------------------------------------------------------
### DEFAULT METHODS (we hit these when a model is NOT supported)
# I'll have it return the message if we caught the error in this way
# Then caller can use try() to check for other types of errors,
# and just print this message otherwise 
# NOT @exported
recover_data.default = function(object, ...) {
    paste("Can't handle an object of class ", dQuote(class(object)[1]), "\n",
          paste(.show_supported(), collapse=""))
}
# NOT @exported
emm_basis.default = function(object, trms, xlev, grid, ...) {
    stop("Can't handle an object of class", dQuote(class(object)[1]), "\n",
         .show_supported())
}

# Private fcn to get a list of supported objects
# does this by looking in namespace [ns] and methods [meth]
# then strips that off leaving extensions
.show_supported = function(ns = "emmeans", meth = "emm_basis") {
    "Use help(\"models\", package = \"emmeans\") for information on supported models."
}


