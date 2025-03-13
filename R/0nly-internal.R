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

# Functions used only internally and of fairly general use
# More of these with specific uses only within a certain context remain there.

## create an lm-compatible set of coefficients in singular cases
## requires both bhat and X to be named; if not, silently returns bhat
.impute.NAs = function(bhat, X) {
    if (is.null(bnm <- names(bhat)) || is.null(colnames(X)))
        return(bhat)
    nm = intersect(bnm, (xnm <- colnames(X)))
    if(length(nm) == length(xnm))
        return(bhat[nm])
    tmp = rep(NA, length(xnm))
    names(tmp) = xnm
    tmp[nm] = bhat[nm]
    tmp
}

## %in%-style operator with partial matching
## e.g.,  ("bonf" %.pin% p.adjust.methods)  is TRUE
"%.pin%" = function (x, table) pmatch(x, table, nomatch = 0L) > 0L

## Like is.numeric() but returns TRUE for character vectors like c("1", "3", "2")
## (but will return FALSE if any x is NA)
.is.num = function(x) {
    nx = suppressWarnings(as.numeric(as.character(x)))
    !any(is.na(nx))
}

# return TRUE if obj@misc$sigma is not NA
.chk.predict = function(obj) {
    sig = obj@misc$sigma
    ok = (is.null(sig) || (!is.null(sig) && !is.na(sig)))
    if(!ok)
        warning("Prediction intervals are not available for this object", call. = FALSE)
    ok
}

### Internal function to implement 'allow.na.levs' option
.chk.fac = function(x) {
    if(get_emm_option("allow.na.levs"))
        factor(x, exclude = NULL)
    else
        factor(x)
}



## Alternative to all.vars, but keeps vars like foo$x and foo[[1]] as-is
##   Passes ... to all.vars
#' @rdname extending-emmeans
#' @order 41
#' @param expr,retain Arguments for \code{.all.vars}, which is an alternative to \code{\link{all.vars}}
#'   that has special provisions for retaining the special characters in \code{retain},
#'   thus allowing model specifications like \code{y ~ data$trt * df[["dose"]]}
#' @export
.all.vars = function(expr, retain = c("\\$", "\\[\\[", "\\]\\]", "'", '"'), ...) {
    if (is.null(expr) || length(expr) == 0)
        return(character(0))
    if (!inherits(expr, "formula")) {
        expr = try(eval(expr), silent = TRUE)
        if(inherits(expr, "try-error")) {
            return(character(0))
        }
    }
    repl = paste(".Av", seq_along(retain), ".", sep = "")
    for (i in seq_along(retain))
        expr = gsub(retain[i], repl[i], expr)
    subs = switch(length(expr), 1, c(1,2), c(2,1,3))
    vars = all.vars(as.formula(paste(expr[subs], collapse = "")), ...)
    retain = gsub("\\\\", "", retain)
    for (i in seq_along(retain))
        vars = gsub(repl[i], retain[i], vars)
    if(length(vars) == 0) vars = "1"   # no vars ---> intercept
    vars
}

### returns TRUE iff there is one or more function call in a formula
.has.fcns = function(form) {
    fcns = setdiff(.all.vars(form, functions = TRUE), 
                   c("~", "+", "-", "*", "/", ":", "(", "|", .all.vars(form)))
    length(fcns) > 0
}

### parse a formula of the form lhs ~ rhs | by into a list
### of variable names in each part
### Returns character(0) for any missing pieces
.parse.by.formula = function(form) {
    allv = .all.vars(form)
    ridx = ifelse(length(form) == 2, 2, 3)
    allrhs = as.vector(form, "character")[ridx]
    allrhs = gsub("\\|", "+ .by. +", allrhs) # '|' --> '.by.'
    allrhs = .all.vars(stats::reformulate(allrhs))
    bidx = grep(".by.", allrhs, fixed = TRUE)
    if (length(bidx) == 0) { # no '|' in formula
        by = character(0)
        rhs = allrhs
    }
    else {
        rhs = allrhs[seq_len(bidx[1] - 1)]
        by = setdiff(allrhs, c(rhs, ".by."))
    }
    lhs = setdiff(allv, allrhs)
    list(lhs = lhs, rhs = rhs, by = by)
}

### Check specs for 'all' and if there, create a list of specs
# Return the specs as they were, or if 'all' is specified, return revised list of specs
# with 'form.rtn' attribute
.parse.specs.for.all = function(object, specs, by) {
    all.key = "."  # The key to use to request all sets of means
    if(is.list(specs))
        return(specs)
    rtn = NULL
    if (inherits(specs, "formula")) {
        rtn = .parse.by.formula(specs)
        if ((length(rtn$rhs) == 0) || (rtn$rhs[1] != all.key))
            return(specs)
        specs = rtn$rhs
        rtn$by = setdiff(rtn$by, specs)
        if(length(rtn$by) > 0)
            by = rtn$by
    }
    # now we have is.character(specs) ...
    if ((length(specs) != 1) || (specs[1] != all.key))
        return(specs)
    # hack the object to bypass estimability checking. 
    # Doesn't matter what stats are as we use only the labels
    if(any(is.na(object@bhat))) {
        k = ncol(object@linfct)
        object@bhat = seq_len(k)
        object@V = diag(k)
        object@nbasis = estimability::all.estble
    }
    stgs = unique(setdiff(joint_tests(object, by = by)[, "model term"], "(confounded)"))
    # I think we need to add nested factors that involve 'by'
    if ((length(by) > 0) && !is.null(nst <- object@model.info$nesting)) {
        for (fac in names(nst))
            if (all(by %in% nst[[fac]]))
                stgs = c(stgs, fac)
        stgs = unique(stgs)
    }
    if(length(stgs) == 0)
        stop("'", all.key, "' specification yielded no terms", call. = FALSE)
    
    result = strsplit(stgs, ":")
    if(!is.null(rtn))
        attr(result, "form.rtn") = rtn
    result
}


# Utility to pick out the args that can be passed to a function
.args.for.fcn = function(fcn, args) {
    oknames = names(as.list(args(fcn)))
    mat = pmatch(names(args), oknames)
    args = args[!is.na(mat)]
    mat = mat[!is.na(mat)]
    names(args) = oknames[mat]
    args
}

# Create a list and give it class class.name
.cls.list <- function(class.name, ...) {
    result <- list(...)
    class(result) <- c(class.name, "list")
    result
}

### Not-so-damn-smart replacement of diag() that will 
### not be so quick to assume I want an identity matrix
### returns matrix(x) when x is a scalar
#' @rdname extending-emmeans
#' @order 42
#' @param x,nrow,ncol Arguments for \code{.diag}, which is an alternative to 
#'   \code{\link{diag}} that lacks its idiosyncrasy of returning an
#'   identity matrix when \code{x} is of length 1.
#' @export
#' 
.diag = function(x, nrow, ncol) {
    if(is.matrix(x))
        diag(x)
    else if((length(x) == 1) && missing(nrow) && missing(ncol)) 
        matrix(x)
    else 
        diag(x, nrow, ncol)
}

# Utility that returns TRUE if getOption("emmeans")[[opt]] is TRUE
.emmGrid.is.true = function(opt) {
    x = get_emm_option(opt)
    if (is.logical(x))  x
    else FALSE
}


# return list of row indexes in tbl for each combination of by
# tbl should be a data.frame
.find.by.rows = function(tbl, by) {
    if (is.null(by))
        return(list(seq_len(nrow(tbl))))
    if (any(is.na(match(by, names(tbl)))))
        stop("'by' variables are not all in the grid")    
    bylevs = tbl[ , by, drop = FALSE]
    by.id = do.call("paste", unname(bylevs))
    uids = unique(by.id)
    result = lapply(uids, function(id) which(by.id == id))
    names(result) = uids
    result
}

# calculate the offset for the given grid
#' @rdname extending-emmeans
#' @order 35
#' @param terms A \code{terms} component
#' @return \code{.get.offset} returns the values, based on \code{grid}, of 
#' any \code{offset} component in \code{terms}
#' @export
.get.offset = function(terms, grid) {
    off.idx = attr(terms, "offset")
    offset = rep(0, nrow(grid))
    tvars = attr(terms, "variables")
    for (i in off.idx)
        offset = offset + eval(tvars[[i+1]], grid)
    offset
}

# combine variables in several `terms` objects
#' @rdname extending-emmeans
#' @order 31
#' @return \code{combine.terms} returns a \code{terms} object resulting
#'   from combining all the terms or formulas in \code{...}.
#' @export
.combine.terms = function(...) {
    trms = list(...)
    vars = unlist(lapply(trms, .all.vars))
    terms(.reformulate(vars, env = environment(trms[[1]])))
}

######################################################################
### Contributed by Jonathon Love, https://github.com/jonathon-love ###
### and adapted by RVL to exclude terms like df$trt or df[["trt"]] ###
######################################################################
# reformulate for us internally in emmeans
# same as stats::reformulate, except it surrounds term labels with backsticks
#
# RVL note: I renamed it .reformulate to avoid certain issues.
#   For example I need reformulate() sometimes to strip off function calls
#   and this .reformulate works quite differently.
#
.reformulate <- function (termlabels, response = NULL, intercept = TRUE, env = parent.frame())
{
    if (!is.character(termlabels) || !length(termlabels))
        stop("'termlabels' must be a character vector of length at least one")
    has.resp = !is.null(response)
    termlabels = sapply(trimws(termlabels), function(x) 
        if (length(grep("\\$|\\[\\[", x)) > 0) x
        else paste0("`", x, "`"))
    termtext = paste(if (has.resp) "response", "~", 
                     paste(termlabels, collapse = "+"), collapse = "")
# prev version:                     paste0("`", trimws(termlabels), "`", collapse = "+"), collapse = "")
    if (!intercept)
        termtext = paste(termtext, "- 1")
    rval = eval(parse(text = termtext, keep.source = FALSE)[[1L]])
    if (has.resp)
        rval[[2L]] = if (is.character(response))
            as.symbol(response)
    else response
    environment(rval) = env
    rval
}

### Find variable names of the form df$x or df[["x"]]
# Returns indexes. In addition, return value has attribute
#   "details"  matrix with 1st row being dataset names, second row is variable names
.find.comp.names = function(vars) {
    comp = grep("\\$|\\[\\[", vars) # untick vars containing "$"
    if (length(comp) > 0) {
        attr(comp, "details") = gsub("\"", "",
            sapply(strsplit(vars[comp], "\\$|\\[\\[|\\]\\]"), function(.) .[1:2]))
    }
    comp
}


# returns a list of all matches to ... or lst with full names from args
.match.dots.list = function(args, ..., lst) {
    if(missing(lst))
        lst = list(...)
    idx = pmatch(names(lst), args, nomatch = 0)
    rtn = lst[idx > 0]
    names(rtn) = args[idx]
    rtn
}

### Find single arg in `...`. If pmatched, return its value, else NULL
.match.dots = function(arg, ..., lst) {
    rtn = .match.dots.list(arg[1], ..., lst = lst)
    if (length(rtn) > 0) rtn[[1]] else NULL
}


# return a list from ..., omitting any that pmatch omit
.zap.args = function(..., omit) {
    args = list(...)
    args[!is.na(pmatch(names(args), omit))] = NULL
    args
}

# return updated object with option list AND dot list
# optionally may exclude any opts in 'exclude.opts'
.update.options = function(object, options, ..., exclude.opts) {
    if (!is.list(options))
        options = as.list(options)
    dot.opts = .match.dots.list(.valid.misc, ...)
    if (!missing(exclude.opts))
        for (nm in exclude.opts)
            dot.opts[[nm]] = NULL
    # entries in both lists are overridden by those in ...
    for (nm in names(dot.opts))
        options[[nm]] = dot.opts[[nm]]
    options[["object"]] = object
    do.call(update.emmGrid, options)
}

# my own model.frame function. Intercepts compound names
# and fixes up the data component accordingly. We do this
# by creating data.frames within data having required variables of simple names 
model.frame = function(formula, data, ...) {
    if (is.null(data))
        return (stats::model.frame(formula, ...))
    idx = .find.comp.names(names(data))
    if (length(idx) > 0) {
        nm = names(data)[idx]
        others = names(data[-idx])
        details =  attr(idx, "details")
        num = suppressWarnings(as.numeric(details[2, ]))
        data = as.list(data)
        for (dfnm in unique(details[1, ])) {
            w = which(details[1, ] == dfnm)
            data[[dfnm]] = as.data.frame(data[nm[w]])
            names(data[[dfnm]]) = details[2, w]
            data[[dfnm]] = .reorder.cols(data[[dfnm]], num[w])
        }
        data[nm] = NULL # toss out stuff we don't need
        # save subst table in environ
        if (get_emm_option("simplify.names")) {
            details[1, ] = nm # top row is now fancy names, bottom is plain names
            all.nms = c(details[2, ], others)
            dup.cnt = sapply(details[2, ], function(.) sum(. == all.nms))
            dup.cnt[!is.na(num)] = 3 # don't simplify numeric ones
            details = details[, dup.cnt == 1, drop = FALSE]
            if (ncol(details) > 0)
                environment(formula)$.simplify.names. = details
        }
    }
    stats::model.frame(formula, data = data, ...)
}

# Utility to simplify names. each elt of top row of tbl is changed to bottom row
.simplify.names = function(nms, tbl) {
    for (j in seq_along(tbl[1,]))
        nms[nms == tbl[1, j]] = tbl[2, j]
    nms
}

# utility to make all names in a summary syntactically valid
.validate.names = function(object) {
    for (a in c("names", "pri.vars", "by.vars"))
        if (!is.null(att <- attr(object, a))) 
            attr(object, a) = make.names(att)
    object
}

# reorder columns of data frame to match numeric positions in num (if any)
.reorder.cols = function(data, num) {
    if (all(is.na(num)))
        return(data)
    m = max(num[!is.na(num)])
    nm = names(data)
    if (m > length(nm)) {
        k = length(nm)
        data = cbind(data, matrix(NA, nrow = nrow(data), ncol = m - k))
        nm = names(data) = c(nm, paste0(".xtra.", seq_len(m - k), "."))
    }
    for (i in num[!is.na(num)]) {
        j = which(nm == as.character(i))
        if (i != j) {
            tmp = nm[i]
            nm[i] = nm[j]
            nm[j] = tmp
        }
    }
    data[, nm, drop = FALSE]
}

# format sigma for use in  messages
.fmt.sigma = function(sigma) {
    if (length(sigma) == 1)
        round(sigma, 4 - floor(log10(sigma)))
    else
        "(various values)"
}

# format a transformation for messages
.fmt.tran = function(misc) {
    tran = misc$tran
    if (is.list(tran)) 
        tran = ifelse(is.null(tran$name), "custom", tran$name)
    if (!is.null(mul <- misc$tran.mult))
        tran = paste0(mul, "*", tran)
    if(!is.null(off <- misc$tran.offset))
        tran = paste0(tran, "(mu + ", off, ")")
    tran
}


# My own utility for requiring a namespace and handling case where it is not available
#   pkg      package name
#   ...      passed to fail
#   fail     function to call if namespace not found
#   quietly passed to requireNamespace()
## I can't decide definitively if I want to suppress S3 masking messages or not...
.requireNS = function(pkg, ..., fail = stop, quietly = TRUE) {
    ### result = suppressMessages(requireNamespace(pkg, quietly = quietly))
    result = requireNamespace(pkg, quietly = TRUE)
    if (!result) fail(...)
    result
}
# of possible use as fail in .requireNS
.nothing = function(...) invisible()


### Utilities for converting symm matrices to and from lower-triangle storage mode
.get.lt = function(X) {
    rtn = X[lower.tri(X, diag = TRUE)]
    attr(rtn, "nrow") = nrow(X)
    rtn
}

.lt2mat = function(lt) {
    if (is.null(n <- attr(lt, "nrow")))
        n = (sqrt(8 * length(lt) + 1) - 1) / 2
    X = matrix(NA, nrow = n, ncol = n)
    lti = which(lower.tri(X, diag = TRUE))
    X[lti] = lt
    X = t(X)
    X[lti] = lt
    X
}

### submodel support...

# Compact a model matrix
# This returns just the R part of X's QR decomposition.
# This is sufficient in lieu of the whole model matrix
# Ideally, X is already a qr object in which case weights is assumed already incorporated
# assign should be right if input is generated by model.matrix() or is $qr slot of lm
#' @rdname extending-emmeans
#' @order 33
#' @param X,weights,assign Arguments for \code{.cmpMM}, which compacts a model
#'   matrix \code{X} into a much smaller matrix that has the same row space.
#'   Specifically, it returns the R portion of its QR decomposition. If \code{X}
#'   is already of class \code{qr}, it is used directly. \code{weights} should be
#'   the weights used in the model fit, and \code{assign} is used for unravelling
#'   any pivoting done by \code{\link{qr}}.
#' @export
.cmpMM = function(X, weights = rep(1, nrow(X)), assign = attr(X$qr, "assign")) {
    if (!get_emm_option("enable.submodel"))
        return(NULL)
    if(!is.qr(X)) {
        if(any(is.na(X)))
            return(NULL)
        X = try({
            X = sweep(X, 1, sqrt(weights), "*") # scale rows by sqrt(weights)
            qr(X)
        }, silent = TRUE)
        if (inherits(X, "try-error"))
            return(NULL)
    }
    R = qr.R(X, complete = FALSE)
    R[, X$pivot] = R
    colnames(R)[X$pivot] = colnames(R)
    attr(R, "assign") = assign
    R
}

# Get 'factors' table and simplify it (remove function calls from term labels)
.smpFT = function(trms) {
    tbl = attr(trms, "factors")
    rn = rownames(tbl)
    newrn = sapply(rn, function(x) all.vars(as.formula(paste("~", x)))[1])
    if (any(newrn != rn)) {
        rownames(tbl) = newrn
        colnames(tbl) = apply(tbl, 2, function(x) paste(newrn[x > 0], collapse = ":"))
    }
    tbl
}

# Alias matrix. Goal here is to find indexes
#    i1 = indices of columns of smaller model
#    i2 = indices of all other terms
# Then for a given linear hypothesis L_R %*% bhat_R (_R subscripts smaller model)
#    we have L_R %*% bhat_R = L_R %*% (A %*% bhat_F) = (L_R %*% A) %*% bhat_F
#    (where _F subscripts full model).
# Moreover, L_R is just columns i1 of the linfct for the effect of interest.
.alias.matrix = function(object, submodel) {
    X = object@model.info$model.matrix  # assumed to have attributes "factors" and "assign"
    if (is.character(X)) {  # model.matrix is a message
        if (nchar(X) > 0) warning(X, call. = FALSE)
        return(NULL)
    }
    assign = attr(X, "assign")
    if (is.null(assign)) {  ### either missing model matrix or assign attribute
        warning("submodel information is not available for this object", call. = FALSE)
        return(NULL)
    }
    tbl = attr(X, "factors")
    if (is.character(submodel)) {
        type2 = pmatch(submodel[1], "type2", nomatch = 0) # 1 if type2, 0 otherwise
        submodel = as.formula(paste("~", paste(names(object@levels), collapse="*")))
    }
    else
        type2 = 0
    # now submodel is a formula
    # create term labels in compatible factor order
    subtbl = attr(terms(update(object@model.info$terms, submodel)), "factors")
    com = intersect(rownames(tbl), rownames(subtbl))
    if(length(com) == 0) { # No matching factors at all
        com = 1; 
        subtbl = matrix(0)
    }
    sublab = apply(subtbl[com, , drop = FALSE], 2, function(x) paste(com[x > 0], collapse = ":"))
    usecols = intersect(colnames(tbl), sublab)
    incl = c(0, which(colnames(tbl) %in% usecols))
    if (type2) {
        incl = max(incl)  # just the last one
        overlap = apply(tbl, 2, function(x) sum(x * tbl[, incl]))
        zaplist = which(overlap < sum(tbl[, incl])) # don't contain our term
        i0 = which(assign %in% c(0, zaplist))
    }
    else
        i0 = integer(0)
    rcols = which(assign %in% incl)
    if (length(i0) > 0) {
        X[, rcols] = qr.resid(qr(X[, i0, drop = FALSE]), X[, rcols, drop = FALSE])
        X[, i0] = 0
    }
    A = qr.coef(qr(X[, rcols, drop = FALSE]), X)  # alias matrix
    A[is.na(A)] = 0  ## NA coefs are really ones constrained to zero
    attr(A, "rcols") = rcols   # cols in submodel
    attr(A, "submodstr") = paste(colnames(tbl)[incl], collapse = " + ")
    A
}
