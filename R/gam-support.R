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

### Support for 'gam' objects
# Note: This is a mess, because both packages 'gam' and 'mgcv' produce these,
# and they are different. Both inherit from glm and lm, though, so recover_data.lm
# still serves for these (I hope)


# Right now I'm assuming gam package. Will add or trap-out mgcv later

# We have two args:
#   nboot   # of bootstrap reps to get variances of smooths
emm_basis.gam = function(object, trms, xlev, grid, nboot = 800, ...) {
    if (!is.null(object$gcv.ubre)) # From mgcv, not gam
        return (emm_basis.gam_mgcv(object, trms, xlev, grid, ...))
    
    result = emm_basis.lm(object, trms, xlev, grid, ...)
    old.smooth = object$smooth
    if (is.null(old.smooth))  # "just an ordinary glm" (My Fair Lady)
        return(result)
    # else we need to add-in some smoothers
    smooth.frame = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    data = object$smooth.frame
    labs = names(data)
    w = object$weights
    resid = object$residuals
    for (i in seq_along(labs)) {
        lab = labs[i]
        sig = apply(smooth.frame[[i]], 1, paste, collapse = ":")
        usig = unique(sig)
        rows = lapply(usig, function(s) which(sig == s))
        xeval = smooth.frame[sapply(rows, "[", 1), lab]
        bsel = matrix(0, nrow = length(sig), ncol = length(usig))
        for (j in seq_along(rows))
            bsel[rows[[j]], j] = 1
        
        cl = attr(data[[i]], "call")
        cl$xeval = substitute(xeval)
        z = resid + old.smooth[, lab]
        bh = as.numeric(eval(cl))
        m = length(bh)
        n = length(result$bhat)
        result$bhat = c(result$bhat, bh)
        result$X = cbind(result$X, bsel)
        boot = replicate(nboot, {
                z = sample(resid, replace = TRUE) + old.smooth[, lab]
                as.numeric(eval(cl))
            })
        covar = if(m == 1) var(boot) 
                else       cov(t(boot))
        result$V = rbind(cbind(result$V, matrix(0, nrow = n, ncol = m)),
                         cbind(matrix(0, nrow = m, ncol = n), covar))
    }
    result
}


### emm_basis method for mgcv::gam objects
emm_basis.gam_mgcv = function(object, trms, xlev, grid, ...) {
    stop("Can't handle mgcv::gam objects")
}