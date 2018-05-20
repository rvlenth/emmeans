### Rudimentary support for brms. 
### Obviously this is way less than is needed, but it does support simpler models

#xxxx' @importFrom brms parse_bf
recover_data.brmsfit = function(object, data, ...) {
    bt = brms::parse_bf(formula(object))
    if (class(bt) != "brmsterms")
        stop("This model is currently not supported.")
    form = bt$dpars$mu$fe
    # we don't have a call component so I'll just put in NULL
    recover_data.call(NULL, terms(form), "na.omit", data = object$data, ...)
}

emm_basis.brmsfit = function(object, trms, xlev, grid, vcov., ...) {
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    contr = lapply(object$data, function(.) attr(., "contrasts"))
    contr = contr[!sapply(contr, is.null)]
    X = model.matrix(trms, m, contrasts.arg = contr)
    nm = gsub("(Intercept)", "Intercept", dimnames(X)[[2]], fixed = TRUE)
    V = vcov(object)[nm, nm, drop = FALSE]
    nbasis = estimability::all.estble
    dfargs = list()
    dffun = function(k, dfargs) Inf
    misc = .std.link.labels(brms::parse_bf(formula(object))$dpars$mu$family, list())
    post.beta = as.matrix(object, pars = paste0("b_", nm), exact = TRUE)
    bhat = apply(post.beta, 2, mean)
    
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, 
         misc=misc, post.beta=post.beta)
}
