### Rudimentary support for brms. 
### Obviously this is way less than is needed, but it does support simpler models


recover_data.brmsfit = function(object, data, ...) {
    fcall = object$call
    form = brms::parse_bf(formula(object))$dpars$mu$fe
    recover_data.call(fcall, terms(form), "na.omit", data = object$data, ...)
}

emm_basis.brmsfit = function(object, trms, xlev, grid, vcov., ...) {
    V = vcov(object)
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    contr = lapply(object$data, function(.) attr(., "contrasts"))
    contr = contr[!sapply(contr, is.null)]
    X = model.matrix(trms, m, contrasts.arg = contr)
    nm = dimnames(X)[[2]]
    nbasis=estimability::all.estble
    dfargs = list()
    dffun = function(k, dfargs) Inf
    misc = .std.link.labels(parse_bf(formula(object))$dpars$mu$family, list())
    # Pick the 1st length(nm) columns -- probably not the best approach...
    post.beta = as.matrix(posterior_samples(object)[, seq_along(nm), drop = FALSE])
    bhat = apply(post.beta, 2, mean)
    
    list(X=X, bhat=bhat, nbasis=nbasis, V=V, dffun=dffun, dfargs=dfargs, 
         misc=misc, post.beta=post.beta)
}
