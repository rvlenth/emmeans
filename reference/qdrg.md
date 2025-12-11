# Quick and dirty reference grid

This function may make it possible to compute a reference grid for a
model object that is otherwise not supported.

## Usage

``` r
qdrg(formula, data, coef, vcov, df, mcmc, object, subset, weights, contrasts,
  link, qr, ordinal, ...)

# Default S3 method
qdrg(formula = stats::formula(object),
  data = try(recover_data.lm(object), silent = TRUE),
  coef = stats::coef(object), vcov = stats::vcov(object),
  df = stats::df.residual(object), mcmc, object, subset,
  weights = stats::weights(object), contrasts = object$contrasts,
  link = ifelse(!is.null(lnk <- object$family$link), lnk, object$link),
  qr = object$qr, ordinal, ...)
```

## Arguments

- formula:

  Formula for the fixed effects

- data:

  Dataset containing the variables in the model

- coef:

  Fixed-effect regression coefficients (must conform to formula)

- vcov:

  Variance-covariance matrix of the fixed effects

- df:

  Error degrees of freedom

- mcmc:

  Posterior sample of fixed-effect coefficients

- object:

  Optional model object. *This rarely works!*; but if provided, we try
  to set other arguments based on an expectation that \`object\` has a
  similar structure to \`lm\` objects. See Details.

- subset:

  Subset of `data` used in fitting the model

- weights:

  Weights used in fitting the model

- contrasts:

  List of contrasts specified in fitting the model

- link:

  Link function (character or list) used, if a generalized linear model.
  (Note: response transformations are auto-detected from `formula`)

- qr:

  QR decomposition of the model matrix; used only if there are `NA`s in
  `coef`.

- ordinal:

  `list` with elements `dim` and `mode`. `ordinal$dim` (integer) is the
  number of levels in an ordinal response. If `ordinal` is provided, the
  intercept terms are modified appropriate to predicting an ordinal
  response, as described in
  [`vignette("models")`](https://rvlenth.github.io/emmeans/articles/models.md),
  Section O, using `ordinal$mode` as the `mode` argument (if not
  provided, `"latent"` is assumed). (All modes are supported except
  \`scale\`) For this to work, we expect the first `ordinal$dim - 1`
  elements of `coef` to be the estimated threshold parameters, followed
  by the coefficients for the linear predictor. Also, if `mode` requires
  back-transforming (e.g., `"prob"` or `"mean.class"`), the user may
  need to supply the `link` for it to work correctly.

- ...:

  Optional arguments passed to
  [`ref_grid`](https://rvlenth.github.io/emmeans/reference/ref_grid.md)

## Value

An `emmGrid` object constructed from the arguments

## Details

Usually, you need to provide either `object`; or `formula`, `coef`,
`vcov`, `data`, and perhaps other parameters. It is often fairly
straightforward to figure out how to get these from the model `object`;
see the documentation for the model class that was fitted. Sometimes one
or more of these quantities contains extra parameters, and if so, you
may need to subset them to make everything conformable. For a given
`formula` and `data`, you can find out what is needed via
`colnames(model.matrix(formula, data))`. (However, for an ordinal model,
we expect the first `ordinal.dim - 1` coefficients to replace
`(Intercept)`. And for a multivariate model, we expect `coef` to be a
matrix with these row names, and `vcov` to have as many rows and columns
as the total number of elements of `coef`.)

If `object` is specified, this function serves as a generic for
dispatching a method for `object`'s class. See the arguments for
`qdrg.default` to see how the arguments are determined by default. For
many `lm`- or `glm`-like models, it may suffice override just one or two
of these defaults. The possibility of a custom `qdrg` method also
provides a minimal way for package developers to provide emmeans
support: it doesn't allow directly applying
[`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.md) on
the model, but at least the user can obtain a reference grid, and then
go from there. Note that when using a \`qdrg\` method, it is best for
the user to specify \`object =\` explicitly in the call, since
\`object\` is not the first argument of \`qdrg()\`. However, if the
first argument is not a formula, the function is retried with it as
`object`.

The functions `qdrg` and `emmobj` are close cousins, in that they both
produce `emmGrid` objects. When starting with summary statistics for an
existing grid, `emmobj` is more useful, while `qdrg` is more useful when
starting from a fitted model.

## Note

For backwards compatibility, an argument `ordinal.dim` is invisibly
supported as part of `...`, and if present, sets
`ordinal = list(dim = ordinal.dim, mode = "latent")`

## Rank deficiencies

Different model-fitting packages take different approaches when the
model matrix is singular, but `qdrg` tries to reconcile them by
comparing the linear functions created by `formula` to `coefs` and
`vcov`. We may then use the estimability package to determine what
quantities are estimable. For reconciling to work properly, `coef`
should be named and `vcov` should have dimnames. To disable this
name-matching action, remove the names from `coef`, e.g., by calling
[`unname()`](https://rdrr.io/r/base/unname.html). No reconciliation is
attempted in multivariate-response cases. For more details on
estimability, see the documentation in the estimability package.

## See also

[`emmobj`](https://rvlenth.github.io/emmeans/reference/emmobj.md) for an
alternative way to construct an `emmGrid`.

## Examples

``` r
# In these examples, use emm_example(..., list = TRUE) # to see just the code

if (require(biglm, quietly = TRUE)) 
    emm_example("qdrg-biglm")
#> 
#> --- Running code from 'system.file("extexamples", "qdrg-biglm.R", package = "emmeans")'
#> 
#> > bigmod <- biglm(log(conc) ~ source + factor(percent), 
#> +     data = pigs)
#> 
#> > rg1 <- qdrg(log(conc) ~ source + factor(percent), 
#> +     data = pigs, coef = coef(bigmod), vcov = vcov(bigmod), df = bigmod$df.residual)
#> 
#> > emmeans(rg1, "source", type = "response")
#>  source response   SE  df asymp.LCL asymp.UCL
#>  fish       29.8 1.09 Inf      27.7      32.0
#>  soy        39.1 1.47 Inf      36.4      42.1
#>  skim       44.6 1.75 Inf      41.2      48.1
#> 
#> Results are averaged over the levels of: percent 
#> Confidence level used: 0.95 
#> Intervals are back-transformed from the log scale 
#> 
    
if(require(coda, quietly = TRUE) && require(lme4, quietly = TRUE)) 
    emm_example("qdrg-coda")
#> 
#> Attaching package: ‘lme4’
#> The following object is masked from ‘package:nlme’:
#> 
#>     lmList
#> 
#> --- Running code from 'system.file("extexamples", "qdrg-coda.R", package = "emmeans")'
#> 
#> > post <- readRDS(system.file("extdata", "cbpplist", 
#> +     package = "emmeans"))$post.beta
#> 
#> > rg2 <- qdrg(~size + period, data = lme4::cbpp, mcmc = post, 
#> +     link = "logit")
#> 
#> > summary(rg2, type = "response")
#>  size period response lower.HPD upper.HPD
#>    15 1        0.1930    0.1214     0.288
#>    15 2        0.0836    0.0398     0.137
#>    15 3        0.0748    0.0369     0.129
#>    15 4        0.0489    0.0138     0.101
#> 
#> Point estimate displayed: median 
#> Results are back-transformed from the logit scale 
#> HPD interval probability: 0.95 
#> 
    
if(require(ordinal, quietly = TRUE)) 
    emm_example("qdrg-ordinal")
#> 
#> --- Running code from 'system.file("extexamples", "qdrg-ordinal.R", package = "emmeans")'
#> 
#> > wine.clm <- clm(rating ~ temp * contact, data = wine)
#> 
#> > ref_grid(wine.clm)
#> 'emmGrid' object with variables:
#>     temp = cold, warm
#>     contact = no, yes
#> 
#> > qdrg(object = wine.clm, ordinal.dim = 5)
#> 'emmGrid' object with variables:
#>     temp = cold, warm
#>     contact = no, yes
#> Transformation: “logit” 
#> 
```
