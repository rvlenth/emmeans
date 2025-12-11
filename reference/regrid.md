# Reconstruct a reference grid with a new transformation or simulations

The typical use of this function is to cause EMMs to be computed on a
different scale, e.g., the back-transformed scale rather than the
linear-predictor scale. In other words, if you want back-transformed
results, do you want to average and then back-transform, or
back-transform and then average?

## Usage

``` r
regrid(object, transform = c("response", "mu", "unlink", "none", "pass",
  links), inv.link.lbl = "response", predict.type,
  bias.adjust = get_emm_option("back.bias.adj"), sigma, N.sim,
  sim = mvtnorm::rmvnorm, ...)
```

## Arguments

- object:

  An object of class `emmGrid`

- transform:

  Character, list, or logical value. If `"response"`, `"mu"`, or `TRUE`,
  the inverse transformation is applied to the estimates in the grid
  (but if there is both a link function and a response transformation,
  `"mu"` back-transforms only the link part); if `"none"` or `FALSE`,
  `object` is re-gridded so that its `bhat` slot contains
  `predict(object)` and its `linfct` slot is the identity. Any internal
  transformation information is preserved. If `transform = "pass"`, the
  object is not re-gridded in any way (this may be useful in conjunction
  with `N.sim`).

  If `transform` is a character value in `links` (which is the set of
  valid arguments for the
  [`make.link`](https://rdrr.io/r/stats/make.link.html) function,
  excepting `"identity"`), or if `transform` is a list of the same form
  as returned by `make.links` or
  [`make.tran`](https://rvlenth.github.io/emmeans/reference/make.tran.md),
  the results are formulated as if the response had been transformed
  with that link function.

- inv.link.lbl:

  Character value. This applies only when `transform` is in `links`, and
  is used to label the predictions if subsequently summarized with
  `type = "response"`.

- predict.type:

  Character value. If provided, the returned object is updated with the
  given type to use by default by `summary.emmGrid` (see
  [`update.emmGrid`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md)).
  This may be useful if, for example, when one specifies
  `transform = "log"` but desires summaries to be produced by default on
  the response scale.

- bias.adjust:

  Logical value for whether to adjust for bias in back-transforming
  (`transform = "response"`). This requires a valid value of `sigma` to
  exist in the object or be specified.

- sigma:

  Error SD assumed for bias correction (when `transform = "response"`
  and a transformation is in effect). If not specified,
  `object@misc$sigma` is used, and a warning is issued if it is not
  found.

- N.sim:

  Integer value. If specified and `object` is based on a frequentist
  model (i.e., does not have a posterior sample), then a fake posterior
  sample is generated using the function `sim`.

- sim:

  A function of three arguments (no names are assumed). If `N.sim` is
  supplied with a frequentist model, this function is called with
  respective arguments `N.sim`, `object@bhat`, and `object@V`. The
  default is the multivariate normal distribution.

- ...:

  Ignored.

## Value

An `emmGrid` object with the requested changes

## Details

The `regrid` function reparameterizes an existing `ref.grid` so that its
`linfct` slot is the identity matrix and its `bhat` slot consists of the
estimates at the grid points. If `transform` is `TRUE`, the inverse
transform is applied to the estimates. Outwardly, when
`transform = "response"`, the result of
[`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
after applying `regrid` is identical to the summary of the original
object using `type="response"`. But subsequent EMMs or contrasts will be
conducted on the new scale – which is the reason this function exists.

This function may also be used to simulate a sample of regression
coefficients for a frequentist model for subsequent use as though it
were a Bayesian model. To do so, specify a value for `N.sim` and a
sample is simulated using the function `sim`. The grid may be further
processed in accordance with the other arguments; or if
`transform = "pass"`, it is simply returned with the only change being
the addition of the simulated sample.

## Note

Another way to use `regrid` is to supply a `regrid` argument to
[`ref_grid`](https://rvlenth.github.io/emmeans/reference/ref_grid.md)
(either directly of indirectly via
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md)), in
which case its value is passed to `regrid` as `transform`. This is often
a simpler approach if the reference grid has not already been
constructed.

## Degrees of freedom

In cases where the degrees of freedom depended on the linear function
being estimated (e.g., Satterthwaite method), the d.f. from the
reference grid are saved, and a kind of “containment” method is
substituted in the returned object, whereby the calculated d.f. for a
new linear function will be the minimum d.f. among those having nonzero
coefficients. This is kind of an *ad hoc* method, and it can
over-estimate the degrees of freedom in some cases. An annotation is
displayed below any subsequent summary results stating that the
degrees-of-freedom method is inherited from the previous method at the
time of re-gridding.

## Examples

``` r
pigs.lm <- lm(log(conc) ~ source + factor(percent), data = pigs)
rg <- ref_grid(pigs.lm)

# This will yield EMMs as GEOMETRIC means of concentrations:
(emm1 <- emmeans(rg, "source", type = "response"))
#>  source response   SE df lower.CL upper.CL
#>  fish       29.8 1.09 23     27.6     32.1
#>  soy        39.1 1.47 23     36.2     42.3
#>  skim       44.6 1.75 23     41.1     48.3
#> 
#> Results are averaged over the levels of: percent 
#> Confidence level used: 0.95 
#> Intervals are back-transformed from the log scale 
pairs(emm1) ## We obtain RATIOS
#>  contrast    ratio     SE df null t.ratio p.value
#>  fish / soy  0.761 0.0403 23    1  -5.153 <0.0001
#>  fish / skim 0.669 0.0362 23    1  -7.428 <0.0001
#>  soy / skim  0.879 0.0466 23    1  -2.442  0.0570
#> 
#> Results are averaged over the levels of: percent 
#> P value adjustment: tukey method for comparing a family of 3 estimates 
#> Tests are performed on the log scale 

# This will yield EMMs as ARITHMETIC means of concentrations:
(emm2 <- emmeans(regrid(rg, transform = "response"), "source"))
#>  source response   SE df lower.CL upper.CL
#>  fish       30.0 1.10 23     27.7     32.2
#>  soy        39.4 1.49 23     36.3     42.5
#>  skim       44.8 1.79 23     41.1     48.5
#> 
#> Results are averaged over the levels of: percent 
#> Confidence level used: 0.95 
pairs(emm2)  ## We obtain DIFFERENCES
#>  contrast    estimate   SE df t.ratio p.value
#>  fish - soy     -9.40 1.86 23  -5.051  0.0001
#>  fish - skim   -14.84 2.10 23  -7.071 <0.0001
#>  soy - skim     -5.44 2.25 23  -2.424  0.0591
#> 
#> Results are averaged over the levels of: percent 
#> P value adjustment: tukey method for comparing a family of 3 estimates 
# Same result, useful if we hadn't already created 'rg'
# emm2 <- emmeans(pigs.lm, "source", regrid = "response")

# Simulate a sample of regression coefficients
set.seed(2.71828)
rgb <- regrid(rg, N.sim = 200, transform = "pass")
#> Simulating a sample of size 200 of regression coefficients.
emmeans(rgb, "source", type = "response")  ## similar to emm1
#>  source response lower.HPD upper.HPD
#>  fish       29.8      27.7      31.9
#>  soy        39.3      36.5      41.8
#>  skim       44.7      41.7      48.3
#> 
#> Results are averaged over the levels of: percent 
#> Point estimate displayed: median 
#> Results are back-transformed from the log scale 
#> HPD interval probability: 0.95 
```
