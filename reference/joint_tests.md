# Compute joint tests of the terms in a model

This function produces an analysis-of-variance-like table based on
linear functions of predictors in a model or `emmGrid` object.
Specifically, the function constructs, for each combination of factors
(or covariates reduced to two or more levels), a set of (interaction)
contrasts via
[`contrast`](https://rvlenth.github.io/emmeans/reference/contrast.md),
and then tests them using
[`test`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
with `joint = TRUE`. Optionally, one or more of the predictors may be
used as `by` variable(s), so that separate tables of tests are produced
for each combination of them.

## Usage

``` r
joint_tests(object, by = NULL, show0df = FALSE, showconf = TRUE,
  cov.reduce = make.meanint(1), ...)

make.meanint(delta = 1, npts = 2)

meanint(x)

make.symmint(ctr, delta = 1, npts = 2)

symmint(ctr)
```

## Arguments

- object:

  a fitted model, `emmGrid`, or `emm_list`. If the latter, its first
  element is used.

- by:

  character names of `by` variables. Separate sets of tests are run for
  each combination of these.

- show0df:

  logical value; if `TRUE`, results with zero numerator degrees of
  freedom are displayed, if `FALSE` they are skipped

- showconf:

  logical value. When we have models with estimability issues (e.g.,
  missing cells), then with `showconf = TRUE`, we test any remaining
  effects that are not purely due to contrasts of a single term. If
  found, they are labeled `(confounded)`. See
  [`vignette("xplanations")`](https://rvlenth.github.io/emmeans/articles/xplanations.md)
  for more information.

- cov.reduce:

  a function. If `object` is a fitted model, it is replaced by
  `ref_grid(object, cov.reduce = cov.reduce, ...)`. For this purpose,
  the functions `meanint` and `symmint` are available for returning an
  interval around the mean or around zero, respectively. Se the section
  below on covariates.

- ...:

  additional arguments passed to `ref_grid` and `emmeans`

- delta, ctr:

  arguments for `make.meanint` and `make.symmint`. `delta` sets the
  distance each side of the center, so that the width of the interval is
  `2*delta`.

- npts:

  number of points to include in the interval

- x:

  argument for `meanint` and `symmint`

## Value

a `summary_emm` object (same as is produced by
[`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)).
All effects for which there are no estimable contrasts are omitted from
the results. There may be an additional row named `(confounded)` which
accounts for additional degrees of freedom for effects not accounted for
in the preceding rows.

The returned object also includes an `"est.fcns"` attribute, which is a
named list containing the linear functions associated with each joint
test. Each row of these is standardized to have length 1. No estimable
functions are included for confounded effects.

`make.meanint` returns the function
`function(x) mean(x) + delta * c(-1, 1)`, and `make.symmint(ctr, delta)`
returns the function `function(x) ctr + delta * c(-1, 1)` (which does
not depend on `x`). The cases with `delta = 1`,
`meanint = make.meanint(1)` and `symmint(ctr) = make.symmint(ctr, 1)`
are retained for back-compatibility reasons. These functions are
available primarily for use with `cov.reduce`.

## Details

In models with only factors, no covariates, these tests correspond to
“type III” tests a la SAS, as long as equal-weighted averaging is used
and there are no estimability issues. When covariates are present and
they interact with factors, the results depend on how the covariate is
handled in constructing the reference grid. See the section on
covariates below. The point that one must always remember is that
`joint_tests` always tests contrasts among EMMs, in the context of the
reference grid, whereas SAS's type III tests are tests of model
coefficients – which may or may not have anything to do with EMMs or
contrasts.

## Note

`joint_tests` is flaky with models having nested fixed effects. In some
cases, terms that could be relevant are not identified, or confounded
with unidentifiable terms.

## Dealing with covariates

A covariate (or any other predictor) must have *more than one value in
the reference grid* in order to test its effect and be included in the
results. Therefore, when `object` is a model, we default to
`cov.reduce = meanint` which sets each covariate at a symmetric interval
about its mean. But when `object` is an existing reference grid, it
often has only one value for covariates, in which case they are excluded
from the joint tests.

While having two points is sufficient when the covariate term has a
linear trend, you need more than two when some kind of curved trend
(polynomial, spline, etc.) is present – else `joint_tests()` will not
show enough degrees of freedom for terms involving the covariate. You
may specify these points manually using `at`, or by including an `npts`
argument in `cov.reduce`, via `make.meanint` or `make.symmint()`. With
some kinds of curved trends, the joint tests of covariate terms may
become somewhat meaningless.

Covariates present further complications in that their values in the
reference grid can affect the joint tests of *other* effects. When
covariates are centered around their means (the default), then the tests
we obtain can be described as joint tests of covariate-adjusted means;
and that is our intended use here. However, some software such as SAS
and [`car::Anova`](https://rdrr.io/pkg/car/man/Anova.html) adopt the
convention of centering covariates around zero; and for that purpose,
one can use `cov.reduce = symmint(0)` when calling with a model object
(or in constructing a reference grid). However, adjusted means with
covariates set at or around zero do not make much sense in the context
of interpreting estimated marginal means, unless the covariate means
really are zero.

See the examples below with the `toy` dataset.

## See also

[`test`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)

## Examples

``` r
pigs.lm <- lm(log(conc) ~ source * factor(percent), data = pigs)

(jt <- joint_tests(pigs.lm))             ## will be same as type III ANOVA
#>  model term     df1 df2 F.ratio p.value
#>  source           2  17  30.256 <0.0001
#>  percent          3  17   8.214  0.0013
#>  source:percent   6  17   0.926  0.5011
#> 

### Estimable functions associated with "percent"
attr(jt, "est.fcns") $ "percent"
#>      (Intercept) sourcesoy sourceskim factor(percent)12 factor(percent)15
#> [1,]           0         0          0         -0.904534          0.000000
#> [2,]           0         0          0          0.000000         -0.904534
#> [3,]           0         0          0          0.000000          0.000000
#>      factor(percent)18 sourcesoy:factor(percent)12 sourceskim:factor(percent)12
#> [1,]          0.000000                  -0.3015113                   -0.3015113
#> [2,]          0.000000                   0.0000000                    0.0000000
#> [3,]         -0.904534                   0.0000000                    0.0000000
#>      sourcesoy:factor(percent)15 sourceskim:factor(percent)15
#> [1,]                   0.0000000                    0.0000000
#> [2,]                  -0.3015113                   -0.3015113
#> [3,]                   0.0000000                    0.0000000
#>      sourcesoy:factor(percent)18 sourceskim:factor(percent)18
#> [1,]                   0.0000000                    0.0000000
#> [2,]                   0.0000000                    0.0000000
#> [3,]                  -0.3015113                   -0.3015113

joint_tests(pigs.lm, weights = "outer")  ## differently weighted
#>  model term     df1 df2 F.ratio p.value
#>  source           2  17  28.601 <0.0001
#>  percent          3  17   7.889  0.0016
#>  source:percent   6  17   0.926  0.5011
#> 

joint_tests(pigs.lm, by = "source")      ## separate joint tests of 'percent'
#> source = fish:
#>  model term df1 df2 F.ratio p.value
#>  percent      3  17   1.712  0.2023
#> 
#> source = soy:
#>  model term df1 df2 F.ratio p.value
#>  percent      3  17   1.290  0.3097
#> 
#> source = skim:
#>  model term df1 df2 F.ratio p.value
#>  percent      3  17   6.676  0.0035
#> 

### Comparisons with type III tests in SAS
toy = data.frame(
    treat = rep(c("A", "B"), c(4, 6)),
    female = c(1, 0, 0, 1,   0, 0, 0, 1, 1, 0 ),
    resp = c(17, 12, 14, 19, 28, 26, 26, 34, 33, 27))
toy.fac = lm(resp ~ treat * factor(female), data = toy)
toy.cov = lm(resp ~ treat * female, data = toy)
# (These two models have identical fitted values and residuals)

# -- SAS output we'd get with toy.fac --
## Source          DF    Type III SS    Mean Square   F Value   Pr > F
## treat            1    488.8928571    488.8928571    404.60   <.0001
## female           1     78.8928571     78.8928571     65.29   0.0002
## treat*female     1      1.7500000      1.7500000      1.45   0.2741
# 
# -- SAS output we'd get with toy.cov --
## Source          DF    Type III SS    Mean Square   F Value   Pr > F
## treat            1    252.0833333    252.0833333    208.62   <.0001
## female           1     78.8928571     78.8928571     65.29   0.0002
## female*treat     1      1.7500000      1.7500000      1.45   0.2741

joint_tests(toy.fac)
#>  model term   df1 df2 F.ratio p.value
#>  treat          1   6 404.601 <0.0001
#>  female         1   6  65.291  0.0002
#>  treat:female   1   6   1.448  0.2741
#> 
joint_tests(toy.cov)   # female is regarded as a 2-level factor by default
#>  model term   df1 df2 F.ratio p.value
#>  treat          1   6 404.601 <0.0001
#>  female         1   6  65.291  0.0002
#>  treat:female   1   6   1.448  0.2741
#> 

## Treat 'female' as a numeric covariate (via cov.keep = 0)
## ... then tests depend on where we center things

# Center around the mean
joint_tests(toy.cov, cov.keep = 0, cov.reduce = make.meanint(delta = 1))
#>  model term   df1 df2 F.ratio p.value
#>  treat          1   6 401.865 <0.0001
#>  female         1   6  65.291  0.0002
#>  treat:female   1   6   1.448  0.2741
#> 
# Center around zero (like SAS's results for toy.cov)
joint_tests(toy.cov, cov.keep = 0, cov.reduce = make.symmint(ctr = 0, delta = 1))
#>  model term   df1 df2 F.ratio p.value
#>  treat          1   6 208.621 <0.0001
#>  female         1   6  65.291  0.0002
#>  treat:female   1   6   1.448  0.2741
#> 
# Center around 0.5 (like SAS's results for toy.fac)
joint_tests(toy.cov, cov.keep = 0, cov.reduce = range)
#>  model term   df1 df2 F.ratio p.value
#>  treat          1   6 404.601 <0.0001
#>  female         1   6  65.291  0.0002
#>  treat:female   1   6   1.448  0.2741
#> 

### Example with empty cells and confounded effects
low3 <- unlist(attr(ubds, "cells")[1:3]) 
ubds.lm <- lm(y ~ A*B*C, data = ubds, subset = -low3)

# Show overall joint tests by C:
ref_grid(ubds.lm, by = "C") |> contrast("consec") |> test(joint = TRUE)
#>  C df1 df2 F.ratio p.value note
#>  1   6  71   8.993 <0.0001    e
#>  2   7  71   9.763 <0.0001    e
#>  3   8  71   8.774 <0.0001     
#> 
#> e: df1 reduced due to non-estimability 

# Break each of the above into smaller components:
joint_tests(ubds.lm, by = "C")
#> C = 1:
#>  model term   df1 df2 F.ratio p.value note
#>  A:B            2  71   5.709  0.0050    e
#>  (confounded)   4  71  10.635 <0.0001     
#> 
#> C = 2:
#>  model term   df1 df2 F.ratio p.value note
#>  A              1  71  14.601  0.0003    e
#>  B              1  71  24.059 <0.0001    e
#>  A:B            3  71   2.778  0.0474    e
#>  (confounded)   2  71   8.656  0.0004     
#> 
#> C = 3:
#>  model term   df1 df2 F.ratio p.value note
#>  A              2  71   7.398  0.0012     
#>  B              2  71  22.157 <0.0001     
#>  A:B            4  71   3.009  0.0237     
#> 
#> e: df1 reduced due to non-estimability 
```
