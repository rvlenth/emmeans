# Estimated marginal means of linear trends

The `emtrends` function is useful when a fitted model involves a
numerical predictor \\x\\ interacting with another predictor `a`
(typically a factor). Such models specify that \\x\\ has a different
trend depending on \\a\\; thus, it may be of interest to estimate and
compare those trends. Analogous to the
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md)
setting, we construct a reference grid of these predicted trends, and
then possibly average them over some of the predictors in the grid.

## Usage

``` r
emtrends(object, specs, var, delta.var = 0.001 * rng, max.degree = 1, ...)
```

## Arguments

- object:

  A supported model object (*not* a reference grid)

- specs:

  Specifications for what marginal trends are desired – as in
  [`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md).
  If `specs` is missing or `NULL`, `emmeans` is not run and the
  reference grid for specified trends is returned.

- var:

  Character value giving the name of a variable with respect to which a
  difference quotient of the linear predictors is computed. In order for
  this to be useful, `var` should be a numeric predictor that interacts
  with at least one factor in `specs`. Then instead of computing EMMs,
  we compute and compare the slopes of the `var` trend over levels of
  the specified other predictor(s). As in EMMs, marginal averages are
  computed for the predictors in `specs` and `by`. See also the
  “Generalizations” section below.

- delta.var:

  The value of *h* to use in forming the difference quotient \\(f(x+h) -
  f(x))/h\\. Changing it (especially changing its sign) may be necessary
  to avoid numerical problems such as logs of negative numbers. The
  default value is 1/1000 of the range of `var` over the dataset.

- max.degree:

  Integer value. The maximum degree of trends to compute (this is capped
  at 5). If greater than 1, an additional factor `degree` is added to
  the grid, with corresponding numerical derivatives of orders
  `1, 2, ..., max.degree` as the estimates.

- ...:

  Additional arguments passed to
  [`ref_grid`](https://rvlenth.github.io/emmeans/reference/ref_grid.md)
  or [`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md)
  as appropriate. See Details.

## Value

An `emmGrid` or `emm_list` object, according to `specs`. See
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md) for
more details on when a list is returned.

## Details

The function works by constructing reference grids for `object` with
various values of `var`, and then calculating difference quotients of
predictions from those reference grids. Finally,
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md) is
called with the given `specs`, thus computing marginal averages as
needed of the difference quotients. Any `...` arguments are passed to
the `ref_grid` and
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md);
examples of such optional arguments include optional arguments (often
`mode`) that apply to specific models; `ref_grid` options such as
`data`, `at`, `cov.reduce`, `mult.names`, `nesting`, or `transform`; and
`emmeans` options such as `weights` (but please avoid `trend` or
`offset`.

## Note

In earlier versions of `emtrends`, the first argument was named `model`
rather than `object`. (The name was changed because of potential
mis-matching with a `mode` argument, which is an option for several
types of models.) For backward compatibility, `model` still works
*provided all arguments are named*.

It is important to understand that trends computed by `emtrends` are
*not* equivalent to polynomial contrasts in a parallel model where `var`
is regarded as a factor. That is because the model `object` here is
assumed to fit a smooth function of `var`, and the estimated trends
reflect *local* behavior at particular value(s) of `var`; whereas when
`var` is modeled as a factor and polynomial contrasts are computed,
those contrasts represent the *global* pattern of changes over *all*
levels of `var`.

See the `pigs.poly` and `pigs.fact` examples below for an illustration.
The linear and quadratic trends depend on the value of `percent`, but
the cubic trend is constant (because that is true of a cubic polynomial,
which is the underlying model). The cubic contrast in the factorial
model has the same P value as for the cubic trend, again because the
cubic trend is the same everywhere.

## Generalizations

Instead of a single predictor, the user may specify some monotone
function of one variable, e.g., `var = "log(dose)"`. If so, the chain
rule is applied. Note that, in this example, if `object` contains
`log(dose)` as a predictor, we will be comparing the slopes estimated by
that model, whereas specifying `var = "dose"` would perform a
transformation of those slopes, making the predicted trends vary
depending on `dose`.

## See also

[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md),
[`ref_grid`](https://rvlenth.github.io/emmeans/reference/ref_grid.md)

## Examples

``` r
fiber.lm <- lm(strength ~ diameter*machine, data=fiber)
# Obtain slopes for each machine ...
( fiber.emt <- emtrends(fiber.lm, "machine", var = "diameter") )
#>  machine diameter.trend    SE df lower.CL upper.CL
#>  A                1.104 0.194  9    0.666     1.54
#>  B                0.857 0.224  9    0.351     1.36
#>  C                0.864 0.208  9    0.394     1.33
#> 
#> Confidence level used: 0.95 
# ... and pairwise comparisons thereof
pairs(fiber.emt)
#>  contrast estimate    SE df t.ratio p.value
#>  A - B     0.24714 0.296  9   0.835  0.6919
#>  A - C     0.24008 0.284  9   0.845  0.6863
#>  B - C    -0.00705 0.306  9  -0.023  0.9997
#> 
#> P value adjustment: tukey method for comparing a family of 3 estimates 

# Suppose we want trends relative to sqrt(diameter)...
emtrends(fiber.lm, ~ machine | diameter, var = "sqrt(diameter)", 
         at = list(diameter = c(20, 30)))
#> diameter = 20:
#>  machine sqrt(diameter).trend   SE df lower.CL upper.CL
#>  A                       9.88 1.73  9     5.96     13.8
#>  B                       7.67 2.00  9     3.14     12.2
#>  C                       7.73 1.86  9     3.52     11.9
#> 
#> diameter = 30:
#>  machine sqrt(diameter).trend   SE df lower.CL upper.CL
#>  A                      12.10 2.12  9     7.30     16.9
#>  B                       9.39 2.45  9     3.84     14.9
#>  C                       9.47 2.28  9     4.31     14.6
#> 
#> Confidence level used: 0.95 

# Obtaining a reference grid
mtcars.lm <- lm(mpg ~ poly(disp, degree = 2) * (factor(cyl) + factor(am)), data = mtcars)

# Center trends at mean disp for each no. of cylinders
mtcTrends.rg <- emtrends(mtcars.lm, var = "disp", 
                          cov.reduce = disp ~ factor(cyl))
summary(mtcTrends.rg)  # estimated trends at grid nodes
#>  disp cyl am disp.trend     SE df
#>   105   4  0    -0.0949 0.0829 20
#>   183   6  0    -0.0024 0.0496 20
#>   353   8  0    -0.0106 0.0105 20
#>   105   4  1    -0.1212 0.0338 20
#>   183   6  1    -0.0217 0.0573 20
#>   353   8  1    -0.0147 0.0645 20
#> 
emmeans(mtcTrends.rg, "am", weights = "prop")
#>  am disp.trend     SE df lower.CL upper.CL
#>   0    -0.0378 0.0312 20   -0.103  0.02733
#>   1    -0.0529 0.0260 20   -0.107  0.00145
#> 
#> Results are averaged over the levels of: cyl 
#> Confidence level used: 0.95 


### Higher-degree trends ...

pigs.poly <- lm(conc ~ poly(percent, degree = 3), data = pigs)
emt <- emtrends(pigs.poly, ~ degree | percent, "percent", max.degree = 3,
                at = list(percent = c(9, 13.5, 18)))
       # note: 'degree' is an extra factor created by 'emtrends'
       
summary(emt, infer = c(TRUE, TRUE))
#> percent =  9.0:
#>  degree    percent.trend     SE df lower.CL upper.CL t.ratio p.value
#>  linear          2.39923 3.6500 25   -5.119    9.917   0.657  0.5170
#>  quadratic      -0.22674 1.1000 25   -2.498    2.044  -0.206  0.8387
#>  cubic           0.00548 0.0825 25   -0.164    0.175   0.066  0.9475
#> 
#> percent = 13.5:
#>  degree    percent.trend     SE df lower.CL upper.CL t.ratio p.value
#>  linear          0.69212 1.5600 25   -2.528    3.912   0.443  0.6618
#>  quadratic      -0.15277 0.1750 25   -0.513    0.207  -0.874  0.3903
#>  cubic           0.00548 0.0825 25   -0.164    0.175   0.066  0.9475
#> 
#> percent = 18.0:
#>  degree    percent.trend     SE df lower.CL upper.CL t.ratio p.value
#>  linear         -0.34928 4.1200 25   -8.830    8.131  -0.085  0.9331
#>  quadratic      -0.07880 1.1500 25   -2.448    2.291  -0.068  0.9459
#>  cubic           0.00548 0.0825 25   -0.164    0.175   0.066  0.9475
#> 
#> Confidence level used: 0.95 

# Compare above results with poly contrasts when 'percent' is modeled as a factor ...
pigs.fact <- lm(conc ~ factor(percent), data = pigs)
emm <- emmeans(pigs.fact, "percent")

contrast(emm, "poly")
#>  contrast  estimate    SE df t.ratio p.value
#>  linear      23.837 14.70 25   1.617  0.1184
#>  quadratic   -5.500  6.29 25  -0.874  0.3903
#>  cubic        0.888 13.40 25   0.066  0.9475
#> 
# Some P values are comparable, some aren't! See Note in documentation
```
