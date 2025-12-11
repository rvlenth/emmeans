# Response-transformation extensions

The `make.tran` function creates the needed information to perform
transformations of the response variable, including inverting the
transformation and estimating variances of back-transformed predictions
via the delta method. `make.tran` is similar to
[`make.link`](https://rdrr.io/r/stats/make.link.html), but it covers
additional transformations. The result can be used as an environment in
which the model is fitted, or as the `tran` argument in
[`update.emmGrid`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md)
(when the given transformation was already applied in an existing
model).

## Usage

``` r
make.tran(type = c("genlog", "power", "boxcox", "sympower", "asin.sqrt",
  "atanh", "bcnPower", "scale"), alpha = 1, beta = 0, param, y, inner, ...)

inverse(y)
```

## Arguments

- type:

  The name of a standard transformation supported by `stat::make.link`,
  or of a special transformation described under Details.

- alpha, beta:

  Numeric parameters needed for special transformations.

- param:

  If non-missing, this specifies either `alpha` or `c(alpha, beta)`
  (provided for backward compatibility). Also, for the same reason, if
  `alpha` is of length more than 1, it is taken as `param`.

- y:

  A numeric response variable used (*and required*) with
  `type = "scale"`, where `scale(y)` determines `alpha` and `beta`.

- inner:

  another transformation. See the section on compound transformations

- ...:

  Additional arguments passed to other functions/methods

## Value

A `list` having at least the same elements as those returned by
[`make.link`](https://rdrr.io/r/stats/make.link.html). The `linkfun`
component is the transformation itself. Each of the functions is
associated with an environment where any parameter values are defined.

`inverse` returns the reciprocal of its argument. It allows the
`"inverse"` link to be auto-detected as a response transformation.

## Note

The `genlog` transformation is technically unneeded, because a response
transformation of the form `log(y + c)` is now auto-detected by
[`ref_grid`](https://rvlenth.github.io/emmeans/reference/ref_grid.md).

We modify certain [`make.link`](https://rdrr.io/r/stats/make.link.html)
results in transformations where there is a restriction on valid
prediction values, so that reasonable inverse predictions are obtained,
no matter what. For example, if a `sqrt` transformation was used but a
predicted value is negative, the inverse transformation is zero rather
than the square of the prediction. A side effect of this is that it is
possible for one or both confidence limits, or even a standard error, to
be zero.

## Details

The `make.tran` function returns a suitable list of functions for
several popular transformations. Besides being usable with `update`, the
user may use this list as an enclosing environment in fitting the model
itself, in which case the transformation is auto-detected when the
special name `linkfun` (the transformation itself) is used as the
response transformation in the call. See the examples below.

The primary purpose of `make.tran` is to support transformations that
require additional parameters, specified as `alpha` and `beta`; these
are the onse shown in the argument-matching list. However, standard
transformations supported by
[`stats::make.link`](https://rdrr.io/r/stats/make.link.html) are also
supported. In the following discussion of ones requiring parameters, we
use \\\alpha\\ and \\\beta\\ to denote `alpha` and `beta`, and \\y\\ to
denote the response variable. The `type` argument specifies the
following transformations:

- `"genlog"`:

  Generalized logarithmic transformation: \\\log\_\beta(y + \alpha)\\,
  where \\y \> -\alpha\\. When \\\beta = 0\\ (the default), we use
  \\\log_e(y + \alpha)\\

- `"power"`:

  Power transformation: \\(y-\beta)^\alpha\\, where \\y \> \beta\\. When
  \\\alpha = 0\\, \\\log(y-\beta)\\ is used instead.

- `"boxcox"`:

  The Box-Cox transformation (unscaled by the geometric mean): \\((y -
  \beta)^\alpha - 1) / \alpha\\, where \\y \> \beta\\. When \\\alpha =
  0\\, \\\log(y - \beta)\\ is used.

- `"sympower"`:

  A symmetrized power transformation on the whole real line: \\\|y -
  \beta\|^\alpha\cdot sign(y - \beta)\\. There are no restrictions on
  \\y\\, but we require \\\alpha \> 0\\ in order for the transformation
  to be monotone and continuous.

- `"asin.sqrt"`:

  Arcsin-square-root transformation: \\\sin^{-1}(y/\alpha)^{1/2}\\.
  Typically, `alpha` will be either 1 (default) or 100.

- `"atanh"`:

  Arctanh transformation: \\\tanh^{-1}(y/\alpha)\\. Typically, `alpha`
  will be either 1 (default) or 100.

- `"bcnPower"`:

  Box-Cox with negatives allowed, as described for the `bcnPower`
  function in the car package. It is defined as the Box-Cox
  transformation \\(z^\alpha - 1) / \alpha\\ of the variable \\z = y +
  (y^2+\beta^2)^{1/2}\\. Note that this requires both parameters and
  that `beta > 0`.

- `"scale"`:

  This one is a little different than the others, in that `alpha` and
  `beta` are ignored; instead, they are determined by calling
  `scale(y, ...)`. The user should give as `y` the response variable in
  the model to be fitted to its scaled version.

Note that with the `"power"`, `"boxcox"`, or `"sympower"`
transformations, the argument `beta` specifies a location shift. In the
`"genpower"` transformation, `beta` specifies the base of the logarithm
– however, quirkily, the default of `beta = 0` is taken to be the
natural logarithm. For example, `make.tran(0.5, 10)` sets up the
\\\log\_{10}(y + \frac12)\\ transformation. In the `"bcnPower"`
transformation, `beta` must be specified as a positive value.

For purposes of back-transformation, the `sqrt(y) + sqrt(y+1)`
transformation is treated exactly the same way as `2*sqrt(y)`, because
both are regarded as estimates of \\2\sqrt\mu\\.

## Cases where `make.tran` may not be needed

For standard transformations with no parameters, we usually don't need
to use `make.tran`; just the name of the transformation is all that is
needed. The functions
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md),
[`ref_grid`](https://rvlenth.github.io/emmeans/reference/ref_grid.md),
and related ones automatically detect response transformations that are
recognized by examining the model formula. These are `log`, `log2`,
`log10`, `log1p`, `sqrt`, `logit`, `probit`, `cauchit`, `cloglog`; as
well as (for a response variable `y`) `asin(sqrt(y))`, `asinh(sqrt(y))`,
`atanh(y)`, and `sqrt(y) + sqrt(y+1)`. In addition, any constant
multiple of these (e.g., `2*sqrt(y)`) is auto-detected and appropriately
scaled (see also the `tran.mult` argument in
[`update.emmGrid`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md)).

A few additional transformations may be specified as character strings
and are auto-detected: `"identity"`, `"1/mu^2"`, `"inverse"`,
`"reciprocal"`, `"log10"`, `"log2"`, `"asin.sqrt"`, `"asinh.sqrt"`, and
`"atanh"`.

## Compound transformations

A transformation that is a function of another function can be created
by specifying `inner` for the other function. For example, the
transformation \\1/\sqrt{y}\\ can be created either by
`make.tran("inverse", inner = "sqrt")` or by `make.tran("power", -0.5)`.
In principle, transformations can be compounded to any depth. Also, if
`type` is `"scale"`, `y` is replaced by `inner$linkfun(y)`, because that
will be the variable that is scaled.

## Examples

``` r
# Fit a model using an oddball transformation:
bctran <- make.tran("boxcox", 0.368)
warp.bc <- with(bctran, 
    lm(linkfun(breaks) ~ wool * tension, data = warpbreaks))
# Obtain back-transformed LS means:    
emmeans(warp.bc, ~ tension | wool, type = "response")
#> wool = A:
#>  tension response   SE df lower.CL upper.CL
#>  L           42.4 4.48 48     34.0     52.0
#>  M           23.1 3.05 48     17.5     29.8
#>  H           23.3 3.07 48     17.7     30.0
#> 
#> wool = B:
#>  tension response   SE df lower.CL upper.CL
#>  L           27.2 3.38 48     20.9     34.6
#>  M           27.9 3.44 48     21.5     35.3
#>  H           18.4 2.65 48     13.6     24.3
#> 
#> Confidence level used: 0.95 
#> Intervals are back-transformed from the Box-Cox (lambda = 0.368) scale 

### Using a scaled response...
# Case where it is auto-detected:
mod <- lm(scale(yield[, 1]) ~ Variety, data = MOats)
emmeans(mod, "Variety", type = "response")
#>  Variety     response   SE df lower.CL upper.CL
#>  Golden Rain     80.0 7.96 15     63.0     97.0
#>  Marvellous      86.7 7.96 15     69.7    103.6
#>  Victory         71.5 7.96 15     54.5     88.5
#> 
#> Confidence level used: 0.95 
#> Intervals are back-transformed from the scale(79.4, 19.4) scale 

# Case where scaling is not auto-detected -- and what to do about it:
copt <- options(contrasts = c("contr.sum", "contr.poly"))
mod.aov <- aov(scale(yield[, 1]) ~ Variety + Error(Block), data = MOats)
emm.aov <- suppressWarnings(emmeans(mod.aov, "Variety", type = "response"))
#> NOTE: Unable to recover scale() parameters. See '? make.tran'

# Scaling was not retrieved, but we can do:
emm.aov <- update(emm.aov, tran = make.tran("scale", y = MOats$yield[, 1]))
emmeans(emm.aov, "Variety", type = "response")
#>  Variety     response   SE  df lower.CL upper.CL
#>  Golden Rain     80.0 7.96 9.8     62.2     97.8
#>  Marvellous      86.7 7.96 9.8     68.9    104.5
#>  Victory         71.5 7.96 9.8     53.7     89.3
#> 
#> Warning: EMMs are biased unless design is perfectly balanced 
#> Confidence level used: 0.95 
#> Intervals are back-transformed from the scale(79.4, 19.4) scale 

### Compound transformations
# The following amount to the same thing:
t1 <- make.tran("inverse", inner = "sqrt")
t2 <- make.tran("power", -0.5)

options(copt)


if (FALSE) { # \dontrun{
### An existing model 'mod' was fitted with a y^(2/3) transformation...
  ptran = make.tran("power", 2/3)
  emmeans(mod, "treatment", tran = ptran)
} # }

pigs.lm <- lm(inverse(conc) ~ source + factor(percent), data = pigs)
emmeans(pigs.lm, "source", type = "response")
#>  source response    SE df lower.CL upper.CL
#>  fish       29.7 0.816 23     28.1     31.5
#>  soy        39.0 1.440 23     36.2     42.2
#>  skim       43.8 1.900 23     40.1     48.1
#> 
#> Results are averaged over the levels of: percent 
#> Confidence level used: 0.95 
#> Intervals are back-transformed from the inverse scale 
```
