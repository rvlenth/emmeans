# Calculate Cohen effect sizes and confidence bounds thereof

Standardized effect sizes are typically calculated using pairwise
differences of estimates, divided by the SD of the population providing
the context for those effects. This function calculates effect sizes
from an `emmGrid` object, and confidence intervals for them, accounting
for uncertainty in both the estimated effects and the population SD.

## Usage

``` r
eff_size(object, sigma, edf, method = "pairwise", ...)
```

## Arguments

- object:

  an
  [`emmGrid`](https://rvlenth.github.io/emmeans/reference/emmGrid-class.md)
  object, typically one defining the EMMs to be contrasted. If instead,
  `class(object) == "emm_list"`, such as is produced by
  `emmeans(model, pairwise ~ treatment)`, a message is displayed; the
  contrasts already therein are used; and `method` is replaced by
  `"identity"`.

- sigma:

  numeric scalar, value of the population SD.

- edf:

  numeric scalar that specifies the equivalent degrees of freedom for
  the `sigma`. This is a way of specifying the uncertainty in `sigma`,
  in that we regard our estimate of `sigma^2` as being proportional to a
  chi-square random variable with `edf` degrees of freedom. (`edf`
  should not be confused with the `df` argument that may be passed via
  `...` to specify the degrees of freedom to use in \\t\\ statistics and
  confidence intervals.)

- method:

  the contrast method to use to define the effects. This is passed to
  [`contrast`](https://rvlenth.github.io/emmeans/reference/contrast.md)
  after the elements of `object` are scaled.

- ...:

  Additional arguments passed to `contrast`

## Value

an
[`emmGrid`](https://rvlenth.github.io/emmeans/reference/emmGrid-class.md)
object containing the effect sizes

## Details

Any `by` variables specified in `object` will remain in force in the
returned effects, unless overridden in the optional arguments.

For models having a single random effect, such as those fitted using
[`lm`](https://rdrr.io/r/stats/lm.html); in that case, the
[`stats::sigma`](https://rdrr.io/r/stats/sigma.html) and
[`stats::df.residual`](https://rdrr.io/r/stats/df.residual.html)
functions may be useful for specifying `sigma` and `edf`. For models
with more than one random effect, `sigma` may be based on some
combination of the random-effect variances.

Specifying `edf` can be rather unintuitive but is also relatively
uncritical; but the smaller the value, the wider the confidence
intervals for effect size. The value of `sqrt(2/edf)` can be interpreted
as the relative accuracy of `sigma`; for example, with `edf = 50`,
\\\sqrt(2/50) = 0.2\\, meaning that `sigma` is accurate to plus or minus
20 percent. Note in an example below, we tried two different `edf`
values as kind of a bracketing/sensitivity-analysis strategy. A value of
`Inf` is allowable, in which case you are assuming that `sigma` is known
exactly. Obviously, this narrows the confidence intervals for the effect
sizes – unrealistically if in fact `sigma` is unknown.

## Note

The effects are always computed on the scale of the *linear-predictor*;
any response transformation or link function is completely ignored. If
you wish to base the effect sizes on the response scale, it is *not*
enough to replace `object` with `regrid(object)`, because this
back-transformation changes the SD required to compute effect sizes.

**Paired data:** Be careful with paired-data situations, where Cohen's d
is typically referenced to the SD of the *paired differences* rather
than the *residual* SD. You may need to enlarge `sigma` by a factor of
`sqrt(2)` to obtain comparable results with other software.

**Disclaimer:** There is substantial disagreement among practitioners on
what is the appropriate `sigma` to use in computing effect sizes; or,
indeed, whether *any* effect-size measure is appropriate for some
situations. The user is completely responsible for specifying
appropriate parameters (or for failing to do so).

Cohen effect sizes do not even exist for generalized linear models or
other models lacking an additive residual error term.

The examples here illustrate a sobering message that effect sizes are
often not nearly as accurate as you may think.

## Computation

This function uses calls to
[`regrid`](https://rvlenth.github.io/emmeans/reference/regrid.md) to put
the estimated marginal means (EMMs) on the log scale. Then an extra
element is added to this grid for the log of `sigma` and its standard
error (where we assume that `sigma` is uncorrelated with the log EMMs).
Then a call to
[`contrast`](https://rvlenth.github.io/emmeans/reference/contrast.md)
subtracts `log{sigma}` from each of the log EMMs, yielding values of
`log(EMM/sigma)`. Finally, the results are re-gridded back to the
original scale and the desired contrasts are computed using `method`. In
the log-scaling part, we actually rescale the absolute values and keep
track of the signs.

## Examples

``` r
fiber.lm <- lm(strength ~ diameter + machine, data = fiber)

emm <- emmeans(fiber.lm, "machine")
eff_size(emm, sigma = sigma(fiber.lm), edf = df.residual(fiber.lm))
#>  contrast effect.size    SE df lower.CL upper.CL
#>  A - B         -0.650 0.650 11   -2.081    0.781
#>  A - C          0.993 0.726 11   -0.604    2.590
#>  B - C          1.643 0.800 11   -0.118    3.405
#> 
#> sigma used for effect sizes: 1.595 
#> Confidence level used: 0.95 

# or equivalently:
eff_size(pairs(emm), sigma(fiber.lm), df.residual(fiber.lm), method = "identity")
#>  contrast effect.size    SE df lower.CL upper.CL
#>  (A - B)       -0.650 0.650 11   -2.081    0.781
#>  (A - C)        0.993 0.726 11   -0.604    2.590
#>  (B - C)        1.643 0.800 11   -0.118    3.405
#> 
#> sigma used for effect sizes: 1.595 
#> Confidence level used: 0.95 


### Mixed model example:
if (require(nlme)) withAutoprint({
  Oats.lme <- lme(yield ~ Variety + factor(nitro), 
                  random = ~ 1 | Block / Variety,
                  data = Oats)
  # Combine variance estimates
  VarCorr(Oats.lme)
  (totSD <- sqrt(214.4724 + 109.6931 + 162.5590))
  # I figure edf is somewhere between 5 (Blocks df) and 51 (Resid df)
  emmV <- emmeans(Oats.lme, ~ Variety)
  eff_size(emmV, sigma = totSD, edf = 5)
  eff_size(emmV, sigma = totSD, edf = 51)
}, spaced = TRUE)
#> Loading required package: nlme
#> 
#> > Oats.lme <- lme(yield ~ Variety + factor(nitro), random = ~1 | Block/Variety, 
#> +     data = Oats)
#> 
#> > VarCorr(Oats.lme)
#>             Variance     StdDev  
#> Block =     pdLogChol(1)         
#> (Intercept) 214.4722     14.64487
#> Variety =   pdLogChol(1)         
#> (Intercept) 109.6928     10.47343
#> Residual    162.5591     12.74987
#> 
#> > (totSD <- sqrt(214.4724 + 109.6931 + 162.559))
#> [1] 22.06183
#> 
#> > emmV <- emmeans(Oats.lme, ~Variety)
#> 
#> > eff_size(emmV, sigma = totSD, edf = 5)
#>  contrast                 effect.size    SE df lower.CL upper.CL
#>  Golden Rain - Marvellous      -0.240 0.330  5   -1.087    0.608
#>  Golden Rain - Victory          0.312 0.336  5   -0.551    1.174
#>  Marvellous - Victory           0.551 0.365  5   -0.387    1.490
#> 
#> Results are averaged over the levels of: nitro 
#> sigma used for effect sizes: 22.06 
#> Degrees-of-freedom method: inherited from containment when re-gridding 
#> Confidence level used: 0.95 
#> 
#> > eff_size(emmV, sigma = totSD, edf = 51)
#>  contrast                 effect.size    SE df lower.CL upper.CL
#>  Golden Rain - Marvellous      -0.240 0.322  5   -1.067    0.587
#>  Golden Rain - Victory          0.312 0.322  5   -0.517    1.140
#>  Marvellous - Victory           0.551 0.325  5   -0.285    1.388
#> 
#> Results are averaged over the levels of: nitro 
#> sigma used for effect sizes: 22.06 
#> Degrees-of-freedom method: inherited from containment when re-gridding 
#> Confidence level used: 0.95 

```
