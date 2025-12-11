# Combine or subset `emmGrid` objects

These functions provide methods for
[`rbind`](https://rdrr.io/r/base/cbind.html) and
[`[`](https://rdrr.io/r/base/Extract.html) that may be used to combine
`emmGrid` objects together, or to extract a subset of cases. The primary
reason for doing this would be to obtain multiplicity-adjusted results
for smaller or larger families of tests or confidence intervals.

## Usage

``` r
# S3 method for class 'emmGrid'
rbind(..., deparse.level = 1, adjust = "bonferroni")

# S3 method for class 'emmGrid'
e1 + e2

# S3 method for class 'emmGrid'
x[i, adjust, drop.levels = TRUE, ...]

# S3 method for class 'emmGrid'
head(x, n = 6, ...)

# S3 method for class 'emmGrid'
tail(x, n = 6, ...)

# S3 method for class 'emmGrid'
subset(x, subset, ...)

# S3 method for class 'emm_list'
rbind(..., which, adjust = "bonferroni")

# S3 method for class 'summary_emm'
rbind(..., which)

force_regular(object)
```

## Arguments

- ...:

  In `rbind`, object(s) of class `emmGrid` or `summary_emm`. In others,
  additional arguments passed to other methods

- deparse.level:

  (required but not used)

- adjust:

  Character value passed to
  [`update.emmGrid`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md)

- e1, e2, x, object:

  Objects of class `emmGrid`

- i:

  Integer vector of indexes

- drop.levels:

  Logical value. If `TRUE`, the `"levels"` slot in the returned object
  is updated to hold only the predictor levels that actually occur

- n:

  integer number of entries to include (or exclude if negative)

- subset:

  logical expression indicating which rows of the grid to keep

- which:

  Integer vector of subset of elements to use; if missing, we use all
  elements.

## Value

A revised object of class `emmGrid`

The result of `e1 + e2` is the same as `rbind(e1, e2)`

The `rbind` method for `emm_list` objects simply combines the `emmGrid`
objects comprising the first element of `...`. Note that the returned
object is not yet summarized, so any `adjust` parameters apply to the
combined `emmGrid`.

The `rbind` method for `summary_emm` objects (or a list thereof) returns
a single `summary_emm` object. This combined object *preserves* any
adjusted P values or confidence limits in the original summaries, since
those quantities have already been computed.

`force_regular` adds extra (invisible) rows to an `emmGrid` object to
make it a regular grid (all combinations of factors). This regular
structure is needed by `emmeans`. An object can become irregular by, for
example, subsetting rows, or by obtaining contrasts of a nested
structure.

## Note

`rbind` throws an error if there are incompatibilities in the objects'
coefficients, covariance structures, etc. But they are allowed to have
different factors; a missing level `'.'` is added to factors as needed.

These functions generally reset `by.vars` to `NULL`; so if you want to
keep any “by” variables, you should follow-up with
[`update.emmGrid`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md).

## Examples

``` r
warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
warp.rg <- ref_grid(warp.lm)

# Do all pairwise comparisons within rows or within columns, 
# all considered as one faily of tests:
w.t <- pairs(emmeans(warp.rg, ~ wool | tension))
t.w <- pairs(emmeans(warp.rg, ~ tension | wool))
rbind(w.t, t.w, adjust = "mvt")
#>  tension wool contrast estimate   SE df t.ratio p.value
#>  L       .    A - B      16.333 5.16 48   3.167  0.0206
#>  M       .    A - B      -4.778 5.16 48  -0.926  0.9119
#>  H       .    A - B       5.778 5.16 48   1.120  0.8258
#>  .       A    L - M      20.556 5.16 48   3.986  0.0019
#>  .       A    L - H      20.000 5.16 48   3.878  0.0027
#>  .       A    M - H      -0.556 5.16 48  -0.108  1.0000
#>  .       B    L - M      -0.556 5.16 48  -0.108  1.0000
#>  .       B    L - H       9.444 5.16 48   1.831  0.3795
#>  .       B    M - H      10.000 5.16 48   1.939  0.3192
#> 
#> P value adjustment: mvt method for 9 tests 
update(w.t + t.w, adjust = "fdr")  ## same as above except for adjustment
#>  tension wool contrast estimate   SE df t.ratio p.value
#>  L       .    A - B      16.333 5.16 48   3.167  0.0080
#>  M       .    A - B      -4.778 5.16 48  -0.926  0.4614
#>  H       .    A - B       5.778 5.16 48   1.120  0.4022
#>  .       A    L - M      20.556 5.16 48   3.986  0.0014
#>  .       A    L - H      20.000 5.16 48   3.878  0.0014
#>  .       A    M - H      -0.556 5.16 48  -0.108  0.9147
#>  .       B    L - M      -0.556 5.16 48  -0.108  0.9147
#>  .       B    L - H       9.444 5.16 48   1.831  0.1319
#>  .       B    M - H      10.000 5.16 48   1.939  0.1314
#> 
#> P value adjustment: fdr method for 9 tests 

# Show only 3 of the 6 cases
summary(warp.rg[c(2, 4, 5)])
#>  wool tension prediction   SE df
#>  B    L             28.2 3.65 48
#>  B    M             28.8 3.65 48
#>  A    H             24.6 3.65 48
#> 

# After-the-fact 'at' specification
subset(warp.rg, wool == "A")  ## or warp.rg |> subset(wool == "A")
#>  wool tension prediction   SE df
#>  A    L             44.6 3.65 48
#>  A    M             24.0 3.65 48
#>  A    H             24.6 3.65 48
#> 


### Working with 'emm_list' objects
mod <- lm(conc ~ source + factor(percent), data = pigs)
all <- emmeans(mod, list(src = pairwise ~ source, pct = consec ~ percent))
rbind(all, which = c(2, 4), adjust = "mvt")
#>  src.contrast pct.contrast          estimate   SE df t.ratio p.value
#>  fish - soy   .                        -9.47 2.33 23  -4.059  0.0028
#>  fish - skim  .                       -15.58 2.39 23  -6.526 <0.0001
#>  soy - skim   .                        -6.11 2.34 23  -2.613  0.0783
#>  .            percent12 - percent9      6.36 2.47 23   2.570  0.0852
#>  .            percent15 - percent12     1.96 2.57 23   0.763  0.9397
#>  .            percent18 - percent15     3.31 3.04 23   1.088  0.7944
#> 
#> Results are averaged over some or all of the levels of: percent, source 
#> P value adjustment: mvt method for 6 tests 

### Irregular object
tmp <- warp.rg[-1]
## emmeans(tmp, "tension")   # will fail because tmp is irregular
emmeans(force_regular(tmp), "tension")   # will show some results
#> Warning: emmeans() results may be corrupted by removal of a nesting structure
#> NOTE: Results may be misleading due to involvement in interactions
#>  tension emmean   SE df lower.CL upper.CL
#>  L       nonEst   NA NA       NA       NA
#>  M         26.4 2.58 48     21.2     31.6
#>  H         21.7 2.58 48     16.5     26.9
#> 
#> Results are averaged over the levels of: wool 
#> Confidence level used: 0.95 
```
