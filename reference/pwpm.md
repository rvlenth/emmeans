# Pairwise P-value matrix (plus other statistics)

This function presents results from `emmeans` and pairwise comparisons
thereof in a compact way. It displays a matrix (or matrices) of
estimates, pairwise differences, and P values. The user may opt to
exclude any of these via arguments `means`, `diffs`, and `pvals`,
respectively. To control the direction of the pairwise differences, use
`reverse`; and to control what appears in the upper and lower
triangle(s), use `flip`. Optional arguments are passed to
`contrast.emmGrid` and/or `summary.emmGrid`, making it possible to
control what estimates and tests are displayed.

## Usage

``` r
pwpm(emm, by, reverse = FALSE, pvals = TRUE, means = TRUE,
  diffs = TRUE, flip = FALSE, digits, ...)
```

## Arguments

- emm:

  An `emmGrid` object

- by:

  Character vector of variable(s) in the grid to condition on. These
  will create different matrices, one for each level or
  level-combination. If missing, `by` is set to `emm@misc$by.vars`. Grid
  factors not in `by` are the *primary* factors: whose levels or level
  combinations are compared pairwise.

- reverse:

  Logical value passed to
  [`pairs.emmGrid`](https://rvlenth.github.io/emmeans/reference/contrast.md).
  Thus, `FALSE` specifies `"pairwise"` comparisons (earlier vs. later),
  and `TRUE` specifies `"revpairwise"` comparisons (later vs. earlier).

- pvals:

  Logical value. If `TRUE`, the pairwise differences of the EMMs are
  included in each matrix according to `flip`.

- means:

  Logical value. If `TRUE`, the estimated marginal means (EMMs) from
  `emm` are included in the matrix diagonal(s).

- diffs:

  Logical value. If `TRUE`, the pairwise differences of the EMMs are
  included in each matrix according to `flip`.

- flip:

  Logical value that determines where P values and differences are
  placed. `FALSE` places the P values in the upper triangle and
  differences in the lower, and `TRUE` does just the opposite.

- digits:

  Integer. Number of digits to display. If missing, an optimal number of
  digits is determined.

- ...:

  Additional arguments passed to
  [`contrast.emmGrid`](https://rvlenth.github.io/emmeans/reference/contrast.md)
  and
  [`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md).
  You should *not* include `method` here, because pairwise comparisons
  are always used.

## Value

A matrix or \`list\` of matrices, one for each \`by\` level.

## Note

If `emm` is the result of a Bayesian analysis, `pwpm` is based on a
frequentist analysis

## See also

A graphical display of essentially the same results is available from
[`pwpp`](https://rvlenth.github.io/emmeans/reference/pwpp.md)

## Examples

``` r
warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
warp.emm <- emmeans(warp.lm, ~ tension | wool)

pwpm(warp.emm)
#> 
#> wool = A
#>        L      M      H
#> L [44.6] 0.0007 0.0009
#> M 20.556 [24.0] 0.9936
#> H 20.000 -0.556 [24.6]
#> 
#> wool = B
#>        L      M      H
#> L [28.2] 0.9936 0.1704
#> M -0.556 [28.8] 0.1389
#> H  9.444 10.000 [18.8]
#> 
#> Row and column labels: tension
#> Upper triangle: P values   adjust = “tukey”
#> Diagonal: [Estimates] (emmean) 
#> Lower triangle: Comparisons (estimate)   earlier vs. later

# use dot options to specify noninferiority tests
pwpm(warp.emm, by = NULL, side = ">", delta = 5, adjust = "none")
#>        L A    M A    H A    L B    M B    H B
#> L A [44.6] <.0001 <.0001 <.0001 <.0001 <.0001
#> M A 20.556 [24.0] 0.1965 0.4404 0.4829 0.0266
#> H A 20.000 -0.556 [24.6] 0.3986 0.4404 0.0210
#> L B 16.333 -4.222 -3.667 [28.2] 0.1965 0.0037
#> M B 15.778 -4.778 -4.222 -0.556 [28.8] 0.0027
#> H B 25.778  5.222  5.778  9.444 10.000 [18.8]
#> 
#> Row and column labels: tension:wool
#> Upper triangle: P values   side = “>”  delta = 5
#> Diagonal: [Estimates] (emmean) 
#> Lower triangle: Comparisons (estimate)   earlier vs. later
```
