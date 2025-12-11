# Miscellaneous methods for `emmGrid` objects

Miscellaneous methods for `emmGrid` objects

## Usage

``` r
# S3 method for class 'emmGrid'
str(object, ...)

# S3 method for class 'emmGrid'
print(x, ..., export = FALSE)

# S3 method for class 'emmGrid'
vcov(object, ..., sep = get_emm_option("sep"))

linfct(object, ...)

# Default S3 method
linfct(object, ...)
```

## Arguments

- object:

  An `emmGrid` object

- ...:

  (required but not used)

- x:

  An `emmGrid` object

- export:

  Logical value. If `FALSE`, the object is printed. If `TRUE`, a list is
  invisibly returned, which contains character elements named `summary`
  and `annotations` that may be saved or displayed as the user sees fit.
  `summary` is a character matrix (or list of such matrices, if a `by`
  variable is in effect). `annotations` is a character vector of the
  annotations that would have been printed below the summary or
  summaries.

- sep:

  separator for pasting levels in creating row and column names for
  [`vcov()`](https://rdrr.io/r/stats/vcov.html) results

## Value

The `vcov` method returns a symmetric matrix of variances and
covariances for `predict.emmGrid(object, type = "lp")`

The `linfct` function and method returns the `linfct` slot of `object`.

## Examples

``` r
warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
warp.emm <- emmeans(warp.lm, ~ tension | wool)
vcov(warp.emm) |> zapsmall()
#>          L A      M A      H A      L B      M B      H B
#> L A 13.29887  0.00000  0.00000  0.00000  0.00000  0.00000
#> M A  0.00000 13.29887  0.00000  0.00000  0.00000  0.00000
#> H A  0.00000  0.00000 13.29887  0.00000  0.00000  0.00000
#> L B  0.00000  0.00000  0.00000 13.29887  0.00000  0.00000
#> M B  0.00000  0.00000  0.00000  0.00000 13.29887  0.00000
#> H B  0.00000  0.00000  0.00000  0.00000  0.00000 13.29887

vcov(pairs(warp.emm), sep = "|") |> zapsmall()
#>           L - M|A  L - H|A   M - H|A   L - M|B  L - H|B   M - H|B
#> L - M|A  26.59774 13.29887 -13.29887   0.00000  0.00000   0.00000
#> L - H|A  13.29887 26.59774  13.29887   0.00000  0.00000   0.00000
#> M - H|A -13.29887 13.29887  26.59774   0.00000  0.00000   0.00000
#> L - M|B   0.00000  0.00000   0.00000  26.59774 13.29887 -13.29887
#> L - H|B   0.00000  0.00000   0.00000  13.29887 26.59774  13.29887
#> M - H|B   0.00000  0.00000   0.00000 -13.29887 13.29887  26.59774
```
