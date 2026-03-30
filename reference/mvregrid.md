# Multivariate regridding

This function is similar to
[`regrid`](https://rvlenth.github.io/emmeans/reference/regrid.md) except
it performs a multivariate transformation. This is useful, for instance,
in multivariate models that have a compositional response.

## Usage

``` r
mvregrid(object, transform = "response", mult.name, newname = mult.name,
  newlevs, fcn, ...)
```

## Arguments

- object:

  An `emmGrid` object

- transform:

  The transformation to use in re-gridding. If `"response"`, we apply
  the inverse of the multivariate transformation in `object@misc$tran`;
  otherwise, we re-grid as if `transform` had been applied to the
  multivariate response. (Note that this will entail first re-gridding
  to the response scale if necessary.)

- mult.name:

  The name of the multivariate factor to be transformed. If missing, we
  use `object@roles$multresp`, and throw an error message if it is
  `NULL` or ambiguous; in that case, the user must repeat the call with
  `mult.name` specified.

- newname:

  the name to be given to the newly transformed variable

- newlevs:

  levels of the newly created factor (must conform to the number of
  columns created by `fcn`). If missing, we use the column names of the
  newly created variable.

- fcn:

  The multivariate function to apply. If character, we look for it in
  the namespace of the compositions package.

- ...:

  Additional arguments passed to `fcn`

## Value

A new `emmGrid` object with the newly created factor as its last factor

## Details

If a multivariate response transformation was used in fitting the model,
its name is auto-detected, and in that case we need not specify `fcn` as
long as its inverse can be found in the namespace of the compositions
package. (That package need not be installed unless `fcn` is a character
value.) For some such models, auto-detection process throws a warning
message, especially if `cbind` is also present in the model formula.

Currently, no bias-adjustment option is available.

## Examples

``` r
if(requireNamespace("compositions"))
    emm_example("mvregrid")
#> Loading required namespace: compositions
#> 
#> --- Running code from 'system.file("extexamples", "mvregrid.R", package = "emmeans")'
#> 
#> > require("compositions")
#> Loading required package: compositions
#> Welcome to compositions, a package for compositional data analysis.
#> Find an intro with "? compositions"
#> 
#> Attaching package: ‘compositions’
#> The following objects are masked from ‘package:stats’:
#> 
#>     anova, cor, cov, dist, var
#> The following object is masked from ‘package:graphics’:
#> 
#>     segments
#> The following objects are masked from ‘package:base’:
#> 
#>     %*%, norm, scale, scale.default
#> 
#> > data(AnimalVegetation)
#> 
#> > AV <- as.data.frame(AnimalVegetation)
#> 
#> > AV$region <- factor(ifelse(AV$regA == 1, "A", "B"))
#> 
#> > AVmod <- lm(ilr(cbind(disc, spick, din, spin)) ~ region, 
#> +     data = AV)
#> 
#> > AVRG <- suppressWarnings(ref_grid(AVmod))
#> 
#> > AVRG
#> 'emmGrid' object with variables:
#>     region = A, B
#>     rep.meas = multivariate response levels: 1, 2, 3
#> Transformation: “ilr” 
#> 
#> > lvls <- c("disc", "spick", "din", "spin")
#> 
#> > confint(mvregrid(AVRG, newname = "comp", newlevs = lvls), 
#> +     by = "region")
#> region = A:
#>  comp  prediction     SE df lower.CL upper.CL
#>  disc       0.322 0.0217 98    0.279    0.365
#>  spick      0.260 0.0134 98    0.233    0.286
#>  din        0.183 0.0152 98    0.153    0.213
#>  spin       0.235 0.0122 98    0.211    0.260
#> 
#> region = B:
#>  comp  prediction     SE df lower.CL upper.CL
#>  disc       0.349 0.0228 98    0.304    0.394
#>  spick      0.229 0.0123 98    0.205    0.254
#>  din        0.184 0.0155 98    0.153    0.214
#>  spin       0.238 0.0122 98    0.214    0.262
#> 
#> Confidence level used: 0.95 
#> 
#> > mvregrid(AVRG, transform = "clr", newname = "comp", 
#> +     newlevs = lvls)
#>  region comp  prediction     SE df
#>  A      disc      0.2722 0.0773 98
#>  B      disc      0.3617 0.0773 98
#>  A      spick     0.0579 0.0527 98
#>  B      spick    -0.0590 0.0527 98
#>  A      din      -0.2902 0.0731 98
#>  B      din      -0.2813 0.0731 98
#>  A      spin     -0.0398 0.0474 98
#>  B      spin     -0.0214 0.0474 98
#> 
#> Results are given on the clr (not the response) scale. 
#> 
#> > arr.lm <- lm(as.matrix(USArrests) ~ state.region)
#> 
#> > arr.emm <- emmeans(arr.lm, ~crime | state.region, 
#> +     mult.name = "crime")
#> 
#> > PC <- princomp(USArrests, cor = TRUE)
#> 
#> > std.and.rot <- function(x, pc, ...) {
#> +     z <- sweep(sweep(x, 2, pc$center, "-"), 2, pc$scale, "/")
#> +     z %*% pc$loadings
#> + }
#> 
#> > mvregrid(arr.emm, newname = "prin.comp", fcn = std.and.rot, 
#> +     pc = PC)
#> state.region = Northeast:
#>  prin.comp  emmean    SE df lower.CL upper.CL
#>  Comp.1    -1.0350 0.477 46 -1.99507  -0.0749
#>  Comp.2    -0.5710 0.274 46 -1.12204  -0.0199
#>  Comp.3     0.4020 0.159 46  0.08119   0.7229
#>  Comp.4    -0.0924 0.140 46 -0.37370   0.1890
#> 
#> state.region = South:
#>  prin.comp  emmean    SE df lower.CL upper.CL
#>  Comp.1     0.7125 0.358 46 -0.00753   1.4326
#>  Comp.2     0.8652 0.205 46  0.45186   1.2785
#>  Comp.3     0.3153 0.120 46  0.07463   0.5559
#>  Comp.4     0.0886 0.105 46 -0.12241   0.2996
#> 
#> state.region = North Central:
#>  prin.comp  emmean    SE df lower.CL upper.CL
#>  Comp.1    -0.8013 0.413 46 -1.63276   0.0301
#>  Comp.2    -0.1986 0.237 46 -0.67589   0.2786
#>  Comp.3    -0.1127 0.138 46 -0.39056   0.1652
#>  Comp.4     0.1027 0.121 46 -0.14094   0.3464
#> 
#> state.region = West:
#>  prin.comp  emmean    SE df lower.CL upper.CL
#>  Comp.1     0.5793 0.397 46 -0.21955   1.3781
#>  Comp.2    -0.4862 0.228 46 -0.94473  -0.0277
#>  Comp.3    -0.5624 0.133 46 -0.82933  -0.2954
#>  Comp.4    -0.1399 0.116 46 -0.37401   0.0942
#> 
#> Confidence level used: 0.95 
#> 
    # Use emm_example("mvregrid", list = TRUE) # to see just the code
```
