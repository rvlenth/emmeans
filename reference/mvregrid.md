# Multivariate regridding

This function is similar to
[`regrid`](https://rvlenth.github.io/emmeans/reference/regrid.md) except
it performs a multivariate transformation. This is useful, for instance,
in multivariate models that have a compositional response.

## Usage

``` r
mvregrid(object, newname = "component", newlevs = seq_len(ncol(newy)),
  mult.name = names(levels)[length(levels)], fcn = paste0(tran, "Inv"),
  ...)
```

## Arguments

- object:

  An `emmGrid` object

- newname:

  The name to give to the newly created multivariate factor

- newlevs:

  Character levels of the newly created factor (must conform to the
  number of columns created by `fcn`)

- mult.name:

  The name of the multivariate factor to be transformed. By default, we
  use the last factor

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
#> > data(AnimalVegetation, package = "compositions")
#> 
#> > AV <- as.data.frame(AnimalVegetation)
#> 
#> > AVmod <- lm(compositions::ilr(cbind(disc, spick, din, 
#> +     spin)) ~ regA, data = AV)
#> 
#> > AVRG <- suppressWarnings(ref_grid(AVmod))
#> 
#> > AVRG
#> 'emmGrid' object with variables:
#>     regA = 0, 1
#>     rep.meas = multivariate response levels: 1, 2, 3
#> Transformation: “::.compositions.ilr” 
#> 
#> > confint(mvregrid(AVRG, newname = "comp", newlevs = c("disc", 
#> +     "spick", "din", "spin"), fcn = "ilrInv"), by = "regA")
#> regA = 0:
#>  comp  prediction     SE df lower.CL upper.CL
#>  disc       0.349 0.0228 98    0.304    0.394
#>  spick      0.229 0.0123 98    0.205    0.254
#>  din        0.184 0.0155 98    0.153    0.214
#>  spin       0.238 0.0122 98    0.214    0.262
#> 
#> regA = 1:
#>  comp  prediction     SE df lower.CL upper.CL
#>  disc       0.322 0.0217 98    0.279    0.365
#>  spick      0.260 0.0134 98    0.233    0.286
#>  din        0.183 0.0152 98    0.153    0.213
#>  spin       0.235 0.0122 98    0.211    0.260
#> 
#> Confidence level used: 0.95 
#> 
    # Use emm_example("mvregrid", list = TRUE) # to see just the code
```
