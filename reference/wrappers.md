# Wrappers for alternative naming of EMMs

These are wrappers for
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md) and
related functions to provide backward compatibility, or for users who
may prefer to use other terminology than “estimated marginal means” –
namely “least-squares means”. These functions also provide the
functionality formerly provided by the lsmeans package, which is now
just a front-end for emmeans.

## Usage

``` r
lsmeans(...)

lstrends(...)

lsmip(...)

lsm(...)

lsmobj(...)

lsm.options(...)

get.lsm.option(x, default = emm_defaults[[x]])
```

## Arguments

- ...:

  Arguments passed to the corresponding `em`*xxxx* function

- x:

  Character name of desired option

- default:

  default value to return if `x` not found

## Value

The result of the call to `em`*xxxx*, suitably modified.

`get.lsm.option` and `lsm.options` remap options from and to
corresponding options in the emmeans options system.

## Details

For each function with `ls`*xxxx* in its name, the same function named
`em`*xxxx* is called. Any estimator names or list items beginning with
“em” are replaced with “ls” before the results are returned

## See also

[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md),
[`emtrends`](https://rvlenth.github.io/emmeans/reference/emtrends.md),
[`emmip`](https://rvlenth.github.io/emmeans/reference/emmip.md),
[`emm`](https://rvlenth.github.io/emmeans/reference/glht-support.md),
[`emmobj`](https://rvlenth.github.io/emmeans/reference/emmobj.md),
[`emm_options`](https://rvlenth.github.io/emmeans/reference/emm_options.md),
[`get_emm_option`](https://rvlenth.github.io/emmeans/reference/emm_options.md)

## Examples

``` r
pigs.lm <- lm(log(conc) ~ source + factor(percent), data = pigs)
lsmeans(pigs.lm, "source")
#>  source lsmean     SE df lower.CL upper.CL
#>  fish     3.39 0.0367 23     3.32     3.47
#>  soy      3.67 0.0374 23     3.59     3.74
#>  skim     3.80 0.0394 23     3.72     3.88
#> 
#> Results are averaged over the levels of: percent 
#> Results are given on the log (not the response) scale. 
#> Confidence level used: 0.95 
```
