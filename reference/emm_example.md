# Run or list additional examples

This function exists so as to provide cleaner-looking examples in help
files when it must be run conditionally on another package. Typically we
want to run the code (`run = TRUE` is the default), or otherwise just
list it on the console (`list = TRUE`).

## Usage

``` r
emm_example(name, run = !list, list = FALSE, ...)
```

## Arguments

- name:

  Character name of file to run. We look for a file with this name (with
  `".R"` appended) in the system files provided with emmeans.

- run:

  Logical choosing whether or not to run the example code

- list:

  Logical choosing whether or not to list the example code

- ...:

  Used only by the developer

## Examples

``` r
# List an example
emm_example("qdrg-biglm", list = TRUE)
#> # Post hoc analysis of a "biglm" object -- not supported by emmeans
#> bigmod <- biglm(log(conc) ~ source + factor(percent), data = pigs)
#> 
#> rg1 <- qdrg(log(conc) ~ source + factor(percent), data = pigs, 
#>     coef = coef(bigmod), vcov = vcov(bigmod), df = bigmod$df.residual)
#> 
#> emmeans(rg1, "source", type = "response")
#> 
#> ## But in this particular case, we could have done it the easy way:
#> ##     rg1 <- qdrg(object = bigmod)
#> 

# Run an example
if (require(biglm))
    emm_example("qdrg-biglm")
#> Loading required package: biglm
#> Loading required package: DBI
#> 
#> --- Running code from 'system.file("extexamples", "qdrg-biglm.R", package = "emmeans")'
#> 
#> > bigmod <- biglm(log(conc) ~ source + factor(percent), 
#> +     data = pigs)
#> 
#> > rg1 <- qdrg(log(conc) ~ source + factor(percent), 
#> +     data = pigs, coef = coef(bigmod), vcov = vcov(bigmod), df = bigmod$df.residual)
#> 
#> > emmeans(rg1, "source", type = "response")
#>  source response   SE  df asymp.LCL asymp.UCL
#>  fish       29.8 1.09 Inf      27.7      32.0
#>  soy        39.1 1.47 Inf      36.4      42.1
#>  skim       44.6 1.75 Inf      41.2      48.1
#> 
#> Results are averaged over the levels of: percent 
#> Confidence level used: 0.95 
#> Intervals are back-transformed from the log scale 
#> 
```
