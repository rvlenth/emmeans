# Oats data in multivariate form

This is the `Oats` dataset provided in the nlme package, but it is
rearranged as one multivariate observation per plot.

## Usage

``` r
MOats
```

## Format

A data frame with 18 observations and 3 variables

- `Variety`:

  a factor with levels `Golden Rain`, `Marvellous`, `Victory`

- `Block`:

  an ordered factor with levels `VI` \< `V` \< `III` \< `IV` \< `II` \<
  `I`

- `yield`:

  a matrix with 4 columns, giving the yields with nitrogen
  concentrations of 0, .2, .4, and .6.

## Source

The dataset [`Oats`](https://rdrr.io/pkg/nlme/man/Oats.html) in the nlme
package.

## Details

These data arise from a split-plot experiment reported by Yates (1935)
and used as an example in Pinheiro and Bates (2000) and other texts. Six
blocks were divided into three whole plots, randomly assigned to the
three varieties of oats. The whole plots were each divided into 4 split
plots and randomized to the four concentrations of nitrogen.

## References

Pinheiro, J. C. and Bates D. M. (2000) *Mixed-Effects Models in S and
S-PLUS*, Springer, New York. (Appendix A.15)

Yates, F. (1935) Complex experiments, *Journal of the Royal Statistical
Society* Suppl. 2, 181-247

## Examples

``` r
MOats.lm <- lm (yield ~ Block + Variety, data = MOats)
MOats.rg <- ref_grid (MOats.lm, mult.name = "nitro")
emmeans(MOats.rg, ~ nitro | Variety)
#> Variety = Golden Rain:
#>  nitro emmean   SE df lower.CL upper.CL
#>  0       80.0 5.54 10     67.7     92.3
#>  0.2     98.5 6.60 10     83.8    113.2
#>  0.4    114.7 8.70 10     95.3    134.0
#>  0.6    124.8 7.30 10    108.6    141.1
#> 
#> Variety = Marvellous:
#>  nitro emmean   SE df lower.CL upper.CL
#>  0       86.7 5.54 10     74.3     99.0
#>  0.2    108.5 6.60 10     93.8    123.2
#>  0.4    117.2 8.70 10     97.8    136.5
#>  0.6    126.8 7.30 10    110.6    143.1
#> 
#> Variety = Victory:
#>  nitro emmean   SE df lower.CL upper.CL
#>  0       71.5 5.54 10     59.2     83.8
#>  0.2     89.7 6.60 10     75.0    104.4
#>  0.4    110.8 8.70 10     91.5    130.2
#>  0.6    118.5 7.30 10    102.2    134.8
#> 
#> Results are averaged over the levels of: Block 
#> Confidence level used: 0.95 
```
