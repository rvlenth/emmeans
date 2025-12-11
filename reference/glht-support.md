# Support for `multcomp::glht`

These functions and methods provide an interface between emmeans and the
[`multcomp::glht`](https://rdrr.io/pkg/multcomp/man/glht.html) function
for simultaneous inference provided by the multcomp package.

## Usage

``` r
emm(...)

as.glht(object, ...)

# S3 method for class 'emmGrid'
as.glht(object, ...)
```

## Arguments

- ...:

  In `emm`, the `specs`, `by`, and `contr` arguments you would normally
  supply to
  [`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md).
  Only `specs` is required. Otherwise, arguments are passed to other
  methods. You may also include a `which` argument; see Details.

- object:

  An object of class `emmGrid` or `emm_list`

## Value

`emm` returns an object of an intermediate class for which there is a
[`multcomp::glht`](https://rdrr.io/pkg/multcomp/man/glht.html) method.

`as.glht` returns an object of class `glht` or `glht_list` according to
whether `object` is of class `emmGrid` or `emm_list`. See Details below
for more on `glht_list`s.

## Note

The multivariate-\\t\\ routines used by `glht` require that all
estimates in the family have the same integer degrees of freedom. In
cases where that is not true, a message is displayed that shows what df
is used. The user may override this via the `df` argument.

## Details for `emm`

`emm` is meant to be called only *from* `"glht"` as its second
(`linfct`) argument. It works similarly to
[`multcomp::mcp`](https://rdrr.io/pkg/multcomp/man/glht.html), except
with `specs` (and optionally `by` and `contr` arguments) provided as in
a call to
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md).

If the specifications in `...` would result in a list (i.e., an
`emm_list` object), then by default, only the last element of that list
is passed to `glht`. However, if `...` contains a `which` argument
consisting of integer values, the list elements with those indexes are
selected and combined and passed on to `glht`. No checking is done on
whether the indexes are valid, and the keyword `which` must be
spelled-out.

## Details for `as.glht`

When no `by` variable is in force, we obtain a `glht` object; otherwise
it is a `glht_list`. The latter is defined in emmeans, not multcomp, and
is simply a `list` of `glht` objects. Appropriate convenience methods
`coef`, `confint`, `plot`, `summary`, and `vcov` are provided, which
simply apply the corresponding `glht` methods to each member.

## Examples

``` r
if(require(multcomp, quietly = TRUE)) 
    emm_example("glht-multcomp") 
#> 
#> Attaching package: ‘TH.data’
#> The following object is masked from ‘package:MASS’:
#> 
#>     geyser
#> 
#> --- Running code from 'system.file("extexamples", "glht-multcomp.R", package = "emmeans")'
#> 
#> > warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
#> 
#> > summary(glht(warp.lm, emm(pairwise ~ tension | wool)))
#> $`wool = A`
#> 
#>   Simultaneous Tests for General Linear Hypotheses
#> 
#> Fit: lm(formula = breaks ~ wool * tension, data = warpbreaks)
#> 
#> Linear Hypotheses:
#>            Estimate Std. Error t value Pr(>|t|)    
#> L - M == 0  20.5556     5.1573   3.986   <0.001 ***
#> L - H == 0  20.0000     5.1573   3.878   <0.001 ***
#> M - H == 0  -0.5556     5.1573  -0.108    0.994    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> (Adjusted p values reported -- single-step method)
#> 
#> 
#> $`wool = B`
#> 
#>   Simultaneous Tests for General Linear Hypotheses
#> 
#> Fit: lm(formula = breaks ~ wool * tension, data = warpbreaks)
#> 
#> Linear Hypotheses:
#>            Estimate Std. Error t value Pr(>|t|)
#> L - M == 0  -0.5556     5.1573  -0.108    0.994
#> L - H == 0   9.4444     5.1573   1.831    0.170
#> M - H == 0  10.0000     5.1573   1.939    0.139
#> (Adjusted p values reported -- single-step method)
#> 
#> 
#> 
#> > summary(glht(warp.lm, emm(pairwise ~ tension | wool, 
#> +     which = 1:2, by = "wool")))
#> $`wool = A`
#> 
#>   Simultaneous Tests for General Linear Hypotheses
#> 
#> Fit: lm(formula = breaks ~ wool * tension, data = warpbreaks)
#> 
#> Linear Hypotheses:
#>               Estimate Std. Error t value Pr(>|t|)    
#> L, . == 0      44.5556     3.6468  12.218  < 0.001 ***
#> M, . == 0      24.0000     3.6468   6.581  < 0.001 ***
#> H, . == 0      24.5556     3.6468   6.734  < 0.001 ***
#> ., L - M == 0  20.5556     5.1573   3.986  0.00113 ** 
#> ., L - H == 0  20.0000     5.1573   3.878  0.00156 ** 
#> ., M - H == 0  -0.5556     5.1573  -0.108  0.99949    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> (Adjusted p values reported -- single-step method)
#> 
#> 
#> $`wool = B`
#> 
#>   Simultaneous Tests for General Linear Hypotheses
#> 
#> Fit: lm(formula = breaks ~ wool * tension, data = warpbreaks)
#> 
#> Linear Hypotheses:
#>               Estimate Std. Error t value Pr(>|t|)    
#> L, . == 0      28.2222     3.6468   7.739   <1e-04 ***
#> M, . == 0      28.7778     3.6468   7.891   <1e-04 ***
#> H, . == 0      18.7778     3.6468   5.149   <1e-04 ***
#> ., L - M == 0  -0.5556     5.1573  -0.108    0.999    
#> ., L - H == 0   9.4444     5.1573   1.831    0.253    
#> ., M - H == 0  10.0000     5.1573   1.939    0.210    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> (Adjusted p values reported -- single-step method)
#> 
#> 
#> 
#> > warp.emm <- emmeans(warp.lm, ~tension | wool)
#> 
#> > summary(as.glht(pairs(warp.emm)))
#> $`wool = A`
#> 
#>   Simultaneous Tests for General Linear Hypotheses
#> 
#> Linear Hypotheses:
#>            Estimate Std. Error t value Pr(>|t|)    
#> L - M == 0  20.5556     5.1573   3.986 0.000625 ***
#> L - H == 0  20.0000     5.1573   3.878 0.000867 ***
#> M - H == 0  -0.5556     5.1573  -0.108 0.993622    
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> (Adjusted p values reported -- single-step method)
#> 
#> 
#> $`wool = B`
#> 
#>   Simultaneous Tests for General Linear Hypotheses
#> 
#> Linear Hypotheses:
#>            Estimate Std. Error t value Pr(>|t|)
#> L - M == 0  -0.5556     5.1573  -0.108    0.994
#> L - H == 0   9.4444     5.1573   1.831    0.170
#> M - H == 0  10.0000     5.1573   1.939    0.139
#> (Adjusted p values reported -- single-step method)
#> 
#> 
#> 
#> > summary(as.glht(pairs(warp.emm), by = NULL))
#> 
#>   Simultaneous Tests for General Linear Hypotheses
#> 
#> Linear Hypotheses:
#>               Estimate Std. Error t value Pr(>|t|)   
#> L - M, A == 0  20.5556     5.1573   3.986  0.00132 **
#> L - H, A == 0  20.0000     5.1573   3.878  0.00181 **
#> M - H, A == 0  -0.5556     5.1573  -0.108  0.99996   
#> L - M, B == 0  -0.5556     5.1573  -0.108  0.99996   
#> L - H, B == 0   9.4444     5.1573   1.831  0.30800   
#> M - H, B == 0  10.0000     5.1573   1.939  0.25533   
#> ---
#> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#> (Adjusted p values reported -- single-step method)
#> 
#> 
    # Use emm_example("glht-multcomp", list = TRUE) # to see just the code
    
```
