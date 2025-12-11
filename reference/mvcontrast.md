# Multivariate contrasts

This function displays tests of multivariate comparisons or contrasts.
The contrasts are constructed at each level of the variable in
`mult.name`, and then we do a multivariate test that the vector of
estimates is equal to `null` (zero by default). The *F* statistic and
degrees of freedom are determined via the Hotelling distribution. that
is, if there are \\m\\ error degrees of freedom and multivariate
dimensionality \\d\\, then the resulting \\F\\ statistic has degrees of
freedom \\(d, m - d + 1)\\ as shown in Hotelling (1931).

## Usage

``` r
mvcontrast(object, method = "eff", mult.name = object@roles$multresp,
  null = 0, by = object@misc$by.vars, adjust = c("sidak",
  p.adjust.methods), show.ests = FALSE, ...)
```

## Arguments

- object:

  An object of class `emmGrid`

- method:

  A contrast method, per
  [`contrast.emmGrid`](https://rvlenth.github.io/emmeans/reference/contrast.md)

- mult.name:

  Character vector of names of the factors whose levels define the
  multivariate means to contrast. If the model itself has a multivariate
  response, that is what is used. Otherwise, `mult.name` *must* be
  specified.

- null:

  Scalar or conformable vector of null-hypothesis values to test against

- by:

  Any `by` variable(s). These should not include the primary variables
  to be contrasted. For convenience, the `by` variable is nulled-out if
  it would result in no primary factors being contrasted.

- adjust:

  Character value of a multiplicity adjustment method (`"none"` for no
  adjustment). The available adjustment methods are more limited that in
  `contrast`, and any default adjustment returned via `method` is
  ignored.

- show.ests:

  Logical flag determining whether the multivariate means are displayed

- ...:

  Additional arguments passed to `contrast`

## Value

An object of class `summary_emm` containing the multivariate test
results; or a list of the estimates and the tests if `show.ests` is
`TRUE`. The test results include the Hotelling \\T^2\\ statistic, \\F\\
ratios, degrees of freedom, and \\P\\ values.

## Note

If some interactions among the primary and `mult.name` factors are
absent, the covariance of the multivariate means is singular; this
situation is accommodated, but the result has reduced degrees of freedom
and a message is displayed. If there are other abnormal conditions such
as non-estimable results, estimates are shown as `NA`.

While designed primarily for testing contrasts, multivariate tests of
the mean vector itself can be implemented via `method = "identity")`
(see the examples).

## References

Hotelling, Harold (1931) "The generalization of Student's ratio",
*Annals of Mathematical Statistics* 2(3), 360–378.
doi:10.1214/aoms/1177732979

## Examples

``` r
MOats.lm <- lm(yield ~ Variety + Block, data = MOats)
MOats.emm <- emmeans(MOats.lm, ~ Variety | rep.meas)
mvcontrast(MOats.emm, "consec", show.ests = TRUE)  # mult.name defaults to rep.meas
#> $estimates
#> contrast = Marvellous - Golden Rain:
#>  rep.meas    estimate    SE df t.ratio p.value
#>  rep.meas0       6.67  7.84 10   0.851  0.4148
#>  rep.meas0.2    10.00  9.34 10   1.071  0.3093
#>  rep.meas0.4     2.50 12.30 10   0.203  0.8430
#>  rep.meas0.6     2.00 10.30 10   0.194  0.8503
#> 
#> contrast = Victory - Marvellous:
#>  rep.meas    estimate    SE df t.ratio p.value
#>  rep.meas0     -15.17  7.84 10  -1.936  0.0817
#>  rep.meas0.2   -18.83  9.34 10  -2.017  0.0713
#>  rep.meas0.4    -6.33 12.30 10  -0.515  0.6177
#>  rep.meas0.6    -8.33 10.30 10  -0.807  0.4385
#> 
#> Results are averaged over the levels of: Block 
#> 
#> $tests
#>  contrast                 T.square df1 df2 F.ratio p.value
#>  Marvellous - Golden Rain    3.082   4   7   0.539  0.9174
#>  Victory - Marvellous        9.181   4   7   1.607  0.4726
#> 
#> P value adjustment: sidak 
#> 

# Test each mean against a specified null vector
mvcontrast(MOats.emm, "identity", name = "Variety", 
           null = c(80, 100, 120, 140), adjust = "none")
#>  Variety     T.square df1 df2 F.ratio p.value
#>  Golden Rain   10.001   4   7   1.750  0.2430
#>  Marvellous    26.628   4   7   4.660  0.0377
#>  Victory       10.232   4   7   1.791  0.2352
#> 
# (Note 'name' is passed to contrast() and overrides default name "contrast")

# 'mult.name' need not refer to a multivariate response
mvcontrast(MOats.emm, "trt.vs.ctrl1", mult.name = "Variety")
#>  contrast                T.square df1 df2 F.ratio p.value
#>  rep.meas0.2 - rep.meas0   21.498   3   8   5.733  0.0634
#>  rep.meas0.4 - rep.meas0   37.578   3   8  10.021  0.0131
#>  rep.meas0.6 - rep.meas0  104.700   3   8  27.920  0.0004
#> 
#> P value adjustment: sidak 
```
