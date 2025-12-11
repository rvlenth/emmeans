# Confidence intervals and tests in emmeans

## Contents

This vignette describes various ways of summarizing `emmGrid` objects.

1.  [`summary()`, `confint()`, and `test()`](#summary)
2.  [Back-transforming to response scale](#tran) (See also the
    [“transformations”
    vignette](https://rvlenth.github.io/emmeans/articles/transformations.md))
3.  [Multiplicity adjustments](#adjust)
4.  [Using “by” variables](#byvars)
5.  [Joint (omnibus) tests](#joint)
6.  [Testing equivalence, noninferiority, nonsuperiority](#equiv)
7.  Graphics (in [“basics”
    vignette](https://rvlenth.github.io/emmeans/articles/basics.html#plots))

[Index of all vignette
topics](https://rvlenth.github.io/emmeans/articles/vignette-topics.md)

## `summary()`, `confint()`, and `test()`

The most important method for `emmGrid` objects is
[`summary()`](https://rdrr.io/r/base/summary.html). For one thing, it is
called by default when you display an
[`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.md)
result. The [`summary()`](https://rdrr.io/r/base/summary.html) function
has a lot of options, and the detailed documentation via
[`help("summary.emmGrid")`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
is worth a look.

For ongoing illustrations, let’s re-create some of the objects in the
[“basics”
vignette](https://rvlenth.github.io/emmeans/articles/basics.md) for the
`pigs` example:

``` r
mod4 <- lm(inverse(conc) ~ source + factor(percent), data = pigs)
RG <- ref_grid(mod4)
EMM.source <- emmeans(RG, "source")
```

Just `summary(<object>)` by itself will produce a summary that varies
somewhat according to context. It does this by setting different
defaults for the `infer` argument, which consists of two logical values,
specifying confidence intervals and tests, respectively. \[The exception
is models fitted using MCMC methods, where
[`summary()`](https://rdrr.io/r/base/summary.html) is diverted to the
[`hpd.summary()`](https://rvlenth.github.io/emmeans/reference/hpd.summary.md)
function, a preferable summary for many Bayesians.\]

The summary of a newly made reference grid will show just estimates and
standard errors, but not confidence intervals or tests (that is,
`infer = c(FALSE, FALSE)`). The summary of an
[`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.md)
result, as we see above, will have intervals, but no tests (i.e.,
`infer = c(TRUE, FALSE)`); and the result of a
[`contrast()`](https://rvlenth.github.io/emmeans/reference/contrast.md)
call (see [comparisons and
contrasts](https://rvlenth.github.io/emmeans/articles/comparisons.md))
will show test statistics and *P* values, but not intervals (i.e.,
`infer = c(FALSE, TRUE)`). There are courtesy methods
[`confint()`](https://rdrr.io/r/stats/confint.html) and
[`test()`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
that just call [`summary()`](https://rdrr.io/r/base/summary.html) with
the appropriate `infer` setting; for example,

``` r
test(EMM.source)
```

``` ro
##  source emmean       SE df t.ratio p.value
##  fish   0.0337 0.000926 23  36.380 <0.0001
##  soy    0.0257 0.000945 23  27.141 <0.0001
##  skim   0.0229 0.000994 23  22.989 <0.0001
## 
## Results are averaged over the levels of: percent 
## Results are given on the inverse (not the response) scale.
```

It is not particularly useful, though, to test these EMMs against the
default of zero – which is why tests are not usually shown. It makes a
lot more sense to test them against some target concentration, say 40.
And suppose we want to do a one-sided test to see if the concentration
is greater than 40. Remembering that the response is inverse-transformed
in this model, and that the inverse transformation reverses the
direction of comparisons, so that a *right*-tailed test on the `conc`
scale becomes a *left*-tailed test on the `inverse(conc)` scale,

``` r
test(EMM.source, null = inverse(40), side = "<")
```

``` ro
##  source emmean       SE df  null t.ratio p.value
##  fish   0.0337 0.000926 23 0.025   9.383  1.0000
##  soy    0.0257 0.000945 23 0.025   0.697  0.7535
##  skim   0.0229 0.000994 23 0.025  -2.156  0.0209
## 
## Results are averaged over the levels of: percent 
## Results are given on the inverse (not the response) scale. 
## P values are left-tailed
```

It is also possible to add calculated columns to the summary, via the
`calc` argument. The calculations can include any columns up through
`df` in the summary, as well as any variable in the object’s `grid`
slot. Among the latter are usually weights in a column named `.wgt.`,
and we can use that to include sample size in the summary:

``` r
confint(EMM.source, calc = c(n = ~.wgt.))
```

``` ro
##  source emmean       SE df  n lower.CL upper.CL
##  fish   0.0337 0.000926 23 10   0.0318   0.0356
##  soy    0.0257 0.000945 23 10   0.0237   0.0276
##  skim   0.0229 0.000994 23  9   0.0208   0.0249
## 
## Results are averaged over the levels of: percent 
## Results are given on the inverse (not the response) scale. 
## Confidence level used: 0.95
```

[Back to Contents](#contents)

## Back-transforming

Transformations and link functions are supported in several ways in
**emmeans**, making this a complex topic worthy of [its own
vignette](https://rvlenth.github.io/emmeans/articles/transformations.md).
Here, we show just the most basic approach. Namely, specifying the
argument `type = "response"` will cause the displayed results to be
back-transformed to the response scale, when a transformation or link
function is incorporated in the model. For example, let’s try the
preceding
[`test()`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
call again:

``` r
test(EMM.source, null = inverse(40), side = "<", type = "response")
```

``` ro
##  source response    SE df null t.ratio p.value
##  fish       29.7 0.816 23   40   9.383  1.0000
##  soy        39.0 1.440 23   40   0.697  0.7535
##  skim       43.8 1.900 23   40  -2.156  0.0209
## 
## Results are averaged over the levels of: percent 
## P values are left-tailed 
## Tests are performed on the inverse scale
```

Note what changes and what doesn’t change. In the
[`test()`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
call, we *still* use the 1/40 as the null value; `null` must always be
specified on the linear-prediction scale, in this case the inverse. In
the output, the displayed estimates, as well as the `null` value, are
shown back-transformed. As well, the standard errors are altered (using
the delta method). However, the *t* ratios and *P* values are identical
to the preceding results. That is, the tests themselves are still
conducted on the linear-predictor scale (as is noted in the output).

Similar statements apply to confidence intervals on the response scale:

``` r
confint(EMM.source, side = "<", level = .90, type = "response")
```

``` ro
##  source response    SE df lower.CL upper.CL
##  fish       29.7 0.816 23     28.6      Inf
##  soy        39.0 1.440 23     37.2      Inf
##  skim       43.8 1.900 23     41.4      Inf
## 
## Results are averaged over the levels of: percent 
## Confidence level used: 0.9 
## Intervals are back-transformed from the inverse scale
```

With `side = "<"`, an *upper* confidence limit is computed on the
inverse scale, then that limit is back-transformed to the response
scale; and since `inverse` reverses everything, those upper confidence
limits become lower ones on the response scale. (We have also
illustrated how to change the confidence level.)

[Back to Contents](#contents)

## Multiplicity adjustments

Both tests and confidence intervals may be adjusted for simultaneous
inference. Such adjustments ensure that the confidence coefficient for a
whole set of intervals is at least the specified level, or to control
for multiplicity in a whole family of tests. This is done via the
`adjust` argument. For
[`ref_grid()`](https://rvlenth.github.io/emmeans/reference/ref_grid.md)
and
[`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.md)
results, the default is `adjust = "none"`. For most
[`contrast()`](https://rvlenth.github.io/emmeans/reference/contrast.md)
results, `adjust` is often something else, depending on what type of
contrasts are created. For example, pairwise comparisons default to
`adjust = "tukey"`, i.e., the Tukey HSD method. The
[`summary()`](https://rdrr.io/r/base/summary.html) function sometimes
*changes* `adjust` if it is inappropriate. For example, with

``` r
confint(EMM.source, adjust = "tukey")
```

    ## Note: adjust = "tukey" was changed to "sidak"
    ## because "tukey" is only appropriate for one set of pairwise comparisons

``` ro
##  source emmean       SE df lower.CL upper.CL
##  fish   0.0337 0.000926 23   0.0313   0.0361
##  soy    0.0257 0.000945 23   0.0232   0.0281
##  skim   0.0229 0.000994 23   0.0203   0.0254
## 
## Results are averaged over the levels of: percent 
## Results are given on the inverse (not the response) scale. 
## Confidence level used: 0.95 
## Conf-level adjustment: sidak method for 3 estimates
```

the adjustment is changed to the Sidak method because the Tukey
adjustment is inappropriate unless you are doing pairwise comparisons.

An adjustment method that is usually appropriate is Bonferroni; however,
it can be quite conservative. Using `adjust = "mvt"` is the closest to
being the “exact” all-around method “single-step” method, as it uses the
multivariate *t* distribution (and the **mvtnorm** package) with the
same covariance structure as the estimates to determine the adjustment.
However, this comes at high computational expense as the computations
are done using simulation techniques. For a large set of tests (and
especially confidence intervals), the computational lag becomes
noticeable if not intolerable.

For tests, `adjust` increases the *P* values over those otherwise
obtained with `adjust = "none"`. Compare the following adjusted tests
with the unadjusted ones previously computed.

``` r
test(EMM.source, null = inverse(40), side = "<", adjust = "bonferroni")
```

``` ro
##  source emmean       SE df  null t.ratio p.value
##  fish   0.0337 0.000926 23 0.025   9.383  1.0000
##  soy    0.0257 0.000945 23 0.025   0.697  1.0000
##  skim   0.0229 0.000994 23 0.025  -2.156  0.0627
## 
## Results are averaged over the levels of: percent 
## Results are given on the inverse (not the response) scale. 
## P value adjustment: bonferroni method for 3 tests 
## P values are left-tailed
```

[Back to Contents](#contents)

## “By” variables

Sometimes you want to break a summary down into smaller pieces; for this
purpose, the `by` argument in
[`summary()`](https://rdrr.io/r/base/summary.html) is useful. For
example,

``` r
confint(RG, by = "source")
```

``` ro
## source = fish:
##  percent prediction      SE df lower.CL upper.CL
##        9     0.0385 0.00135 23   0.0357   0.0413
##       12     0.0333 0.00125 23   0.0307   0.0359
##       15     0.0326 0.00138 23   0.0297   0.0354
##       18     0.0304 0.00138 23   0.0275   0.0332
## 
## source = soy:
##  percent prediction      SE df lower.CL upper.CL
##        9     0.0305 0.00126 23   0.0279   0.0331
##       12     0.0253 0.00124 23   0.0227   0.0278
##       15     0.0245 0.00128 23   0.0219   0.0272
##       18     0.0223 0.00162 23   0.0190   0.0257
## 
## source = skim:
##  percent prediction      SE df lower.CL upper.CL
##        9     0.0277 0.00127 23   0.0251   0.0303
##       12     0.0225 0.00125 23   0.0199   0.0250
##       15     0.0217 0.00139 23   0.0189   0.0246
##       18     0.0195 0.00163 23   0.0162   0.0229
## 
## Results are given on the inverse (not the response) scale. 
## Confidence level used: 0.95
```

If there is also an `adjust` in force when `by` variables are used, by
default, the adjustment is made *separately* on each `by` group; e.g.,
in the above, we would be adjusting for sets of 4 intervals, not all 12
together (but see “cross-adjustments” below.)

There can be a `by` specification in
[`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.md)
(or equivalently, a `|` in the formula); and if so, it is passed on to
[`summary()`](https://rdrr.io/r/base/summary.html) and used unless
overridden by another `by`. Here are examples, not run:

``` r
emmeans(mod4, ~ percent | source)     ### same results as above
summary(.Last.value, by = "percent")       ### grouped the other way
```

Specifying `by = NULL` will remove all grouping.

### Adjustments across `by` groups

As was mentioned, each `by` group is regarded as a separate family with
regards to the `adjust` procedure. For example, consider a model with
interaction for the `warpbreaks` data, and construct pairwise
comparisons of `tension` by `wool`:

``` r
warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
warp.pw <- pairs(emmeans(warp.lm, ~ tension | wool))
warp.pw
```

``` ro
## wool = A:
##  contrast estimate   SE df t.ratio p.value
##  L - M      20.556 5.16 48   3.986  0.0007
##  L - H      20.000 5.16 48   3.878  0.0009
##  M - H      -0.556 5.16 48  -0.108  0.9936
## 
## wool = B:
##  contrast estimate   SE df t.ratio p.value
##  L - M      -0.556 5.16 48  -0.108  0.9936
##  L - H       9.444 5.16 48   1.831  0.1704
##  M - H      10.000 5.16 48   1.939  0.1389
## 
## P value adjustment: tukey method for comparing a family of 3 estimates
```

We have two sets of 3 comparisons, and the (default) Tukey adjustment is
made *separately* in each group.

However, sometimes we want the multiplicity adjustment to be broader.
This broadening can be done in two ways. One is to remove the `by`
variable, which then treats all results as one family. In our example:

``` r
test(warp.pw, by = NULL, adjust = "bonferroni")
```

``` ro
##  contrast wool estimate   SE df t.ratio p.value
##  L - M    A      20.556 5.16 48   3.986  0.0014
##  L - H    A      20.000 5.16 48   3.878  0.0019
##  M - H    A      -0.556 5.16 48  -0.108  1.0000
##  L - M    B      -0.556 5.16 48  -0.108  1.0000
##  L - H    B       9.444 5.16 48   1.831  0.4396
##  M - H    B      10.000 5.16 48   1.939  0.3504
## 
## P value adjustment: bonferroni method for 6 tests
```

This accomplishes the goal of putting all the comparisons in one family
of 6 comparisons. Note that the Tukey adjustment may not be used here
because we no longer have *one* set of pairwise comparisons.

An alternative is to specify `cross.adjust`, which specifies an
additional adjustment method to apply to corresponding sets of
within-group adjusted *P* values:

``` r
test(warp.pw, adjust = "tukey", cross.adjust = "bonferroni")
```

``` ro
## wool = A:
##  contrast estimate   SE df t.ratio p.value
##  L - M      20.556 5.16 48   3.986  0.0013
##  L - H      20.000 5.16 48   3.878  0.0018
##  M - H      -0.556 5.16 48  -0.108  1.0000
## 
## wool = B:
##  contrast estimate   SE df t.ratio p.value
##  L - M      -0.556 5.16 48  -0.108  1.0000
##  L - H       9.444 5.16 48   1.831  0.3407
##  M - H      10.000 5.16 48   1.939  0.2777
## 
## P value adjustment: tukey method for comparing a family of 3 estimates 
## Cross-group P-value adjustment: bonferroni method for 2 tests
```

These adjustments are less conservative than the previous result, but it
is still a conservative adjustment to the set of 6 tests. Had we also
specified `adjust = "bonferroni"`, we would have obtained the same
adjusted *P* values as we obtained with `by = NULL`.

### Simple comparisons

There is also a `simple` argument for
[`contrast()`](https://rvlenth.github.io/emmeans/reference/contrast.md)
that is in essence the inverse of `by`; the contrasts are run using
everything *except* the specified variables as `by` variables. To
illustrate, let’s consider the model for `pigs` that includes the
interaction (so that the levels of one factor compare differently at
levels of the other factor).

``` r
mod5 <- lm(inverse(conc) ~ source * factor(percent), data = pigs)
RG5 <- ref_grid(mod5)
contrast(RG5, "consec", simple = "percent")
```

``` ro
## source = fish:
##  contrast               estimate      SE df t.ratio p.value
##  percent12 - percent9  -6.64e-03 0.00285 17  -2.328  0.0833
##  percent15 - percent12 -6.68e-05 0.00285 17  -0.023  1.0000
##  percent18 - percent15 -1.40e-03 0.00285 17  -0.489  0.9283
## 
## source = soy:
##  contrast               estimate      SE df t.ratio p.value
##  percent12 - percent9  -4.01e-03 0.00255 17  -1.572  0.3168
##  percent15 - percent12  2.61e-04 0.00255 17   0.102  0.9993
##  percent18 - percent15 -2.18e-03 0.00361 17  -0.605  0.8872
## 
## source = skim:
##  contrast               estimate      SE df t.ratio p.value
##  percent12 - percent9  -5.26e-03 0.00255 17  -2.061  0.1401
##  percent15 - percent12 -2.86e-03 0.00285 17  -1.001  0.6526
##  percent18 - percent15 -3.76e-03 0.00383 17  -0.982  0.6650
## 
## Note: contrasts are still on the inverse scale. Consider using
##       regrid() if you want contrasts of back-transformed estimates. 
## P value adjustment: mvt method for 3 tests
```

In fact, we may do *all* one-factor comparisons by specifying
`simple = "each"`. This typically produces a lot of output, so use it
with care.

[Back to Contents](#contents)

## Joint tests

From the above, we already know how to test individual results. For
pairwise comparisons (details in [the “comparisons”
vignette](https://rvlenth.github.io/emmeans/articles/comparisons.md)),
we might do

``` r
PRS.source <- pairs(EMM.source)
PRS.source
```

``` ro
##  contrast    estimate      SE df t.ratio p.value
##  fish - soy   0.00803 0.00134 23   6.009 <0.0001
##  fish - skim  0.01083 0.00137 23   7.922 <0.0001
##  soy - skim   0.00280 0.00134 23   2.092  0.1136
## 
## Results are averaged over the levels of: percent 
## Note: contrasts are still on the inverse scale. Consider using
##       regrid() if you want contrasts of back-transformed estimates. 
## P value adjustment: tukey method for comparing a family of 3 estimates
```

But suppose we want an *omnibus* test that all these comparisons are
zero. Easy enough, using the `joint` argument in `test` (note: the
`joint` argument is *not* available in
[`summary()`](https://rdrr.io/r/base/summary.html); only in
[`test()`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)):

``` r
test(PRS.source, joint = TRUE)
```

``` ro
##  df1 df2 F.ratio p.value note
##    2  23  34.009 <0.0001  d  
## 
## d: df1 reduced due to linear dependence
```

Notice that there are three comparisons, but only 2 d.f. for the test,
as cautioned in the message.

The test produced with `joint = TRUE` is a “type III” test (assuming the
default equal weights are used to obtain the EMMs). See more on these
types of tests for higher-order effects in the [“interactions” vignette
section on
contrasts](https://rvlenth.github.io/emmeans/articles/interactions.html#contrasts).

For convenience, there is also a
[`joint_tests()`](https://rvlenth.github.io/emmeans/reference/joint_tests.md)
function that performs joint tests of contrasts among each term in a
model or `emmGrid` object.

``` r
joint_tests(RG5)
```

``` ro
##  model term     df1 df2 F.ratio p.value
##  source           2  17  30.309 <0.0001
##  percent          3  17   8.441  0.0012
##  source:percent   6  17   0.481  0.8135
```

The tests of main effects are of families of contrasts; those for
interaction effects are for interaction contrasts. These results are
essentially the same as a “Type-III ANOVA”, but may differ in situations
where there are empty cells or other non-estimability issues, or if
generalizations are present such as unequal weighting. (Another
distinction is that sums of squares and mean squares are not shown; that
is because these really are tests of contrasts among predictions, and
they may or may not correspond to model sums of squares.)

One may use `by` variables with `joint_tests`. For example:

``` r
joint_tests(RG5, by = "source")
```

``` ro
## source = fish:
##  model term df1 df2 F.ratio p.value
##  percent      3  17   2.967  0.0614
## 
## source = soy:
##  model term df1 df2 F.ratio p.value
##  percent      3  17   1.376  0.2840
## 
## source = skim:
##  model term df1 df2 F.ratio p.value
##  percent      3  17   4.835  0.0130
```

In some models, it is possible to specify `submodel = "type2"`, thereby
obtaining something akin to a Type II analysis of variance. See the
[messy-data
vignette](https://rvlenth.github.io/emmeans/articles/messy-data.html#type2submodel)
for an example.

[Back to Contents](#contents)

## Testing equivalence, noninferiority, and nonsuperiority

The `delta` argument in
[`summary()`](https://rdrr.io/r/base/summary.html) or
[`test()`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
allows the user to specify a threshold value to use in a test of
equivalence, non-inferiority, or non-superiority. An equivalence test is
kind of a backwards significance test, where small *P* values are
associated with small differences relative to a specified threshold
value `delta`. The help page for `summary.emmGrid` gives the details of
these tests. Suppose in the present example, we consider two sources to
be equivalent if they are within 0.005 of each other. We can test this
as follows:

``` r
test(PRS.source, delta = 0.005, adjust = "none")
```

``` ro
##  contrast    estimate      SE df t.ratio p.value
##  fish - soy   0.00803 0.00134 23   2.268  0.9835
##  fish - skim  0.01083 0.00137 23   4.266  0.9999
##  soy - skim   0.00280 0.00134 23  -1.641  0.0572
## 
## Results are averaged over the levels of: percent 
## Note: contrasts are still on the inverse scale. Consider using
##       regrid() if you want contrasts of back-transformed estimates. 
## Statistics are tests of equivalence with a threshold of 0.005 
## P values are left-tailed
```

Using the 0.005 threshold, the *P* value is quite small for comparing
soy and skim, providing some statistical evidence that their difference
is enough smaller than the threshold to consider them equivalent.

[Back to Contents](#contents)

## Graphics

Graphical displays of `emmGrid` objects are described in the [“basics”
vignette](https://rvlenth.github.io/emmeans/articles/basics.html#plots)

[Index of all vignette
topics](https://rvlenth.github.io/emmeans/articles/vignette-topics.md)
