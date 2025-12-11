# Summaries, predictions, intervals, and tests for `emmGrid` objects

These are the primary methods for obtaining numerical or tabular results
from an `emmGrid` object. `summary.emmGrid` is the general function for
summarizing `emmGrid` objects. It also serves as the print method for
these objects; so for convenience,
[`summary()`](https://rdrr.io/r/base/summary.html) arguments may be
included in calls to functions such as
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md) and
[`contrast`](https://rvlenth.github.io/emmeans/reference/contrast.md)
that construct `emmGrid` objects. Note that by default, summaries for
Bayesian models are diverted to
[`hpd.summary`](https://rvlenth.github.io/emmeans/reference/hpd.summary.md).

## Usage

``` r
# S3 method for class 'emmGrid'
summary(object, infer, level, adjust, by,
  cross.adjust = "none", type, df, calc, null, delta, side, frequentist,
  bias.adjust = get_emm_option("back.bias.adj"), sigma, ...)

# S3 method for class 'emmGrid'
confint(object, parm, level = 0.95, ...)

test(object, null, ...)

# S3 method for class 'emmGrid'
test(object, null = 0, joint = FALSE, verbose = FALSE,
  rows, by, status = FALSE, ...)

# S3 method for class 'emmGrid'
predict(object, type, interval = c("none", "confidence",
  "prediction"), level = 0.95,
  bias.adjust = get_emm_option("back.bias.adj"), sigma, ...)

# S3 method for class 'emmGrid'
as.data.frame(x, row.names = NULL, optional,
  check.names = TRUE, destroy.annotations = FALSE, ...)

# S3 method for class 'summary_emm'
x[..., as.df = FALSE]
```

## Arguments

- object:

  An object of class `"emmGrid"` (see
  [emmGrid-class](https://rvlenth.github.io/emmeans/reference/emmGrid-class.md))

- infer:

  A vector of one or two logical values. The first determines whether
  confidence intervals are displayed, and the second determines whether
  *t* tests and *P* values are displayed. If only one value is provided,
  it is used for both.

- level:

  Numerical value between 0 and 1. Confidence level for confidence
  intervals, if `infer[1]` is `TRUE`.

- adjust:

  Character value naming the method used to adjust \\p\\ values or
  confidence limits; or to adjust comparison arrows in `plot`. See the
  P-value adjustments section below.

- by:

  Character name(s) of variables to use for grouping into separate
  tables. This affects the family of tests considered in adjusted *P*
  values.

- cross.adjust:

  Character: \\p\\-value adjustment method to additionally apply
  *across* the `by` groups. See the section on P-value adjustments for
  details.

- type:

  Character: type of prediction desired. This only has an effect if
  there is a known transformation or link function. `"response"`
  specifies that the inverse transformation be applied. `"mu"` (or
  equivalently, `"unlink"`) is usually the same as `"response"`, but in
  the case where the model has both a link function and a response
  transformation, only the link part is back-transformed. Other valid
  values are `"link"`, `"lp"`, and `"linear.predictor"`; these are
  equivalent, and request that results be shown for the linear
  predictor, with no back-transformation. The default is `"link"`,
  unless the `"predict.type"` option is in force; see
  [`emm_options`](https://rvlenth.github.io/emmeans/reference/emm_options.md),
  and also the section below on transformations and links.

- df:

  Numeric. If non-missing, a constant number of degrees of freedom to
  use in constructing confidence intervals and *P* values (`NA`
  specifies asymptotic results).

- calc:

  Named list of character value(s) or formula(s). The expressions in
  `char` are evaluated and appended to the summary, just after the `df`
  column. The expression may include any names up through `df` in the
  summary, any additional names in `object@grid` (such as `.wgt.` or
  `.offset.`), or any earlier elements of `calc`.

- null:

  Numeric. Null hypothesis value(s), on the linear-predictor scale,
  against which estimates are tested. May be a single value used for
  all, or a numeric vector of length equal to the number of tests in
  each family (i.e., `by` group in the displayed table).

- delta:

  Numeric value (on the linear-predictor scale). If zero, ordinary tests
  of significance are performed. If positive, this specifies a threshold
  for testing equivalence (using the TOST or two-one-sided-test method),
  non-inferiority, or non-superiority, depending on `side`. See Details
  for how the test statistics are defined.

- side:

  Numeric or character value specifying whether the test is left-tailed
  (`-1`, `"-"`, `"<"`, `"left"`, or `"nonsuperiority"`); right-tailed
  (`1`, `"+"`, `">"`, `"right"`, or `"noninferiority"`); or two-sided
  (`0`, `2`, `"!="`, `"two-sided"`, `"both"`, `"equivalence"`, or
  `"="`). See the special section below for more details.

- frequentist:

  Ignored except if a Bayesian model was fitted. If missing or `FALSE`,
  the object is passed to
  [`hpd.summary`](https://rvlenth.github.io/emmeans/reference/hpd.summary.md).
  Otherwise, a logical value of `TRUE` will have it return a frequentist
  summary.

- bias.adjust:

  Logical value for whether to adjust for bias in back-transforming
  (`type = "response"`). This requires a valid value of `sigma` to exist
  in the object or be specified.

- sigma:

  Error SD assumed for bias correction (when `type = "response"` and a
  transformation is in effect), or for constructing prediction
  intervals. If not specified, `object@misc$sigma` is used, and a
  warning is issued if it is not found or not valid. *Note:* `sigma` may
  be a vector, but be careful that it correctly corresponds (perhaps
  after recycling) to the order of the reference grid.

- ...:

  Optional arguments such as `scheffe.rank` (see “P-value adjustments”).
  In `confint.emmGrid`, `predict.emmGrid`, and `test.emmGrid`, these
  arguments are passed to `summary.emmGrid`.

- parm:

  (Required argument for `confint` methods, but not used)

- joint:

  Logical value. If `FALSE`, the arguments are passed to
  `summary.emmGrid` with `infer=c(FALSE, TRUE)`. If `joint = TRUE`, a
  joint test of the hypothesis L beta = null is performed, where L is
  `linfct(object)` and beta is the vector of fixed effects estimated by
  `object@betahat`. This will be either an *F* test or a chi-square
  (Wald) test depending on whether degrees of freedom are available. See
  also
  [`joint_tests`](https://rvlenth.github.io/emmeans/reference/joint_tests.md).

- verbose:

  Logical value. If `TRUE` and `joint = TRUE`, a table of the effects
  being tested is printed.

- rows:

  Integer values. The rows of L to be tested in the joint test. If
  missing, all rows of L are used. If not missing, `by` variables are
  ignored.

- status:

  logical. If `TRUE`, a `note` column showing status flags (for rank
  deficiencies and estimability issues) is displayed even when empty. If
  `FALSE`, the column is included only if there are such issues.

- interval:

  Type of interval desired (partial matching is allowed): `"none"` for
  no intervals, otherwise confidence or prediction intervals with given
  arguments, via `confint.emmGrid`. Note: prediction intervals are not
  available unless the model family is `"gaussian"`.

- x:

  object of the given class

- row.names:

  passed to [`as.data.frame`](https://rdrr.io/r/base/as.data.frame.html)

- optional:

  required argument, but ignored in `as.data.frame.emmGrid`

- check.names:

  passed to [`data.frame`](https://rdrr.io/r/base/data.frame.html)

- destroy.annotations:

  Logical value. If `FALSE`, an object of class `summary_emm` is
  returned (which inherits from `data.frame`), but if displayed, details
  like confidence levels, P-value adjustments, transformations, etc. are
  also shown. But unlike the result of `summary`, the number of digits
  displayed is obtained from `getOption("digits")` rather than using the
  optimal digits algorithm we usually use. Thus, it is formatted more
  like a regular data frame, but with any annotations and groupings
  still intact. If `TRUE` (not recommended), a “plain vanilla” data
  frame is returned, based on `row.names` and `check.names`.

- as.df:

  Logical value. With `x[..., as.df = TRUE]`, the result is object is
  coerced to a [`data.frame`](https://rdrr.io/r/base/data.frame.html)
  before the subscripting is applied. With `as.df = FALSE`, the result
  is returned as a `summary_emm` object when possible.

## Value

`summary.emmGrid`, `confint.emmGrid`, and `test.emmGrid` return an
object of class `"summary_emm"`, which is an extension of
[`data.frame`](https://rdrr.io/r/base/data.frame.html) but with a
special `print` method that displays it with custom formatting. For
models fitted using MCMC methods, the call is diverted to
[`hpd.summary`](https://rvlenth.github.io/emmeans/reference/hpd.summary.md)
(with `prob` set to `level`, if specified); one may alternatively use
general MCMC summarization tools with the results of `as.mcmc`.

`predict` returns a vector of predictions for each row of `object@grid`.

The `as.data.frame` method returns an object that inherits from
`"data.frame"`.

## Details

`confint.emmGrid` is equivalent to
`summary.emmGrid with infer = c(TRUE, FALSE)`. The function
`test.emmGrid`, when called with `joint = FALSE`, is equivalent to
`summary.emmGrid` with `infer = c(FALSE, TRUE)`.

With `joint = TRUE`, `test.emmGrid` calculates the Wald test of the
hypothesis `linfct %*% bhat = null`, where `linfct` and `bhat` refer to
slots in `object` (possibly subsetted according to `by` or `rows`). An
error is thrown if any row of `linfct` is non-estimable. It is
permissible for the rows of `linfct` to be linearly dependent, as long
as `null == 0`, in which case a reduced set of contrasts is tested.
Linear dependence and nonzero `null` cause an error. The returned object
has an additional `"est.fcns"` attribute, which is a list of the linear
functions associated with the joint test.

## Note

In doing testing and a transformation and/or link is in force, any
`null` and/or `delta` values specified must always be on the scale of
the linear predictor, regardless of the setting for \`type\`. If
`type = "response"`, the null value displayed in the summary table will
be back-transformed from the value supplied by the user. But the
displayed `delta` will not be changed, because there (often) is not a
natural way to back-transform it.

When we have `type = "response"`, and `bias.adj = TRUE`, the `null`
value displayed in the output is both back-transformed and
bias-adjusted, leading to a rather non-intuitive-looking null value.
However, since the tests themselves are performed on the link scale,
this is the response value at which a \*P\* value of 1 would be
obtained.

The default `show` method for `emmGrid` objects (with the exception of
newly created reference grids) is `print(summary())`. Thus, with
ordinary usage of
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md) and
such, it is unnecessary to call `summary` unless there is a need to
specify other than its default options.

If a data frame is needed, `summary`, `confint`, and `test` serve this
need. `as.data.frame` routes to `summary` by default; calling it with
`destroy.annotations = TRUE` is not recommended for exactly that reason.
If you want to see more digits in the output, use
`print(summary(object), digits = ...)`; and if you *always* want to see
more digits, use `emm_options(opt.digits = FALSE)`.

## Defaults

The `misc` slot in `object` may contain default values for `by`, `calc`,
`infer`, `level`, `adjust`, `type`, `null`, `side`, and `delta`. These
defaults vary depending on the code that created the object. The
[`update`](https://rdrr.io/r/stats/update.html) method may be used to
change these defaults. In addition, any options set using
`emm_options(summary = ...)` will trump those stored in the object's
`misc` slot.

## Transformations and links

With `type = "response"`, the transformation assumed can be found in
`object@misc$tran`, and its label, for the summary is in
`object@misc$inv.lbl`. Any \\t\\ or \\z\\ tests are still performed on
the scale of the linear predictor, not the inverse-transformed one.
Similarly, confidence intervals are computed on the linear-predictor
scale, then inverse-transformed.

Be aware that only univariate transformations and links are supported in
this way. Some multivariate transformations are supported by
[`mvregrid`](https://rvlenth.github.io/emmeans/reference/mvregrid.md).

## Bias adjustment when back-transforming

When `bias.adjust` is `TRUE`, then back-transformed estimates are
adjusted by adding \\0.5 h''(u)\sigma^2\\, where \\h\\ is the inverse
transformation and \\u\\ is the linear predictor. This is based on a
second-order Taylor expansion. There are better or exact adjustments for
certain specific cases, and these may be incorporated in future updates.

Note: In certain models, e.g., those with non-gaussian families, `sigma`
is initialized as `NA`, and so by default, bias adjustment is skipped
and a warning is issued. You may override this by specifying a value for
`sigma`. However, *with ordinary generalized linear models, bias
adjustment is inappropriate* and you should not try to do it. With GEEs
and GLMMs, you probably should *not* use `sigma(model)`, and instead you
should create an appropriate value using the estimated random effects,
e.g., from `VarCorr(model)`. An example is provided in the
“transformations” vignette.

*A word of caution:* This bias-adjustment method is merely a one-term
correction, and it multiplies `sigma^2`. When `sigma` is reasonably
small relative to the scale of the response, this is fine; but when
`sigma` is large, this bias adjustment can produce bizarre and
out-of-range results. Consider for example a logit or probit model where
the back-transformed response is constrained to the interval \[0, 1\];
if `sigma` is more than 1, the adjustment can easily pull the estimate
outside of these constraints.

## P-value adjustments

The `adjust` argument specifies a multiplicity adjustment for tests or
confidence intervals. This adjustment always is applied *separately* to
each table or sub-table that you see in the printed output (see
[`rbind.emmGrid`](https://rvlenth.github.io/emmeans/reference/rbind.emmGrid.md)
for how to combine tables). If there are non-estimable cases in a `by`
group, those cases are *excluded* before determining the adjustment;
that means there could be different adjustments in different groups.

The valid values of `adjust` are as follows:

- `"tukey"`:

  Uses the Studentized range distribution with the number of means in
  the family. (Available for two-sided cases only.)

- `"scheffe"`:

  Computes \\p\\ values from the \\F\\ distribution, according to the
  Scheffe critical value of \\\sqrt{rF(\alpha; r, d)}\\, where \\d\\ is
  the error degrees of freedom and \\r\\ is the rank of the set of
  linear functions under consideration. By default, the value of `r` is
  computed from `linfct(object)` for each by group; however, if the user
  specifies an argument matching `scheffe.rank`, its value will be used
  instead. Ordinarily, if there are \\k\\ means involved, then \\r = k -
  1\\ for a full set of contrasts involving all \\k\\ means, and \\r =
  k\\ for the means themselves. (The Scheffe adjustment is available for
  two-sided cases only.)

- `"sidak"`:

  Makes adjustments as if the estimates were independent (a conservative
  adjustment in many cases).

- `"bonferroni"`:

  Multiplies \\p\\ values, or divides significance levels by the number
  of estimates. This is a conservative adjustment.

- `"dunnettx"`:

  Uses our own*ad hoc* approximation to the Dunnett distribution for a
  family of estimates having pairwise correlations of \\0.5\\ (as is
  true when comparing treatments with a control with equal sample
  sizes). The accuracy of the approximation improves with the number of
  simultaneous estimates, and is much faster than `"mvt"`. (Available
  for two-sided cases only.)

- `"mvt"`:

  Uses the multivariate \\t\\ distribution to assess the probability or
  critical value for the maximum of \\k\\ estimates. This method
  produces the same \\p\\ values and intervals as the default `summary`
  or `confint` methods to the results of
  [`as.glht`](https://rvlenth.github.io/emmeans/reference/glht-support.md).
  In the context of pairwise comparisons or comparisons with a control,
  this produces “exact” Tukey or Dunnett adjustments, respectively.
  However, the algorithm (from the mvtnorm package) uses a Monte Carlo
  method, so results are not exactly repeatable unless the same
  random-number seed is used (see
  [`set.seed`](https://rdrr.io/r/base/Random.html)). As the family size
  increases, the required computation time will become noticeable or
  even intolerable, making the `"tukey"`, `"dunnettx"`, or others more
  attractive.

- `"none"`:

  Makes no adjustments to the \\p\\ values.

For tests, not confidence intervals, the Bonferroni-inequality-based
adjustment methods in
[`p.adjust`](https://rdrr.io/r/stats/p.adjust.html) are also available
(currently, these include `"holm"`, `"hochberg"`, `"hommel"`,
`"bonferroni"`, `"BH"`, `"BY"`, `"fdr"`, and `"none"`). If a
`p.adjust.methods` method other than `"bonferroni"` or `"none"` is
specified for confidence limits, the straight Bonferroni adjustment is
used instead. Also, if an adjustment method is not appropriate (e.g.,
using `"tukey"` with one-sided tests, or with results that are not
pairwise comparisons), a more appropriate method (usually `"sidak"`) is
substituted.

In some cases, confidence and \\p\\-value adjustments are only
approximate – especially when the degrees of freedom or standard errors
vary greatly within the family of tests. The `"mvt"` method is always
the correct one-step adjustment, but it can be very slow. One may use
[`as.glht`](https://rvlenth.github.io/emmeans/reference/glht-support.md)
with methods in the multcomp package to obtain non-conservative
multi-step adjustments to tests.

*Warning:* Non-estimable cases are *included* in the family to which
adjustments are applied. You may wish to subset the object using the
`[]` operator to work around this problem.

The `cross.adjust` argument is a way of specifying a multiplicity
adjustment across the `by` groups (otherwise by default, each group is
treated as a separate family in regards to multiplicity adjustments). It
applies only to \\p\\ values. Valid options are one of the
`p.adjust.methods` or `"sidak"`. This argument is ignored unless it is
other than `"none"`, there is more than one `by` group, and they are all
the same size. Under those conditions, we first use `adjust` to
determine the within-group adjusted \\p\\ values. Imagine each group's
adjusted \\p\\ values arranged in side-by-side columns, thus forming a
matrix with the number of columns equal to the number of `by` groups.
Then we use the `cross.adjust` method to further adjust the adjusted
\\p\\ values in each row of this matrix. Note that an *overall*
Bonferroni (or Sidak) adjustment is obtainable by specifying *both*
`adjust` and `cross.adjust` as `"bonferroni"` (or `"sidak"`). However,
less conservative (but yet conservative) overall adjustments are
available when it is possible to use an “exact” within-group method
(e.g., `adjust = "tukey"` for pairwise comparisons) and `cross.adjust`
as a conservative adjustment. \[`cross.adjust` methods other than
`"none"`, `"bonferroni"`, or `"sidak"` do not seem advisable, but other
`p.adjust` methods are available if you can make sense of them.\]

## Tests of significance, nonsuperiority, noninferiority, or equivalence

When `delta = 0`, test statistics are the usual tests of significance.
They are of the form `(estimate - null)/SE`. Notationally:

- Significance:

  \\H_0: \theta = \theta_0\\ versus  
  \\H_1: \theta \< \theta_0\\ (left-sided), or  
  \\H_1: \theta \> \theta_0\\ (right-sided), or  
  \\H_1: \theta \ne \theta_0\\ (two-sided)  
  The test statistic is  
  \\t = (Q - \theta_0)/SE\\  
  where \\Q\\ is our estimate of \\\theta\\; then left, right, or
  two-sided \\p\\ values are produced, depending on `side`.

When `delta` is positive, the test statistic depends on `side` as
follows.

- Left-sided (nonsuperiority):

  \\H_0: \theta \ge \theta_0 + \delta\\ versus \\H_1: \theta \<
  \theta_0 + \delta\\  
  \\t = (Q - \theta_0 - \delta)/SE\\  
  The \\p\\ value is the lower-tail probability.

- Right-sided (noninferiority):

  \\H_0: \theta \le \theta_0 - \delta\\ versus \\H_1: \theta \>
  \theta_0 - \delta\\  
  \\t = (Q - \theta_0 + \delta)/SE\\  
  The \\p\\ value is the upper-tail probability.

- Two-sided (equivalence):

  \\H_0: \|\theta - \theta_0\| \ge \delta\\ versus \\H_1: \|\theta -
  \theta_0\| \< \delta\\  
  \\t = (\|Q - \theta_0\| - \delta)/SE\\  
  The \\p\\ value is the *lower*-tail probability.  
  Note that \\t\\ is the maximum of \\t\_{nonsup}\\ and
  \\-t\_{noninf}\\. This is equivalent to choosing the less significant
  result in the two-one-sided-test (TOST) procedure.

## Non-estimable cases

When the model is rank-deficient, each row `x` of `object`'s `linfct`
slot is checked for estimability. If `sum(x*bhat)` is found to be
non-estimable, then the string `NonEst` is displayed for the estimate,
and associated statistics are set to `NA`. The estimability check is
performed using the orthonormal basis `N` in the `nbasis` slot for the
null space of the rows of the model matrix. Estimability fails when
\\\|\|Nx\|\|^2 / \|\|x\|\|^2\\ exceeds `tol`, which by default is
`1e-8`. You may change it via
[`emm_options`](https://rvlenth.github.io/emmeans/reference/emm_options.md)
by setting `estble.tol` to the desired value.

See the warning above that non-estimable cases are still included when
determining the family size for *P*-value adjustments.

## Warning about potential misuse of P values

Some in the statistical and scientific community argue that the term
“statistical significance” should be completely abandoned, and that
criteria such as “p \< 0.05” never be used to assess the importance of
an effect. These practices can be too misleading and are prone to abuse.
See [the “basics”
vignette](https://rvlenth.github.io/emmeans/doc/basics.html#pvalues) for
more discussion.

## See also

[`hpd.summary`](https://rvlenth.github.io/emmeans/reference/hpd.summary.md)

## Examples

``` r
warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
warp.emm <- emmeans(warp.lm, ~ tension | wool)
warp.emm    # implicitly runs 'summary'
#> wool = A:
#>  tension emmean   SE df lower.CL upper.CL
#>  L         44.6 3.65 48     37.2     51.9
#>  M         24.0 3.65 48     16.7     31.3
#>  H         24.6 3.65 48     17.2     31.9
#> 
#> wool = B:
#>  tension emmean   SE df lower.CL upper.CL
#>  L         28.2 3.65 48     20.9     35.6
#>  M         28.8 3.65 48     21.4     36.1
#>  H         18.8 3.65 48     11.4     26.1
#> 
#> Confidence level used: 0.95 

confint(warp.emm, by = NULL, level = .90)
#>  tension wool emmean   SE df lower.CL upper.CL
#>  L       A      44.6 3.65 48     38.4     50.7
#>  M       A      24.0 3.65 48     17.9     30.1
#>  H       A      24.6 3.65 48     18.4     30.7
#>  L       B      28.2 3.65 48     22.1     34.3
#>  M       B      28.8 3.65 48     22.7     34.9
#>  H       B      18.8 3.65 48     12.7     24.9
#> 
#> Confidence level used: 0.9 

# --------------------------------------------------------------
pigs.lm <- lm(log(conc) ~ source + factor(percent), data = pigs)
pigs.emm <- emmeans(pigs.lm, "percent", type = "response")
summary(pigs.emm)    # (inherits type = "response")
#>  percent response   SE df lower.CL upper.CL
#>        9     31.4 1.28 23     28.8     34.1
#>       12     37.5 1.44 23     34.7     40.6
#>       15     39.0 1.70 23     35.6     42.7
#>       18     42.3 2.24 23     37.9     47.2
#> 
#> Results are averaged over the levels of: source 
#> Confidence level used: 0.95 
#> Intervals are back-transformed from the log scale 
summary(pigs.emm, calc = c(n = ".wgt."))  # Show sample size
#>  percent response   SE df n lower.CL upper.CL
#>        9     31.4 1.28 23 8     28.8     34.1
#>       12     37.5 1.44 23 9     34.7     40.6
#>       15     39.0 1.70 23 7     35.6     42.7
#>       18     42.3 2.24 23 5     37.9     47.2
#> 
#> Results are averaged over the levels of: source 
#> Confidence level used: 0.95 
#> Intervals are back-transformed from the log scale 

# For which percents is EMM non-inferior to 35, based on a 10% threshold?
# Note the test is done on the log scale even though we have type = "response"
test(pigs.emm, null = log(35), delta = log(1.10), side = ">")
#>  percent response   SE df null t.ratio p.value
#>        9     31.4 1.28 23   35  -0.360  0.6390
#>       12     37.5 1.44 23   35   4.295  0.0001
#>       15     39.0 1.70 23   35   4.635 <0.0001
#>       18     42.3 2.24 23   35   5.384 <0.0001
#> 
#> Results are averaged over the levels of: source 
#> Statistics are tests of noninferiority with a threshold of 0.09531 
#> P values are right-tailed 
#> Tests are performed on the log scale 

con <- contrast(pigs.emm, "consec")
test(con)
#>  contrast              ratio     SE df null t.ratio p.value
#>  percent12 / percent9   1.20 0.0671 23    1   3.202  0.0111
#>  percent15 / percent12  1.04 0.0604 23    1   0.650  0.8613
#>  percent18 / percent15  1.09 0.0750 23    1   1.194  0.5201
#> 
#> Results are averaged over the levels of: source 
#> P value adjustment: mvt method for 3 tests 
#> Tests are performed on the log scale 

test(con, joint = TRUE)
#>  df1 df2 F.ratio p.value
#>    3  23   7.981  0.0008
#> 

# default Scheffe adjustment - rank = 3
summary(con, infer = c(TRUE, TRUE), adjust = "scheffe")
#>  contrast              ratio     SE df lower.CL upper.CL null t.ratio p.value
#>  percent12 / percent9   1.20 0.0671 23    1.011     1.42    1   3.202  0.0343
#>  percent15 / percent12  1.04 0.0604 23    0.872     1.24    1   0.650  0.9344
#>  percent18 / percent15  1.09 0.0750 23    0.882     1.34    1   1.194  0.7027
#> 
#> Results are averaged over the levels of: source 
#> Confidence level used: 0.95 
#> Conf-level adjustment: scheffe method with rank 3 
#> Intervals are back-transformed from the log scale 
#> P value adjustment: scheffe method with rank 3 
#> Tests are performed on the log scale 

# Consider as some of many possible contrasts among the six cell means
summary(con, infer = c(TRUE, TRUE), adjust = "scheffe", scheffe.rank = 5)
#>  contrast              ratio     SE df lower.CL upper.CL null t.ratio p.value
#>  percent12 / percent9   1.20 0.0671 23    0.976     1.47    1   3.202  0.1090
#>  percent15 / percent12  1.04 0.0604 23    0.841     1.28    1   0.650  0.9940
#>  percent18 / percent15  1.09 0.0750 23    0.845     1.40    1   1.194  0.9165
#> 
#> Results are averaged over the levels of: source 
#> Confidence level used: 0.95 
#> Conf-level adjustment: scheffe method with rank 5 
#> Intervals are back-transformed from the log scale 
#> P value adjustment: scheffe method with rank 5 
#> Tests are performed on the log scale 

# Show estimates to more digits
print(test(con), digits = 7)
#>  contrast                 ratio         SE df null t.ratio p.value
#>  percent12 / percent9  1.196684 0.06710564 23    1   3.202  0.0109
#>  percent15 / percent12 1.038570 0.06042501 23    1   0.650  0.8614
#>  percent18 / percent15 1.085945 0.07499759 23    1   1.194  0.5201
#> 
#> Results are averaged over the levels of: source 
#> P value adjustment: mvt method for 3 tests 
#> Tests are performed on the log scale 

# --------------------------------------------------------------
# Cross-adjusting P values
prs <- pairs(warp.emm)   # pairwise comparisons of tension, by wool
test(prs, adjust = "tukey", cross.adjust = "bonferroni")
#> wool = A:
#>  contrast estimate   SE df t.ratio p.value
#>  L - M      20.556 5.16 48   3.986  0.0013
#>  L - H      20.000 5.16 48   3.878  0.0018
#>  M - H      -0.556 5.16 48  -0.108  1.0000
#> 
#> wool = B:
#>  contrast estimate   SE df t.ratio p.value
#>  L - M      -0.556 5.16 48  -0.108  1.0000
#>  L - H       9.444 5.16 48   1.831  0.3407
#>  M - H      10.000 5.16 48   1.939  0.2777
#> 
#> P value adjustment: tukey method for comparing a family of 3 estimates 
#> Cross-group P-value adjustment: bonferroni method for 2 tests 

# Same comparisons taken as one big family (more conservative)
test(prs, adjust = "bonferroni", by = NULL)
#>  contrast wool estimate   SE df t.ratio p.value
#>  L - M    A      20.556 5.16 48   3.986  0.0014
#>  L - H    A      20.000 5.16 48   3.878  0.0019
#>  M - H    A      -0.556 5.16 48  -0.108  1.0000
#>  L - M    B      -0.556 5.16 48  -0.108  1.0000
#>  L - H    B       9.444 5.16 48   1.831  0.4396
#>  M - H    B      10.000 5.16 48   1.939  0.3504
#> 
#> P value adjustment: bonferroni method for 6 tests 
```
