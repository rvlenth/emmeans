# Utilities and options for emmeans

## Contents

1.  [Updating an `emmGrid` object](#update)
2.  [Setting options](#options)
    1.  [Setting and viewing defaults](#defaults)
    2.  [Optimal digits to display](#digits)
    3.  [Startup options](#startup)
3.  [Combining and subsetting `emmGrid` objects](#rbind)
4.  [Accessing results to use elsewhere](#data)
5.  [Adding grouping factors](#groups)
6.  [Re-labeling and re-leveling an `emmGrid`](#relevel)

[Index of all vignette
topics](https://rvlenth.github.io/emmeans/articles/vignette-topics.md)

## Updating an `emmGrid` object

Several internal settings are saved when functions like
[`ref_grid()`](https://rvlenth.github.io/emmeans/reference/ref_grid.md),
[`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.md),
[`contrast()`](https://rvlenth.github.io/emmeans/reference/contrast.md),
etc. are run. Those settings can be manipulated via the
[`update()`](https://rdrr.io/r/stats/update.html) method for `emmGrid`s.
To illustrate, consider the `pigs` dataset and model yet again:

``` r
pigs.lm <- lm(log(conc) ~ source + factor(percent), data = pigs)
pigs.emm <- emmeans(pigs.lm, "source")
pigs.emm
```

``` ro
##  source emmean     SE df lower.CL upper.CL
##  fish     3.39 0.0367 23     3.32     3.47
##  soy      3.67 0.0374 23     3.59     3.74
##  skim     3.80 0.0394 23     3.72     3.88
## 
## Results are averaged over the levels of: percent 
## Results are given on the log (not the response) scale. 
## Confidence level used: 0.95
```

We see confidence intervals but not tests, by default. This happens as a
result of internal settings in `pigs.emm.s` that are passed to
[`summary()`](https://rdrr.io/r/base/summary.html) when the object is
displayed. If we are going to work with this object a lot, we might want
to change its internal settings rather than having to rely on explicitly
calling [`summary()`](https://rdrr.io/r/base/summary.html) with several
arguments. If so, just update the internal settings to what is desired;
for example:

``` r
pigs.emm.s <- update(pigs.emm, infer = c(TRUE, TRUE), null = log(35),
                     calc = c(n = ".wgt."))
pigs.emm.s
```

``` ro
##  source emmean     SE df  n lower.CL upper.CL null t.ratio p.value
##  fish     3.39 0.0367 23 10     3.32     3.47 3.56  -4.385  0.0002
##  soy      3.67 0.0374 23 10     3.59     3.74 3.56   2.988  0.0066
##  skim     3.80 0.0394 23  9     3.72     3.88 3.56   6.130 <0.0001
## 
## Results are averaged over the levels of: percent 
## Results are given on the log (not the response) scale. 
## Confidence level used: 0.95
```

Note that by adding of `calc`, we have set a default to calculate and
display the sample size when the object is summarized. See
[`help("update.emmGrid")`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md)
for details on the keywords that can be changed. Mostly, they are the
same as the names of arguments in the functions that construct these
objects.

Of course, we can always get what we want via calls to
[`test()`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md),
[`confint()`](https://rdrr.io/r/stats/confint.html) or
[`summary()`](https://rdrr.io/r/base/summary.html) with appropriate
arguments. But the [`update()`](https://rdrr.io/r/stats/update.html)
function is more useful in sophisticated manipulations of objects, or
called implicitly via the `...` or `options` argument in
[`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.md)
and other functions. Those options are passed to
[`update()`](https://rdrr.io/r/stats/update.html) just before the object
is returned. For example, we could have done the above update within the
[`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.md)
call as follows (results are not shown because they are the same as
before):

``` r
emmeans(pigs.lm, "source", infer = c(TRUE, TRUE), null = log(35),
        calc = c(n = ".wgt."))
```

[Back to contents](#contents)

## Setting options

Speaking of the `options` argument, note that the default in
[`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.md) is
`options = get_emm_option("emmeans")`. Let’s see what that is:

``` r
get_emm_option("emmeans")
```

``` ro
## $infer
## [1]  TRUE FALSE
```

So, by default, confidence intervals, but not tests, are displayed when
the result is summarized. The reverse is true for results of
[`contrast()`](https://rvlenth.github.io/emmeans/reference/contrast.md)
(and also the default for
[`pairs()`](https://rdrr.io/r/graphics/pairs.html) which calls
[`contrast()`](https://rvlenth.github.io/emmeans/reference/contrast.md)):

``` r
get_emm_option("contrast")
```

``` ro
## $infer
## [1] FALSE  TRUE
```

There are also defaults for a newly constructed reference grid:

``` r
get_emm_option("ref_grid")
```

``` ro
## $is.new.rg
## [1] TRUE
## 
## $infer
## [1] FALSE FALSE
```

The default is to display neither intervals nor tests when summarizing.
In addition, the flag `is.new.rg` is set to `TRUE`, and that is why one
sees a [`str()`](https://rdrr.io/r/utils/str.html) listing rather than a
summary as the default when the object is simply shown by typing its
name at the console.

### Setting and viewing defaults

The user may have other preferences. She may want to see both intervals
and tests whenever contrasts are produced; and perhaps she also wants to
always default to the response scale when transformations or links are
present. We can change the defaults by setting the corresponding
options; and that is done via the
[`emm_options()`](https://rvlenth.github.io/emmeans/reference/emm_options.md)
function:

``` r
emm_options(emmeans = list(type = "response"),
            contrast = list(infer = c(TRUE, TRUE)))
```

Now, new
[`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.md)
results and contrasts follow the new defaults:

``` r
pigs.anal.p <- emmeans(pigs.lm, consec ~ percent)
pigs.anal.p
```

``` ro
## $emmeans
##  percent response   SE df lower.CL upper.CL
##        9     31.4 1.28 23     28.8     34.1
##       12     37.5 1.44 23     34.7     40.6
##       15     39.0 1.70 23     35.6     42.7
##       18     42.3 2.24 23     37.9     47.2
## 
## Results are averaged over the levels of: source 
## Confidence level used: 0.95 
## Intervals are back-transformed from the log scale 
## 
## $contrasts
##  contrast              ratio     SE df lower.CL upper.CL null t.ratio p.value
##  percent12 / percent9   1.20 0.0671 23    1.038     1.38    1   3.202  0.0110
##  percent15 / percent12  1.04 0.0604 23    0.896     1.20    1   0.650  0.8613
##  percent18 / percent15  1.09 0.0750 23    0.911     1.29    1   1.194  0.5201
## 
## Results are averaged over the levels of: source 
## Confidence level used: 0.95 
## Conf-level adjustment: mvt method for 3 estimates 
## Intervals are back-transformed from the log scale 
## P value adjustment: mvt method for 3 tests 
## Tests are performed on the log scale
```

Observe that the contrasts “inherited” the `type = "response"` default
from the EMMs.

NOTE: Setting the above options does *not* change how existing `emmGrid`
objects are displayed; it only affects ones constructed in the future.

There is one more option – `summary` – that overrides all other display
defaults for both existing and future objects. For example, specifying
`emm_options(summary = list(infer = c(TRUE, TRUE)))` will result in both
intervals and tests being displayed, regardless of their internal
defaults, unless `infer` is explicitly specified in a call to
[`summary()`](https://rdrr.io/r/base/summary.html).

To temporarily revert to factory defaults in a single call to
[`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.md) or
[`contrast()`](https://rvlenth.github.io/emmeans/reference/contrast.md)
or [`pairs()`](https://rdrr.io/r/graphics/pairs.html), specify
`options = NULL` in the call. To reset everything to factory defaults
(which we do presently), null-out all of the **emmeans** package
options:

``` r
options(emmeans = NULL)
```

### Optimal digits to display

When an `emmGrid` object is summarized and displayed, the factory
default is to display it with just enough digits as is justified by the
standard errors or HPD intervals of the estimates displayed. You may use
the `"opt.digits"` option to change this. If it is `TRUE` (the default),
we display only enough digits as is justified (but at least 3). If it is
set to `FALSE`, the number of digits is set using the R system’s
default, `getOption("digits")`; this is often much more precision than
is justified. To illustrate, here is the summary of `pigs.emm` displayed
without optimizing digits. Compare it with the first summary in this
vignette.

``` r
emm_options(opt.digits = FALSE)
pigs.emm
```

``` ro
##  source   emmean         SE df lower.CL upper.CL
##  fish   3.394492 0.03668122 23 3.318612 3.470373
##  soy    3.667260 0.03744798 23 3.589793 3.744727
##  skim   3.796770 0.03938283 23 3.715300 3.878240
## 
## Results are averaged over the levels of: percent 
## Results are given on the log (not the response) scale. 
## Confidence level used: 0.95
```

``` r
emm_options(opt.digits = TRUE)  # revert to optimal digits
```

By the way, setting this option does *not* round the calculated values
computed by
[`summary.emmGrid()`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
or saved in a `summary)emm` object; it simply controls the precision
displayed by `print.summary_emm()`.

### Startup options

The options accessed by
[`emm_options()`](https://rvlenth.github.io/emmeans/reference/emm_options.md)
and
[`get_emm_option()`](https://rvlenth.github.io/emmeans/reference/emm_options.md)
are stored in a list named `emmeans` within R’s options environment.
Therefore, if you desire options other than the defaults provided on a
regular basis, this can be easily arranged by specifying them in your
startup script for R. For example, if you want to default to
Satterthwaite degrees of freedom for `lmer` models, and display
confidence intervals rather than tests for contrasts, your `.Rprofile`
file could contain the line

``` r
options(emmeans = list(lmer.df = "satterthwaite", 
                       contrast = list(infer = c(TRUE, FALSE))))
```

[Back to contents](#contents)

## Combining and subsetting `emmGrid` objects

Two or more `emmGrid` objects may be combined using the
[`rbind()`](https://rdrr.io/r/base/cbind.html) or `+` methods. The most
common reason (or perhaps the only good reason) to do this is to combine
EMMs or contrasts into one family for purposes of applying a
multiplicity adjustment to tests or intervals. A user may want to
combine the three pairwise comparisons of sources with the three
comparisons above of consecutive percents into a single family of six
tests with a suitable multiplicity adjustment. This is done quite
simply:

``` r
rbind(pairs(pigs.emm.s), pigs.anal.p[[2]])
```

``` ro
##  contrast              estimate     SE df t.ratio p.value
##  fish - soy             -0.2728 0.0529 23  -5.153  0.0002
##  fish - skim            -0.4023 0.0542 23  -7.428 <0.0001
##  soy - skim             -0.1295 0.0530 23  -2.442  0.1364
##  percent12 - percent9    0.1796 0.0561 23   3.202  0.0238
##  percent15 - percent12   0.0378 0.0582 23   0.650  1.0000
##  percent18 - percent15   0.0825 0.0691 23   1.194  1.0000
## 
## Results are averaged over some or all of the levels of: percent, source 
## Results are given on the log (not the response) scale. 
## P value adjustment: bonferroni method for 6 tests
```

The default adjustment is `"bonferroni"`; we could have specified
something different via the `adjust` argument. An equivalent way to
combine `emmGrid`s is via the addition operator. Any options may be
provided by [`update()`](https://rdrr.io/r/stats/update.html). Below, we
combine the same results into a family but ask for the “exact”
multiplicity adjustment.

``` r
update(pigs.anal.p[[2]] + pairs(pigs.emm.s), adjust = "mvt")
```

``` ro
##  contrast              ratio     SE df lower.CL upper.CL null t.ratio p.value
##  percent12 / percent9  1.197 0.0671 23    1.022    1.402    1   3.202  0.0213
##  percent15 / percent12 1.039 0.0604 23    0.881    1.224    1   0.650  0.9680
##  percent18 / percent15 1.086 0.0750 23    0.894    1.320    1   1.194  0.7306
##  fish / soy            0.761 0.0403 23    0.656    0.884    1  -5.153  0.0002
##  fish / skim           0.669 0.0362 23    0.574    0.779    1  -7.428 <0.0001
##  soy / skim            0.879 0.0466 23    0.756    1.020    1  -2.442  0.1111
## 
## Results are averaged over some or all of the levels of: source, percent 
## Confidence level used: 0.95 
## Conf-level adjustment: mvt method for 6 estimates 
## Intervals are back-transformed from the log scale 
## P value adjustment: mvt method for 6 tests 
## Tests are performed on the log scale
```

Also evident in comparing these results is that settings are obtained
from the first object combined. So in the second output, where they are
combined in reverse order, we get both confidence intervals and tests,
and transformation to the response scale.

###### 

To subset an `emmGrid` object, just use the subscripting operator `[]`.
For instance,

``` r
pigs.emm[2:3]
```

``` ro
##  source emmean     SE df lower.CL upper.CL
##  soy      3.67 0.0374 23     3.59     3.74
##  skim     3.80 0.0394 23     3.72     3.88
## 
## Results are averaged over the levels of: percent 
## Results are given on the log (not the response) scale. 
## Confidence level used: 0.95
```

## Accessing results to use elsewhere

Sometimes, users want to use the results of an analysis (say, an
[`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.md)
call) in other computations. The
[`summary()`](https://rdrr.io/r/base/summary.html) method creates a
`summary_emm` object that inherits from the `data.frame` class; so one
may use the variables therein just as those in a data frame.

An `emmGrid` object has its own internal structure and we can’t directly
access the values we see displayed. If follow-up computations are
needed, use [`summary()`](https://rdrr.io/r/base/summary.html) (or
[`confint()`](https://rdrr.io/r/stats/confint.html) or
[`test()`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)),
creates a `summary_emm` object which inherits from `data.frame` – making
it possible to access the values. For illustration, let’s add the widths
of the confidence intervals in our example.

``` r
CIs <- confint(pigs.emm)
CIs$CI.width <- with(CIs, upper.CL - lower.CL)
CIs
```

``` ro
##  source emmean     SE df lower.CL upper.CL CI.width
##  fish     3.39 0.0367 23     3.32     3.47    0.152
##  soy      3.67 0.0374 23     3.59     3.74    0.155
##  skim     3.80 0.0394 23     3.72     3.88    0.163
## 
## Results are averaged over the levels of: percent 
## Results are given on the log (not the response) scale. 
## Confidence level used: 0.95
```

By the way, the values stored internally are kept to full precision,
more than is typically displayed:

``` r
CIs$emmean
```

``` ro
## [1] 3.394492 3.667260 3.796770
```

If you want to display more digits, specify so using the `print` method:

``` r
print(CIs, digits = 5)
```

``` ro
##  source emmean       SE df lower.CL upper.CL CI.width
##  fish   3.3945 0.036681 23   3.3186   3.4704  0.15176
##  soy    3.6673 0.037448 23   3.5898   3.7447  0.15493
##  skim   3.7968 0.039383 23   3.7153   3.8782  0.16294
## 
## Results are averaged over the levels of: percent 
## Results are given on the log (not the response) scale. 
## Confidence level used: 0.95
```

[Back to contents](#contents)

## Adding grouping factors

Sometimes, users want to group levels of a factor into a smaller number
of groups. Those groups may then be, say, averaged separately and
compared, or used as a `by` factor. The
[`add_grouping()`](https://rvlenth.github.io/emmeans/reference/manip-factors.md)
function serves this purpose. The function takes four arguments: the
object, the name of the grouping factor to be created, the name of the
reference factor that is being grouped, and a vector of level names of
the grouping factor corresponding to levels of the reference factor.
Suppose for example that we want to distinguish animal and non-animal
sources of protein in the `pigs` example:

``` r
pigs.emm.ss <- add_grouping(pigs.emm.s, "type", "source",
                            c("animal", "vegetable", "animal"))
str(pigs.emm.ss)
```

``` ro
## 'emmGrid' object with variables:
##     source = fish, soy, skim
##     type = animal, vegetable
## Nesting structure:  source %in% type
## Transformation: "log"
```

Note that the new object has a nesting structure (see more about this in
the [“messy-data”
vignette](https://rvlenth.github.io/emmeans/articles/messy-data.html#nesting)),
with the reference factor nested in the new grouping factor. Now we can
obtain means and comparisons for each group

``` r
emmeans(pigs.emm.ss, pairwise ~ type)
```

``` ro
## $emmeans
##  type      emmean     SE df  n lower.CL upper.CL
##  animal      3.60 0.0267 23 19     3.54     3.65
##  vegetable   3.67 0.0374 23 10     3.59     3.74
## 
## Results are averaged over the levels of: percent, source 
## Results are given on the log (not the response) scale. 
## Confidence level used: 0.95 
## 
## $contrasts
##  contrast           estimate     SE df t.ratio p.value
##  animal - vegetable  -0.0716 0.0455 23  -1.573  0.1295
## 
## Results are averaged over the levels of: percent, source 
## Results are given on the log (not the response) scale.
```

[Back to contents](#contents)

## Re-labeling or re-leveling an `emmGrid`

Sometimes it is desirable to re-label the rows of an `emmGrid`, or cast
it in terms of other factor(s). This can be done via the `levels`
argument in [`update()`](https://rdrr.io/r/stats/update.html).

As an example, sometimes a fitted model has a treatment factor that
comprises combinations of other factors. In subsequent analysis, we may
well want to break it down into the individual factors’ contributions.
Consider, for example, the `warpbreaks` data provided with R. We will
define a single factor and fit a non homogeneous-variance model:

``` r
warp <- transform(warpbreaks, treat = interaction(wool, tension))
library(nlme)
warp.gls <- gls(breaks ~ treat, weights = varIdent(form = ~ 1|treat), data = warp)
( warp.emm <- emmeans(warp.gls, "treat") )
```

``` ro
##  treat emmean   SE   df lower.CL upper.CL
##  A.L     44.6 6.03 7.97     30.6     58.5
##  B.L     28.2 3.29 8.00     20.6     35.8
##  A.M     24.0 2.89 8.00     17.3     30.7
##  B.M     28.8 3.14 8.00     21.5     36.0
##  A.H     24.6 3.42 8.00     16.7     32.5
##  B.H     18.8 1.63 8.00     15.0     22.5
## 
## Degrees-of-freedom method: satterthwaite 
## Confidence level used: 0.95
```

But now we want to re-cast this `emmGrid` into one that has separate
factors for `wool` and `tension`. We can do this as follows:

``` r
warp.fac <- update(warp.emm, levels = list(
                wool = c("A", "B"), tension = c("L", "M", "H")))
str(warp.fac)
```

``` ro
## 'emmGrid' object with variables:
##     wool = A, B
##     tension = L, M, H
```

So now we can do various contrasts involving the separate factors:

``` r
contrast(warp.fac, "consec", by = "wool")
```

``` ro
## wool = A:
##  contrast estimate   SE   df t.ratio p.value
##  M - L     -20.556 6.69 11.4  -3.074  0.0203
##  H - M       0.556 4.48 15.6   0.124  0.9899
## 
## wool = B:
##  contrast estimate   SE   df t.ratio p.value
##  M - L       0.556 4.55 16.0   0.122  0.9881
##  H - M     -10.000 3.54 12.0  -2.824  0.0269
## 
## Degrees-of-freedom method: satterthwaite 
## P value adjustment: mvt method for 2 tests
```

Note: When re-leveling to more than one factor, you have to be careful
to anticipate that the levels will be expanded using
[`expand.grid()`](https://rdrr.io/r/base/expand.grid.html): the first
factor in the list varies the fastest and the last varies the slowest.
That was the case in our example, but in others, it may not be. Had the
levels of `treat` been ordered as `A.L, A.M, A.H, B.L, B.M, B.H`, then
we would have had to specify the levels of `tension` first and the
levels of `wool` second.

[Back to contents](#contents)

[Index of all vignette
topics](https://rvlenth.github.io/emmeans/articles/vignette-topics.md)
