# Set or change emmeans options

Use `emm_options` to set or change various options that are used in the
emmeans package. These options are set separately for different contexts
in which `emmGrid` objects are created, in a named list of option lists.

## Usage

``` r
emm_options(..., disable)

get_emm_option(x, default = emm_defaults[[x]])

with_emm_options(..., expr)

emm_defaults
```

## Format

An object of class `list` of length 25.

## Arguments

- ...:

  Option names and values (see Details)

- disable:

  If non-missing, this will reset all options to their defaults if
  `disable` tests `TRUE` (but first save them for possible later
  restoration). Otherwise, all previously saved options are restored.
  This is important for bug reporting; please see the section below on
  reproducible bugs. When `disable` is specified, the other arguments
  are ignored.

- x:

  Character value - the name of an option to be queried

- default:

  Value to return if `x` is not found

- expr:

  Expression to evaluate. If missing, the last element of `...` is used.

## Value

`emm_options` returns the current options (same as the result of
`getOption("emmeans")`) – invisibly, unless called with no arguments.

`get_emm_option` returns the currently stored option for `x`, or its
default value if not found.

`with_emm_options()` temporarily sets the options in `...`, then
evaluates `try(expr)` and returns the result.

## Details

emmeans's options are stored as a list in the system option `"emmeans"`.
Thus, `emm_options(foo = bar)` is the same as
`options(emmeans = list(..., foo = bar))` where `...` represents any
previously existing options. The list `emm_defaults` contains the
default values in case the corresponding element of system option
`emmeans` is `NULL`.

Currently, the following main list entries are supported:

- `ref_grid`:

  A named `list` of defaults for objects created by
  [`ref_grid`](https://rvlenth.github.io/emmeans/reference/ref_grid.md).
  This could affect other objects as well. For example, if `emmeans` is
  called with a fitted model object, it calls `ref_grid` and this option
  will affect the resulting `emmGrid` object.

- `emmeans`:

  A named `list` of defaults for objects created by
  [`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md) or
  [`emtrends`](https://rvlenth.github.io/emmeans/reference/emtrends.md).

- `contrast`:

  A named `list` of defaults for objects created by
  [`contrast.emmGrid`](https://rvlenth.github.io/emmeans/reference/contrast.md)
  or
  [`pairs.emmGrid`](https://rvlenth.github.io/emmeans/reference/contrast.md).

- `summary`:

  A named `list` of defaults used by the methods
  [`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md),
  [`predict.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md),
  [`test.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md),
  [`confint.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md),
  and [`emmip`](https://rvlenth.github.io/emmeans/reference/emmip.md).
  The only option that can affect the latter four is `"predict.type"`.

- `allow.na.levs`:

  A logical value that if `TRUE` (the default), allows `NA` to be among
  the levels of a factor. Older versions of emmeans did not allow this.
  So if problems come up (say in a messy dataset that includes
  incomplete cases), try setting this to `FALSE`.

- `sep`:

  A character value to use as a separator in labeling factor
  combinations. Such labels are potentially used in several places such
  as
  [`contrast`](https://rvlenth.github.io/emmeans/reference/contrast.md)
  and
  [`plot.emmGrid`](https://rvlenth.github.io/emmeans/reference/plot.md)
  when combinations of factors are compared or plotted. The default is
  `" "`.

- `parens`:

  Character vector that determines which labels are parenthesized when
  they are contrasted. The first element is a regular expression, and
  the second and third elements are used as left and right parentheses.
  See details for the `parens` argument in
  [`contrast`](https://rvlenth.github.io/emmeans/reference/contrast.md).
  The default will parenthesize labels containing the four arithmetic
  operators, using round parentheses.

- `cov.keep`:

  The default value of `cov.keep` in
  [`ref_grid`](https://rvlenth.github.io/emmeans/reference/ref_grid.md).
  Defaults to `"2"`, i.e., two-level covariates are treated like
  factors.

- `graphics.engine`:

  A character value matching `c("ggplot", "lattice")`, setting the
  default engine to use in
  [`emmip`](https://rvlenth.github.io/emmeans/reference/emmip.md) and
  [`plot.emmGrid`](https://rvlenth.github.io/emmeans/reference/plot.md).
  Defaults to `"ggplot"`.

- `gg.theme`:

  An integer. Set it to 1 or 2 if you want the appearance of plots used
  in emmeans version 1.x.x or 2.x.x. Or set it to the actual theme that
  you want to use, e.g., `emm_options(gg.theme = ggplot2::theme_dark())`

- `msg.interaction`:

  A logical value controlling whether or not a message is displayed when
  `emmeans` averages over a factor involved in an interaction. It is
  probably not appropriate to do this, unless the interaction is weak.
  Defaults to `TRUE`.

- `msg.nesting`:

  A logical value controlling whether or not to display a message when a
  nesting structure is auto-detected. The existence of such a structure
  affects computations of EMMs. Sometimes, a nesting structure is
  falsely detected – namely when a user has omitted some main effects
  but included them in interactions. This does not change the model fit,
  but it produces a different parameterization that is picked up when
  the reference grid is constructed. Defaults to `TRUE`.

- `rg.limit`:

  An integer value setting a limit on the number of rows in a newly
  constructed reference grid. This is checked based on the number of
  levels of the factors involved; but it excludes the levels of any
  multivariate responses because those are not yet known. The reference
  grid consists of all possible combinations of the predictors, and this
  can become huge if there are several factors. An error is thrown if
  this limit is exceeded. One can use the `nuisance` argument of
  [`ref_grid`](https://rvlenth.github.io/emmeans/reference/ref_grid.md)
  to collapse on nuisance factors, thus making the grid smaller.
  Defaults to 10,000.

- `simplify.names`:

  A logical value controlling whether to simplify (when possible) names
  in the model formula that refer to datasets – for example, should we
  simplify a predictor name like “`data$trt`” to just “`trt`”? Defaults
  to `TRUE`.

- `opt.digits`:

  A logical value controlling the precision with which summaries are
  printed. If `TRUE` (default), the number of digits displayed is just
  enough to reasonably distinguish estimates from the ends of their
  confidence intervals; but always at least 3 digits. If `FALSE`, the
  system value `getOption("digits")` is used.

- `pval.digits`:

  An integer indicating the precision with which p-values are printed.
  Integers between `2` and `6` are acceptable; the default is `4`.
  Floating point numbers are truncated. If a number outside the accepted
  range is used, the nearest acceptable integer will be used instead.
  Non-numeric arguments are ignored.

- `back.bias.adj`:

  A logical value controlling whether we try to adjust bias when
  back-transforming. If `FALSE`, we use naive back transformation. If
  `TRUE` *and `sigma` is available and valid*, a second-order adjustment
  is applied to estimate the mean on the response scale. A warning is
  issued if no valid `sigma` is available

- `enable.submodel`:

  A logical value. If `TRUE`, enables support for selected model classes
  to implement the `submodel` option. If `FALSE`, this support is
  disabled. Setting this option to `FALSE` could save excess memory
  consumption.

Some other options have more specific purposes:

- `estble.tol`:

  Tolerance for determining estimability in rank-deficient cases. If
  absent, the value in `emm_defaults$estble.tol)` is used.

- `save.ref_grid`:

  Logical value of `TRUE` if you wish the latest reference grid created
  to be saved in `.Last.ref_grid`. The default is `FALSE`.

- Options for `lme4::lmerMod` models:

  Options `lmer.df`, `disable.pbkrtest`, `pbkrtest.limit`,
  `disable.lmerTest`, and `lmerTest.limit` options affect how degrees of
  freedom are computed for `lmerMod` objects produced by the lme4
  package). See that section of the "models" vignette for details.

- `post.ci.method`:

  method for estimating confidence intervals for Bayesian models.
  Options are `"HPD"` (default) or `"quantile"`.

## Reproducible bugs

Most options set display attributes and such that are not likely to be
associated with bugs in the code. However, some other options (e.g.,
`cov.keep`) are essentially configuration settings that may affect
how/whether the code runs, and the settings for these options may cause
subtle effects that may be hard to reproduce. Therefore, when sending a
bug report, please create a reproducible example and make sure the bug
occurs with all options set at their defaults. This is done by preceding
it with `emm_options(disable = TRUE)`.

By the way, `disable` works like a stack (LIFO buffer), in that
`disable = TRUE` is equivalent to
`emm_options(saved.opts = emm_options())` and
`emm_options(disable = FALSE)` is equivalent to
`options(emmeans = get_emm_option("saved.opts"))`. To completely erase
all options, use `options(emmeans = NULL)`

## Startup options

emmeans's options are kept in the system options in a list names
`emmeans`. Thus, if the user wants certain options set on startup, this
can be arranged by setting the `emmeans` option in one's startup script.
For example, the user's `.Rprofile` file could contain:

        options(emmeans = list(summary = list(predict.type = "response"),
                               lmer.df = "satterthwaite"))

## See also

[`update.emmGrid`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md)

## Examples

``` r
if (FALSE) { # \dontrun{
emm_options(ref_grid = list(level = .90),
            contrast = list(infer = c(TRUE,FALSE)),
            estble.tol = 1e-6)
# Sets default confidence level to .90 for objects created by ref.grid
# AS WELL AS emmeans called with a model object (since it creates a 
# reference grid). In addition, when we call 'contrast', 'pairs', etc.,
# confidence intervals rather than tests are displayed by default.
} # }

if (FALSE) { # \dontrun{
emm_options(disable.pbkrtest = TRUE)
# This forces use of asymptotic methods for lmerMod objects.
# Set to FALSE or NULL to re-enable using pbkrtest.
} # }

# See tolerance being used for determining estimability
get_emm_option("estble.tol")
#> [1] 1e-08

if (FALSE) { # \dontrun{
# Set all options to their defaults
emm_options(disable = TRUE)
# ... and perhaps follow with code for a minimal reproducible bug,
#     which may include emm_options() clls if they are pertinent ...

# restore options that had existed previously
emm_options(disable = FALSE)
} # }


# Temporarily alter the formatting of labels in an interaction contrast:
pigs.lm <- lm(inverse(conc) ~ source * factor(percent), data = pigs)
with_emm_options(parens = c("\\-", "<", ">"), sep = ": ",
    contrast(contrast(ref_grid(pigs.lm), "consec", by = "percent"),
        "consec", by = NULL))
#>  contrast                                          estimate      SE df t.ratio
#>  <skim - soy: percent9> - <soy - fish: percent9>    0.00919 0.00460 17   1.997
#>  <soy - fish: percent12> - <skim - soy: percent9>  -0.00656 0.00361 17  -1.817
#>  <skim - soy: percent12> - <soy - fish: percent12>  0.00531 0.00442 17   1.201
#>  <soy - fish: percent15> - <skim - soy: percent12> -0.00498 0.00383 17  -1.301
#>  <skim - soy: percent15> - <soy - fish: percent15>  0.00186 0.00477 17   0.391
#>  <soy - fish: percent18> - <skim - soy: percent15> -0.00265 0.00460 17  -0.577
#>  <skim - soy: percent18> - <soy - fish: percent18>  0.00108 0.00722 17   0.150
#>  p.value
#>   0.2827
#>   0.3713
#>   0.7520
#>   0.6892
#>   0.9983
#>   0.9866
#>   1.0000
#> 
#> P value adjustment: mvt method for 7 tests 
```
