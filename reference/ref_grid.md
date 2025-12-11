# Create a reference grid from a fitted model

Using a fitted model object, determine a reference grid for which
estimated marginal means are defined. The resulting `ref_grid` object
encapsulates all the information needed to calculate EMMs and make
inferences on them.

## Usage

``` r
ref_grid(object, at, cov.reduce = mean,
  cov.keep = get_emm_option("cov.keep"), mult.names, mult.levs,
  options = get_emm_option("ref_grid"), data, df, type, regrid, nesting,
  offset, sigma, counterfactuals, nuisance = character(0), non.nuisance,
  wt.nuis = "equal", rg.limit = get_emm_option("rg.limit"), ...)
```

## Arguments

- object:

  An object produced by a supported model-fitting function, such as
  `lm`. Many models are supported. See
  [[`vignette("models", "emmeans")`](https://rvlenth.github.io/emmeans/articles/models.md)](https://rvlenth.github.io/emmeans/doc/models.md).

- at:

  Optional named list of levels for the corresponding variables

- cov.reduce:

  A function, logical value, or formula; or a named list of these. Each
  covariate *not* specified in `cov.keep` or `at` is reduced according
  to these specifications. See the section below on “Using `cov.reduce`
  and `cov.keep`”.

- cov.keep:

  Character vector: names of covariates that are *not* to be reduced;
  these are treated as factors and used in weighting calculations.
  `cov.keep` may also include integer value(s), and if so, the maximum
  of these is used to set a threshold such that any covariate having no
  more than that many unique values is automatically included in
  `cov.keep`.

- mult.names:

  Character value: the name(s) to give to the pseudo-factor(s) whose
  levels delineate the elements of a multivariate response. If this is
  provided, it overrides the default name(s) used for `class(object)`
  when it has a multivariate response (e.g., the default is `"rep.meas"`
  for `"mlm"` objects).

- mult.levs:

  A named list of levels for the dimensions of a multivariate response.
  If there is more than one element, the combinations of levels are
  used, in [`expand.grid`](https://rdrr.io/r/base/expand.grid.html)
  order. The (total) number of levels must match the number of
  dimensions. If `mult.name` is specified, this argument is ignored.

- options:

  If non-`NULL`, a named `list` of arguments to pass to
  [`update.emmGrid`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md),
  just after the object is constructed.

- data:

  A `data.frame` to use to obtain information about the predictors (e.g.
  factor levels). If missing, then
  [`recover_data`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
  is used to attempt to reconstruct the data. See the note with
  [`recover_data`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
  for an important precaution.

- df:

  Numeric value. This is equivalent to specifying `options(df = df)`.
  See
  [`update.emmGrid`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md).

- type:

  Character value. If provided, this is saved as the `"predict.type"`
  setting. See
  [`update.emmGrid`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md)
  and the section below on prediction types and transformations.

- regrid:

  Character, logical, or list. If non-missing, the reference grid is
  reconstructed via
  [`regrid`](https://rvlenth.github.io/emmeans/reference/regrid.md) with
  the argument `transform = regrid`. See the section below on prediction
  types and transformations. *Note:* This argument was named `transform`
  in version 1.7.2 and earlier. For compatibility with old code,
  `transform` is still accepted if found among `...`, as long as it
  doesn't match `tran`.

- nesting:

  If the model has nested fixed effects, this may be specified here via
  a character vector or named `list` specifying the nesting structure.
  Specifying `nesting` overrides any nesting structure that is
  automatically detected. See the section below on Recovering or
  Overriding Model Information.

- offset:

  Numeric scalar value (if a vector, only the first element is used).
  This may be used to add an offset, or override offsets based on the
  model. A common usage would be to specify `offset = 0` for a Poisson
  regression model, so that predictions from the reference grid become
  rates relative to the offset that had been specified in the model.

- sigma:

  Numeric value to use for subsequent predictions or back-transformation
  bias adjustments. If not specified, we use `sigma(object)`, if
  available, and `NULL` otherwise. Note: This applies only when the
  family is `"gaussian"`; for other families, `sigma` is set to `NA` and
  cannot be overridden.

- counterfactuals:

  `counterfactuals` specifies character names of counterfactual factors.
  If this is non-missing, a reference grid is created consisting of
  combinations of counterfactual levels and the actual levels of those
  same factors. This grid is always converted to the response
  transformation scale and averaged over the actual factor levels. See
  the section below on counterfactuals.

- nuisance, non.nuisance, wt.nuis:

  If `nuisance` is a vector of predictor names, those predictors are
  omitted from the reference grid. Instead, the result will be as if we
  had averaged over the levels of those factors, with either equal or
  proportional weights as specified in `wt.nuis` (see the `weights`
  argument in
  [`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md)).
  The factors in `nuisance` must not interact with other factors, not
  even other nuisance factors. Specifying nuisance factors can save
  considerable storage and computation time, and help avoid exceeding
  the maximum reference-grid size (`get_emm_option("rg.limit")`).
  (*Note:* For certain models where the `emm_basis` method returns a
  re-gridded parameterization, nuisance factors cannot be used, and an
  error is thrown.)

- rg.limit:

  Integer limit on the number of reference-grid rows to allow (checked
  before any multivariate responses are included).

- ...:

  Optional arguments passed to
  [`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md),
  [`emm_basis`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md),
  and
  [`recover_data`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md),
  such as `params`, `vcov.` (see **Covariance matrix** below), or
  options such as `mode` for specific model types (see
  [vignette("models",
  "emmeans")](https://rvlenth.github.io/emmeans/doc/models.md)).

## Value

An object of the S4 class `"emmGrid"` (see
[`emmGrid-class`](https://rvlenth.github.io/emmeans/reference/emmGrid-class.md)).
These objects encapsulate everything needed to do calculations and
inferences for estimated marginal means, and contain nothing that
depends on the model-fitting procedure.

## Details

To users, the `ref_grid` function itself is important because most of
its arguments are in effect arguments of
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md) and
related functions, in that those functions pass their `...` arguments to
`ref_grid`.

The reference grid consists of combinations of independent variables
over which predictions are made. Estimated marginal means are defined as
these predictions, or marginal averages thereof. The grid is determined
by first reconstructing the data used in fitting the model (see
[`recover_data`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)),
or by using the `data.frame` provided in `data`.

By “independent variables,” we mean (in most cases) the results of
[`all.vars()`](https://rdrr.io/r/base/allnames.html) applied to the
fixed-effects part of the right-hand side of the model formula. Any
random effects are excluded. Thus, if the model formula in an
[`lme4::lmer`](https://rdrr.io/pkg/lme4/man/lmer.html) call is
`yield ~ fert + seed*density + log(rain) + (1|block/plot)`, the
independent variables are `fert`, `seed`, `density`, and `rain` (not
`log(rain)`). In multivariate models, the dimension of the multivariate
response is also considered an independent variable.

The default reference grid is determined by the observed levels of any
factors, the ordered unique values of character-valued predictors, and
the results of `cov.reduce` for numeric predictors. These may be
overridden using `at`. See also the section below on
recovering/overriding model information.

## Note

The system default for `cov.keep` causes models containing indicator
variables to be handled differently than in emmeans version 1.4.1 or
earlier. To replicate older analyses, change the default via
`emm_options(cov.keep = character(0))`.

Some earlier versions of emmeans offer a `covnest` argument. This is now
obsolete; if `covnest` is specified, it is harmlessly ignored. Cases
where it was needed are now handled appropriately via the code
associated with `cov.keep`.

## Using `cov.reduce` and `cov.keep`

The `cov.keep` argument was not available in emmeans versions 1.4.1 and
earlier. Any covariates named in this list are treated as if they are
factors: all the unique levels are kept in the reference grid. The user
may also specify an integer value, in which case any covariate having no
more than that number of unique values is implicitly included in
`cov.keep`. The default for `cov.keep` is set and retrieved via the
[`emm_options`](https://rvlenth.github.io/emmeans/reference/emm_options.md)
framework, and the system default is `"2"`, meaning that covariates
having only two unique values are automatically treated as two-level
factors. See also the Note below on backward compatibility.

There is a subtle distinction between including a covariate in
`cov.keep` and specifying its values manually in `at`: Covariates
included in `cov.keep` are treated as factors for purposes of weighting,
while specifying levels in `at` will not include the covariate in
weighting. See the `mtcars.lm` example below for an illustration.

`cov.reduce` may be a function, logical value, formula, or a named list
of these. If a single function, it is applied to each covariate. If
logical and `TRUE`, `mean` is used. If logical and `FALSE`, it is
equivalent to including all covariates in `cov.keep`. Use of
`cov.reduce = FALSE` is inadvisable because it can result in a huge
reference grid; it is far better to use `cov.keep`.

If a formula (which must be two-sided), then a model is fitted to that
formula using [`lm`](https://rdrr.io/r/stats/lm.html); then in the
reference grid, its response variable is set to the results of
[`predict`](https://rdrr.io/r/stats/predict.html) for that model, with
the reference grid as `newdata`. (This is done *after* the reference
grid is determined.) A formula is appropriate here when you think
experimental conditions affect the covariate as well as the response.

To allow for situations where a simple
[`lm()`](https://rdrr.io/r/stats/lm.html) call as described above won't
be adequate, a formula of the form `ext ~ fcnname` is also supported,
where the left-hand side may be `ext`, `extern`, or `external` (and must
*not* be a predictor name) and the right-hand side is the name of an
existing function. The function is called with one argument, a data
frame with columns for each variable in the reference grid. The function
is expected to use that frame as new data to be used to obtain
predictions for one or more models; and it should return a named list or
data frame with replacement values for one or more of the covariates.

If `cov.reduce` is a named list, then the above criteria are used to
determine what to do with covariates named in the list. (However,
formula elements do not need to be named, as those names are determined
from the formulas' left-hand sides.) Any unresolved covariates are
reduced using `"mean"`.

Any `cov.reduce` of `cov.keep` specification for a covariate also named
in `at` is ignored.

## Interdependent covariates

Care must be taken when covariate values depend on one another. For
example, when a polynomial model was fitted using predictors `x`, `x2`
(equal to `x^2`), and `x3` (equal to `x^3`), the reference grid will by
default set `x2` and `x3` to their means, which is inconsistent. The
user should instead use the `at` argument to set these to the square and
cube of `mean(x)`. Better yet, fit the model using a formula involving
`poly(x, 3)` or `I(x^2)` and `I(x^3)`; then there is only `x` appearing
as a covariate; it will be set to its mean, and the model matrix will
have the correct corresponding quadratic and cubic terms.

## Matrix covariates

Support for covariates that appear in the dataset as matrices is very
limited. If the matrix has but one column, it is treated like an
ordinary covariate. Otherwise, with more than one column, each column is
reduced to a single reference value – the result of applying
`cov.reduce` to each column (averaged together if that produces more
than one value); you may not specify values in `at`; and they are not
treated as variables in the reference grid, except for purposes of
obtaining predictions.

## Recovering or overriding model information

Ability to support a particular class of `object` depends on the
existence of `recover_data` and `emm_basis` methods – see
[extending-emmeans](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
for details. The call `methods("recover_data")` will help identify
these.

**Data.** In certain models, (e.g., results of
[`glmer.nb`](https://rdrr.io/pkg/lme4/man/glmer.nb.html)), it is not
possible to identify the original dataset. In such cases, we can work
around this by setting `data` equal to the dataset used in fitting the
model, or a suitable subset. Only the complete cases in `data` are used,
so it may be necessary to exclude some unused variables. Using `data`
can also help save computing, especially when the dataset is large. In
any case, `data` must represent all factor levels used in fitting the
model. It *cannot* be used as an alternative to `at`. (Note: If there is
a pattern of `NAs` that caused one or more factor levels to be excluded
when fitting the model, then `data` should also exclude those levels.)

**Covariance matrix.** By default, the variance-covariance matrix for
the fixed effects is obtained from `object`, usually via its
[`vcov`](https://rdrr.io/r/stats/vcov.html) method. However, the user
may override this via a `vcov.` argument, specifying a matrix or a
function. If a matrix, it must be square and of the same dimension and
parameter order of the fixed effects. If a function, must return a
suitable matrix when it is called with arguments `(object, ...)`. Be
careful with possible unintended conflicts with arguments in `...`; for
example,
[`sandwich::vcovHAC()`](https://sandwich.R-Forge.R-project.org/reference/vcovHAC.html)
has optional arguments `adjust` and `weights` that may be intended for
[`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.md)
but will also be passed to `vcov.()`.

**Nested factors.** Having a nesting structure affects marginal
averaging in `emmeans` in that it is done separately for each level (or
combination thereof) of the grouping factors. `ref_grid` tries to
discern which factors are nested in other factors, but it is not always
obvious, and if it misses some, the user must specify this structure via
`nesting`; or later using
[`update.emmGrid`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md).
The `nesting` argument may be a character vector, a named `list`, or
`NULL`. If a `list`, each name should be the name of a single factor in
the grid, and its entry a character vector of the name(s) of its
grouping factor(s). `nested` may also be a character value of the form
`"factor1 %in% (factor2*factor3)"` (the parentheses are optional). If
there is more than one such specification, they may be appended
separated by commas, or as separate elements of a character vector. For
example, these specifications are equivalent:
`nesting = list(state = "country", city = c("state", "country")`,
`nesting = "state %in% country, city %in% (state*country)"`, and
`nesting = c("state %in% country", "city %in% state*country")`.

## Predictors with subscripts and data-set references

When the fitted model contains subscripts or explicit references to data
sets, the reference grid may optionally be post-processed to simplify
the variable names, depending on the `simplify.names` option (see
[`emm_options`](https://rvlenth.github.io/emmeans/reference/emm_options.md)),
which by default is `TRUE`. For example, if the model formula is
`data1$resp ~ data1$trt + data2[[3]] + data2[["cov"]]`, the simplified
predictor names (for use, e.g., in the `specs` for
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md))
will be `trt`, `data2[[3]]`, and `cov`. Numerical subscripts are not
simplified; nor are variables having simplified names that coincide,
such as if `data2$trt` were also in the model.

Please note that this simplification is performed *after* the reference
grid is constructed. Thus, non-simplified names must be used in the `at`
argument (e.g., `` at = list(`data2["cov"]` = 2:4) ``.

If you don't want names simplified, use
`emm_options(simplify.names = FALSE)`.

## Prediction types and transformations

Transformations can exist because of a link function in a generalized
linear model, or as a response transformation, or even both. In many
cases, they are auto-detected, for example a model formula of the form
`sqrt(y) ~ ...`. Even transformations containing multiplicative or
additive constants, such as `2*sqrt(y + pi) ~ ...`, are auto-detected. A
response transformation of `y + 1 ~ ...` is *not* auto-detected, but
`I(y + 1) ~ ...` is interpreted as `identity(y + 1) ~ ...`. A warning is
issued if it gets too complicated. Complex transformations like the
Box-Cox transformation are not auto-detected; but see the help page for
[`make.tran`](https://rvlenth.github.io/emmeans/reference/make.tran.md)
for information on some advanced methods.

There is a subtle difference between specifying `type = "response"` and
`regrid = "response"`. While the summary statistics for the grid itself
are the same, subsequent use in
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md) will
yield different results if there is a response transformation or link
function. With `type = "response"`, EMMs are computed by averaging
together predictions on the *linear-predictor* scale and then
back-transforming to the response scale; while with
`regrid = "response"`, the predictions are already on the response scale
so that the EMMs will be the arithmetic means of those response-scale
predictions. To add further to the possibilities, *geometric* means of
the response-scale predictions are obtainable via
`regrid = "log", type = "response"`. See also the help page for
[`regrid`](https://rvlenth.github.io/emmeans/reference/regrid.md).

*Order-of-processing issues:* The `regrid` argument, if present, is
acted on immediately after the reference grid is constructed, while some
of the `...` arguments may be used to update the object at the very end.
Thus, code like `ref_grid(mod, tran = "sqrt", regrid = "response")` will
not work correctly if the intention was to specify the response
transformation, because the re-grid is done *before* it processes
`tran = "sqrt"`. To get the intended result, do
`regrid(ref_grid(mod, tran = "sqrt"), transform = "response")`.

## Counterfactuals

If `counterfactuals` is specified, the rows of the entire dataset become
part of the reference grid, and the other reference levels are confined
to those named in `counterfactuals`. In this type of analysis (called
G-computation), we substitute (or impute) each combination of
counterfactual levels into the entire dataset. Thus, predictions from
this grid are those of each observation under each of the counterfactual
levels. For this to make sense, we require an assumption of
exchangeability of these levels.

This grid is always converted to the response scale, as G-computation on
the linear-predictor scale produces the same results as ordinary
weighted EMMs. If we have counterfactual factors `A, B`, the reference
grid also includes factors `actual_A, actual_B` which are used to track
which observations originally had the `A, B` levels before they were
changed by the counterfactuals code. We average the response-scale
predictions for each combination of actual levels and imputed levels
(and multivariate levels, if any). See additional discussion of how
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md)
handles counterfactuals under that documentation.

Currently, counterfactuals are not supported when the reference grid
requires post-processing (e.g., ordinal models with `mode = "prob"`).
Cases where we have nested factor levels can be complicated if mixed-in
with counterfactuals, and we make no guarantees. Note that past
implementations included arguments `wt.counter` and `avg.counter`, which
are now deprecated and are just ignored if specified.

## Optional side effect

If the `save.ref_grid` option is set to `TRUE` (see
[`emm_options`](https://rvlenth.github.io/emmeans/reference/emm_options.md)),
The most recent result of `ref_grid`, whether called directly or
indirectly via
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md),
[`emtrends`](https://rvlenth.github.io/emmeans/reference/emtrends.md),
or some other function that calls one of these, is saved in the user's
environment as `.Last.ref_grid`. This facilitates checking what
reference grid was used, or reusing the same reference grid for further
calculations. This automatic saving is disabled by default, but may be
enabled via `emm_options(save.ref_grid = TRUE)`.

## See also

Reference grids are of class
[`emmGrid`](https://rvlenth.github.io/emmeans/reference/emmGrid-class.md),
and several methods exist for them – for example
[`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md).
Reference grids are fundamental to
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md).
Supported models are detailed in
[[`vignette("models", "emmeans")`](https://rvlenth.github.io/emmeans/articles/models.md)](https://rvlenth.github.io/emmeans/doc/models.md).
See
[`update.emmGrid`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md)
for details of arguments that can be in `options` (or in `...`).

## Examples

``` r
fiber.lm <- lm(strength ~ machine*diameter, data = fiber)
ref_grid(fiber.lm)
#>  machine diameter prediction    SE df
#>  A           24.1       40.2 0.777  9
#>  B           24.1       41.6 0.858  9
#>  C           24.1       38.5 0.966  9
#> 

ref_grid(fiber.lm, at = list(diameter = c(15, 25)))
#>  machine diameter prediction    SE df
#>  A             15       30.1 2.110  9
#>  B             15       33.8 2.570  9
#>  C             15       30.6 1.490  9
#>  A             25       41.2 0.750  9
#>  B             25       42.3 0.782  9
#>  C             25       39.3 1.090  9
#> 

if (FALSE) { # \dontrun{
# We could substitute the sandwich estimator vcovHAC(fiber.lm)
# as follows:
summary(ref_grid(fiber.lm, vcov. = sandwich::vcovHAC))
} # }

# If we thought that the machines affect the diameters
# (admittedly not plausible in this example), then we should use:
ref_grid(fiber.lm, cov.reduce = diameter ~ machine)
#>  machine diameter prediction    SE df
#>  A           25.2       41.4 0.749  9
#>  B           26.0       43.2 0.749  9
#>  C           21.2       36.0 0.749  9
#> 

### Model with indicator variables as predictors:
mtcars.lm <- lm(mpg ~ disp + wt + vs * am, data = mtcars)
(rg.default <- ref_grid(mtcars.lm))
#>  disp   wt vs am prediction   SE df
#>   231 3.22  0  0       19.1 1.26 26
#>   231 3.22  1  0       20.0 1.18 26
#>   231 3.22  0  1       18.4 1.14 26
#>   231 3.22  1  1       23.3 1.54 26
#> 
(rg.nokeep <- ref_grid(mtcars.lm, cov.keep = character(0)))
#>  disp   wt    vs    am prediction    SE df
#>   231 3.22 0.438 0.406       19.9 0.484 26
#> 
(rg.at <- ref_grid(mtcars.lm, at = list(vs = 0:1, am = 0:1)))
#>  disp   wt vs am prediction   SE df
#>   231 3.22  0  0       19.1 1.26 26
#>   231 3.22  1  0       20.0 1.18 26
#>   231 3.22  0  1       18.4 1.14 26
#>   231 3.22  1  1       23.3 1.54 26
#> 

# Two of these have the same grid but different weights:
rg.default@grid
#>       disp      wt vs am .wgt.
#> 1 230.7219 3.21725  0  0    12
#> 2 230.7219 3.21725  1  0     7
#> 3 230.7219 3.21725  0  1     6
#> 4 230.7219 3.21725  1  1     7
rg.at@grid
#>       disp      wt vs am .wgt.
#> 1 230.7219 3.21725  0  0     1
#> 2 230.7219 3.21725  1  0     1
#> 3 230.7219 3.21725  0  1     1
#> 4 230.7219 3.21725  1  1     1

### Using cov.reduce formulas...
# Above suggests we can vary disp indep. of other factors - unrealistic
rg.alt <- ref_grid(mtcars.lm, at = list(wt = c(2.5, 3, 3.5)),
    cov.reduce = disp ~ vs * wt)
rg.alt@grid
#>        disp  wt vs am .wgt.
#> 1  185.6376 2.5  0  0    12
#> 2  236.7553 3.0  0  0    12
#> 3  287.8730 3.5  0  0    12
#> 4  125.1602 2.5  1  0     7
#> 5  157.9451 3.0  1  0     7
#> 6  190.7300 3.5  1  0     7
#> 7  185.6376 2.5  0  1     6
#> 8  236.7553 3.0  0  1     6
#> 9  287.8730 3.5  0  1     6
#> 10 125.1602 2.5  1  1     7
#> 11 157.9451 3.0  1  1     7
#> 12 190.7300 3.5  1  1     7

# Alternative to above where we model sqrt(disp)
disp.mod <- lm(sqrt(disp) ~ vs * wt, data = mtcars)
disp.fun <- function(dat)
    list(disp = predict(disp.mod, newdata = dat)^2)
rg.alt2 <- ref_grid(mtcars.lm, at = list(wt = c(2.5, 3, 3.5)),
    cov.reduce = external ~ disp.fun)
#> Error in get(as.character(dep.x[[xnm]][[3]]), inherits = TRUE): object 'disp.fun' not found
rg.alt2@grid
#> Error: object 'rg.alt2' not found


# Multivariate example
MOats.lm = lm(yield ~ Block + Variety, data = MOats)
ref_grid(MOats.lm, mult.names = "nitro")
#>  Block Variety     nitro prediction    SE df
#>  VI    Golden Rain 0           80.3  9.05 10
#>  V     Golden Rain 0           68.9  9.05 10
#>  III   Golden Rain 0           72.9  9.05 10
#>  IV    Golden Rain 0           69.9  9.05 10
#>  II    Golden Rain 0           76.3  9.05 10
#>  I     Golden Rain 0          111.6  9.05 10
#>  VI    Marvellous  0           86.9  9.05 10
#>  V     Marvellous  0           75.6  9.05 10
#>  III   Marvellous  0           79.6  9.05 10
#>  IV    Marvellous  0           76.6  9.05 10
#>  II    Marvellous  0           82.9  9.05 10
#>  I     Marvellous  0          118.3  9.05 10
#>  VI    Victory     0           71.8  9.05 10
#>  V     Victory     0           60.4  9.05 10
#>  III   Victory     0           64.4  9.05 10
#>  IV    Victory     0           61.4  9.05 10
#>  II    Victory     0           67.8  9.05 10
#>  I     Victory     0          103.1  9.05 10
#>  VI    Golden Rain 0.2         84.6 10.80 10
#>  V     Golden Rain 0.2         80.3 10.80 10
#>  III   Golden Rain 0.2         97.9 10.80 10
#>  IV    Golden Rain 0.2         93.3 10.80 10
#>  II    Golden Rain 0.2        107.3 10.80 10
#>  I     Golden Rain 0.2        127.6 10.80 10
#>  VI    Marvellous  0.2         94.6 10.80 10
#>  V     Marvellous  0.2         90.3 10.80 10
#>  III   Marvellous  0.2        107.9 10.80 10
#>  IV    Marvellous  0.2        103.3 10.80 10
#>  II    Marvellous  0.2        117.3 10.80 10
#>  I     Marvellous  0.2        137.6 10.80 10
#>  VI    Victory     0.2         75.8 10.80 10
#>  V     Victory     0.2         71.4 10.80 10
#>  III   Victory     0.2         89.1 10.80 10
#>  IV    Victory     0.2         84.4 10.80 10
#>  II    Victory     0.2         98.4 10.80 10
#>  I     Victory     0.2        118.8 10.80 10
#>  VI    Golden Rain 0.4        108.1 14.20 10
#>  V     Golden Rain 0.4        101.4 14.20 10
#>  III   Golden Rain 0.4        111.4 14.20 10
#>  IV    Golden Rain 0.4        106.1 14.20 10
#>  II    Golden Rain 0.4        115.1 14.20 10
#>  I     Golden Rain 0.4        145.8 14.20 10
#>  VI    Marvellous  0.4        110.6 14.20 10
#>  V     Marvellous  0.4        103.9 14.20 10
#>  III   Marvellous  0.4        113.9 14.20 10
#>  IV    Marvellous  0.4        108.6 14.20 10
#>  II    Marvellous  0.4        117.6 14.20 10
#>  I     Marvellous  0.4        148.3 14.20 10
#>  VI    Victory     0.4        104.3 14.20 10
#>  V     Victory     0.4         97.6 14.20 10
#>  III   Victory     0.4        107.6 14.20 10
#>  IV    Victory     0.4        102.3 14.20 10
#>  II    Victory     0.4        111.3 14.20 10
#>  I     Victory     0.4        141.9 14.20 10
#>  VI    Golden Rain 0.6        114.1 11.90 10
#>  V     Golden Rain 0.6        115.1 11.90 10
#>  III   Golden Rain 0.6        103.4 11.90 10
#>  IV    Golden Rain 0.6        125.4 11.90 10
#>  II    Golden Rain 0.6        132.4 11.90 10
#>  I     Golden Rain 0.6        158.4 11.90 10
#>  VI    Marvellous  0.6        116.1 11.90 10
#>  V     Marvellous  0.6        117.1 11.90 10
#>  III   Marvellous  0.6        105.4 11.90 10
#>  IV    Marvellous  0.6        127.4 11.90 10
#>  II    Marvellous  0.6        134.4 11.90 10
#>  I     Marvellous  0.6        160.4 11.90 10
#>  VI    Victory     0.6        107.8 11.90 10
#>  V     Victory     0.6        108.8 11.90 10
#>  III   Victory     0.6         97.1 11.90 10
#>  IV    Victory     0.6        119.1 11.90 10
#>  II    Victory     0.6        126.1 11.90 10
#>  I     Victory     0.6        152.1 11.90 10
#> 
# Silly illustration of how to use 'mult.levs' to make comb's of two factors
ref_grid(MOats.lm, mult.levs = list(T=LETTERS[1:2], U=letters[1:2]))
#>  Block Variety     T U prediction    SE df
#>  VI    Golden Rain A a       80.3  9.05 10
#>  V     Golden Rain A a       68.9  9.05 10
#>  III   Golden Rain A a       72.9  9.05 10
#>  IV    Golden Rain A a       69.9  9.05 10
#>  II    Golden Rain A a       76.3  9.05 10
#>  I     Golden Rain A a      111.6  9.05 10
#>  VI    Marvellous  A a       86.9  9.05 10
#>  V     Marvellous  A a       75.6  9.05 10
#>  III   Marvellous  A a       79.6  9.05 10
#>  IV    Marvellous  A a       76.6  9.05 10
#>  II    Marvellous  A a       82.9  9.05 10
#>  I     Marvellous  A a      118.3  9.05 10
#>  VI    Victory     A a       71.8  9.05 10
#>  V     Victory     A a       60.4  9.05 10
#>  III   Victory     A a       64.4  9.05 10
#>  IV    Victory     A a       61.4  9.05 10
#>  II    Victory     A a       67.8  9.05 10
#>  I     Victory     A a      103.1  9.05 10
#>  VI    Golden Rain B a       84.6 10.80 10
#>  V     Golden Rain B a       80.3 10.80 10
#>  III   Golden Rain B a       97.9 10.80 10
#>  IV    Golden Rain B a       93.3 10.80 10
#>  II    Golden Rain B a      107.3 10.80 10
#>  I     Golden Rain B a      127.6 10.80 10
#>  VI    Marvellous  B a       94.6 10.80 10
#>  V     Marvellous  B a       90.3 10.80 10
#>  III   Marvellous  B a      107.9 10.80 10
#>  IV    Marvellous  B a      103.3 10.80 10
#>  II    Marvellous  B a      117.3 10.80 10
#>  I     Marvellous  B a      137.6 10.80 10
#>  VI    Victory     B a       75.8 10.80 10
#>  V     Victory     B a       71.4 10.80 10
#>  III   Victory     B a       89.1 10.80 10
#>  IV    Victory     B a       84.4 10.80 10
#>  II    Victory     B a       98.4 10.80 10
#>  I     Victory     B a      118.8 10.80 10
#>  VI    Golden Rain A b      108.1 14.20 10
#>  V     Golden Rain A b      101.4 14.20 10
#>  III   Golden Rain A b      111.4 14.20 10
#>  IV    Golden Rain A b      106.1 14.20 10
#>  II    Golden Rain A b      115.1 14.20 10
#>  I     Golden Rain A b      145.8 14.20 10
#>  VI    Marvellous  A b      110.6 14.20 10
#>  V     Marvellous  A b      103.9 14.20 10
#>  III   Marvellous  A b      113.9 14.20 10
#>  IV    Marvellous  A b      108.6 14.20 10
#>  II    Marvellous  A b      117.6 14.20 10
#>  I     Marvellous  A b      148.3 14.20 10
#>  VI    Victory     A b      104.3 14.20 10
#>  V     Victory     A b       97.6 14.20 10
#>  III   Victory     A b      107.6 14.20 10
#>  IV    Victory     A b      102.3 14.20 10
#>  II    Victory     A b      111.3 14.20 10
#>  I     Victory     A b      141.9 14.20 10
#>  VI    Golden Rain B b      114.1 11.90 10
#>  V     Golden Rain B b      115.1 11.90 10
#>  III   Golden Rain B b      103.4 11.90 10
#>  IV    Golden Rain B b      125.4 11.90 10
#>  II    Golden Rain B b      132.4 11.90 10
#>  I     Golden Rain B b      158.4 11.90 10
#>  VI    Marvellous  B b      116.1 11.90 10
#>  V     Marvellous  B b      117.1 11.90 10
#>  III   Marvellous  B b      105.4 11.90 10
#>  IV    Marvellous  B b      127.4 11.90 10
#>  II    Marvellous  B b      134.4 11.90 10
#>  I     Marvellous  B b      160.4 11.90 10
#>  VI    Victory     B b      107.8 11.90 10
#>  V     Victory     B b      108.8 11.90 10
#>  III   Victory     B b       97.1 11.90 10
#>  IV    Victory     B b      119.1 11.90 10
#>  II    Victory     B b      126.1 11.90 10
#>  I     Victory     B b      152.1 11.90 10
#> 

# Comparing estimates with and without counterfactuals
neuralgia.glm <- glm(Pain ~ Treatment + Sex + Age + Duration, 
                     family = binomial(), data = neuralgia)
emmeans(neuralgia.glm, "Treatment", type = "response")
#>  Treatment  prob     SE  df asymp.LCL asymp.UCL
#>  A         0.196 0.1050 Inf    0.0617     0.475
#>  B         0.126 0.0822 Inf    0.0323     0.384
#>  P         0.855 0.0852 Inf    0.6053     0.958
#> 
#> Results are averaged over the levels of: Sex 
#> Confidence level used: 0.95 
#> Intervals are back-transformed from the logit scale 

emmeans(neuralgia.glm, "Treatment", counterfactuals = "Treatment")
#>  Treatment  prob     SE  df asymp.LCL asymp.UCL
#>  A         0.283 0.0811 Inf    0.1243     0.442
#>  B         0.221 0.0710 Inf    0.0813     0.360
#>  P         0.754 0.0814 Inf    0.5944     0.914
#> 
#> Results are averaged over the levels of: actual_Treatment 
#> Confidence level used: 0.95 


# Using 'params'
require("splines")
#> Loading required package: splines
my.knots = c(2.5, 3, 3.5)
mod = lm(Sepal.Length ~ Species * ns(Sepal.Width, knots = my.knots), data = iris)
## my.knots is not a predictor, so need to name it in 'params'
ref_grid(mod, params = "my.knots") 
#>  Species    Sepal.Width prediction     SE  df
#>  setosa            3.06       4.71 0.1100 135
#>  versicolor        3.06       6.30 0.1070 135
#>  virginica         3.06       6.73 0.0909 135
#> 
```
