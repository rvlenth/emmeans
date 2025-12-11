# Update an `emmGrid` object

Objects of class `emmGrid` contain several settings that affect such
things as what arguments to pass to
[`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md).
The `update` method allows safer management of these settings than by
direct modification of its slots.

## Usage

``` r
# S3 method for class 'emmGrid'
update(object, ..., silent = FALSE)

# S3 method for class 'emmGrid'
levels(x) <- value

# S3 method for class 'summary_emm'
update(object, by.vars, mesg, ...)
```

## Arguments

- object:

  An `emmGrid` object

- ...:

  Options to be set. These must match a list of known options (see
  Details)

- silent:

  Logical value. If `FALSE` (the default), a message is displayed if any
  options are not matched. If `TRUE`, no messages are shown.

- x:

  an `emmGrid` object

- value:

  `list` or replacement levels. See the documentation for
  `update.emmGrid` with the `levels` argument, as well as the section
  below on “Replaciong levels”

- by.vars, mesg:

  Attributes that can be altered in `update.summary_emm`

## Value

an updated `emmGrid` object.

`levels<-` replaces the levels of the object in-place. See the section
on replacing levels for details.

## Note

When it makes sense, an option set by `update` will persist into future
results based on that object. But some options are disabled as well. For
example, a `calc` option will be nulled-out if `contrast` is called,
because it probably will not make sense to do the same calculations on
the contrast results, and in fact the variable(s) needed may not even
still exist. `factor(percent)`.

## Details

The names in `...` are partially matched against those that are valid,
and if a match is found, it adds or replaces the current setting. The
valid names are

- `tran`, `tran2`:

  (`list` or `character`) specifies the transformation which, when
  inverted, determines the results displayed by
  [`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md),
  [`predict.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md),
  or [`emmip`](https://rvlenth.github.io/emmeans/reference/emmip.md)
  when `type="response"`. The value may be the name of a standard
  transformation from
  [`make.link`](https://rdrr.io/r/stats/make.link.html) or additional
  ones supported by name, such as `"log2"`; or, for a custom
  transformation, a `list` containing at least the functions `linkinv`
  (the inverse of the transformation) and `mu.eta` (the derivative
  thereof). The
  [`make.tran`](https://rvlenth.github.io/emmeans/reference/make.tran.md)
  function returns such lists for a number of popular transformations.
  See the help page of
  [`make.tran`](https://rvlenth.github.io/emmeans/reference/make.tran.md)
  for details as well as information on the additional named
  transformations that are supported. `tran2` is just like `tran` except
  it is a second transformation (i.e., a response transformation in a
  generalized linear model).

- `tran.mult`:

  Multiple for `tran`. For example, for the response transformation
  `2*sqrt(y)` (or `sqrt(y) + sqrt(y + 1)`, for that matter), we should
  have `tran = "sqrt"` and `tran.mult = 2`. If absent, a multiple of 1
  is assumed.

- `tran.offset`:

  Additive constant before a transformation is applied. For example, a
  response transformation of `log(y + pi)` has `tran.offset = pi`. If no
  value is present, an offset of 0 is assumed.

- `estName`:

  (`character`) is the column label used for displaying predictions or
  EMMs.

- `inv.lbl`:

  (`character)`) is the column label to use for predictions or EMMs when
  `type="response"`.

- `by.vars`:

  (`character)` vector or `NULL`) the variables used for grouping in the
  summary, and also for defining subfamilies in a call to
  [`contrast`](https://rvlenth.github.io/emmeans/reference/contrast.md).

- `pri.vars`:

  (`character` vector) are the names of the grid variables that are not
  in `by.vars`. Thus, the combinations of their levels are used as
  columns in each table produced by
  [`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md).

- `alpha`:

  (numeric) is the default significance level for tests, in
  [`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
  as well as
  [`plot.emmGrid`](https://rvlenth.github.io/emmeans/reference/plot.md)
  when `CIs = TRUE`. Be cautious that methods that depend on specifying
  `alpha` are prone to abuse. See the discussion in
  [[`vignette("basics", "emmeans")`](https://rvlenth.github.io/emmeans/articles/basics.md)](https://rvlenth.github.io/emmeans/doc/basics.html#pvalues).

- `adjust`:

  (`character)`) is the default for the `adjust` argument in
  [`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md).

- `cross.adjust`:

  (`character)`) is the default for the `cross.adjust` argument in
  [`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
  (used for adjusting between groups).

- `famSize`:

  (integer) is the number of means involved in a family of inferences;
  used in Tukey adjustment

- `infer`:

  (`logical` vector of length 2) is the default value of `infer` in
  [`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md).

- `level`:

  (numeric) is the default confidence level, `level`, in
  [`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md).
  *Note:* You must specify all five letters of ‘level’ to distinguish it
  from the slot name ‘levels’.

- `df`:

  (numeric) overrides the default degrees of freedom with a specified
  single value.

- `calc`:

  (list) additional calculated columns. See
  [`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md).

- `null`:

  (numeric) null hypothesis for `summary` or `test` (taken to be zero if
  missing).

- `side`:

  (numeric or character) `side` specification for for `summary` or
  `test` (taken to be zero if missing).

- `sigma`:

  (numeric) Error SD to use in predictions and for bias-adjusted
  back-transformations

- `delta`:

  (numeric) `delta` specification for `summary` or `test` (taken to be
  zero if missing).

- `predict.type` or `type`:

  (character) sets the default method of displaying predictions in
  [`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md),
  [`predict.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md),
  and [`emmip`](https://rvlenth.github.io/emmeans/reference/emmip.md).
  Valid values are `"link"` (with synonyms `"lp"` and `"linear"`), or
  `"response"`.

- `bias.adjust`, `frequentist`:

  (logical) These are used by `summary` if the value of these arguments
  are not specified.

- `estType`:

  (`character`) is used internally to determine what `adjust` methods
  are appropriate. It should match one of
  `c("prediction", "contrast", "pairs")`. As an example of why this is
  needed, the Tukey adjustment should only be used for pairwise
  comparisons (`estType = "pairs"`); if `estType` is some other string,
  Tukey adjustments are not allowed.

- `avgd.over`:

  (`character)` vector) are the names of the variables whose levels are
  averaged over in obtaining marginal averages of predictions, i.e.,
  estimated marginal means. Changing this might produce a misleading
  printout, but setting it to `character(0)` will suppress the “averaged
  over” message in the summary.

- `initMesg`:

  (`character`) is a string that is added to the beginning of any
  annotations that appear below the
  [`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
  display.

- `methDesc`:

  (`character`) is a string that may be used for creating names for a
  list of `emmGrid` objects.

- `nesting`:

  (Character or named `list`) specifies the nesting structure. See
  “Recovering or overriding model information” in the documentation for
  [`ref_grid`](https://rvlenth.github.io/emmeans/reference/ref_grid.md).
  The current nesting structure is displayed by
  [`str.emmGrid`](https://rvlenth.github.io/emmeans/reference/emmGrid-methods.md).

- `levels`:

  named `list` of new levels for the elements of the current `emmGrid`.
  The list name(s) are used as new variable names, and if needed, the
  list is expanded using `expand.grid`. These results replace current
  variable names and levels. This specification changes the `levels`,
  `grid`, `roles`, and `misc` slots in the updated `emmGrid`, and resets
  `pri.vars`, `by.vars`, `adjust`, `famSize`, and `avgd.over`. In
  addition, if there is nesting of factors, that may be altered; a
  warning is issued if it involves something other than mere name
  changes. *Note:* All six letters of `levels` is needed in order to
  distinguish it from `level`.

- `submodel`:

  `formula` or `character` value specifying a submodel (requires this
  feature being supported by underlying methods for the model class).
  When specified, the `linfct` slot is replaced by its aliases for the
  specified sub-model. Any factors in the sub-model that do not appear
  in the model matrix are ignored, as are any interactions that are not
  in the main model, and any factors associate with multivariate
  responses. The estimates displayed are then computed as if the
  sub-model had been fitted. (However, the standard errors will be based
  on the error variance(s) of the full model.) *Note:* The formula
  should refer only to predictor names, *excluding* any function calls
  (such as `factor` or `poly`) that appear in the original model
  formula. See the example.

  The character values allowed should partially match `"minimal"` or
  `"type2"`. With `"minimal"`, the sub-model is taken to be the one only
  involving the surviving factors in `object` (the ones averaged over
  being omitted). Specifying `"type2"` is the same as `"minimal"` except
  only the highest-order term in the submodel is retained, and all
  effects not containing it are orthogonalized-out. Thus, in a purely
  linear situation such as an `lm` model, the joint test of the modified
  object is in essence a type-2 test as in
  [`car::Anova`](https://rdrr.io/pkg/car/man/Anova.html).

  Please note that it is possible (or even likely) that there will be
  disparity between the `grid` and `linfct` slots when a submodel is
  used. This is because `grid` contains the *claimed* values of the
  predictors and `linfct` contains *aliases* of them computed from the
  submodel.

  For some objects such as generalized linear models, specifying
  `submodel` will typically not produce the same estimates or type-2
  tests as would be obtained by actually fitting a separate model with
  those specifications. The reason is that those models are fitted by
  iterative-reweighting methods, whereas the `submodel` calculations
  preserve the final weights used in fitting the full model.

- (any other slot name):

  If the name matches an element of `slotNames(object)` other than
  `levels`, that slot is replaced by the supplied value, if it is of the
  required class (otherwise an error occurs).

  The user must be very careful in replacing slots because they are
  interrelated; for example, the lengths and dimensions of `grid`,
  `linfct`, `bhat`, and `V` must conform.

## Replacing levels

The `levels<-` method uses `update.emmGrid` to replace the levels of one
or more factors. This method allows selectively replacing the levels of
just one factor (via subsetting operators), whereas
`update(x, levels = list(...))` requires a list of *all* factors and
their levels. If any factors are to be renamed, we must replace all
levels and include the new names in the replacements. See the examples.

## Method for `summary_emm` objects

This method exists so that we can change the way a summary is displayed,
by changing the by variables or the annotations.

## See also

[`emm_options`](https://rvlenth.github.io/emmeans/reference/emm_options.md)

## Examples

``` r
# Using an already-transformed response:
pigs.lm <- lm(log(conc) ~ source * factor(percent), data = pigs)

# Reference grid that knows about the transformation
# and asks to include the sample size in any summaries:
pigs.rg <- update(ref_grid(pigs.lm), tran = "log", 
                    predict.type = "response",
                    calc = c(n = ~.wgt.))
emmeans(pigs.rg, "source")
#> NOTE: Results may be misleading due to involvement in interactions
#>  source response   SE df  n lower.CL upper.CL
#>  fish       29.9 1.12 17 10     27.6     32.3
#>  soy        38.9 1.60 17 10     35.7     42.4
#>  skim       46.1 1.97 17  9     42.1     50.4
#> 
#> Results are averaged over the levels of: percent 
#> Confidence level used: 0.95 
#> Intervals are back-transformed from the log scale 

# Obtain estimates for the additive model
# [Note that the submodel refers to 'percent', not 'factor(percent)']
emmeans(pigs.rg, "source", submodel = ~ source + percent)
#> NOTE: Results may be misleading due to involvement in interactions
#>  source response   SE df  n lower.CL upper.CL
#>  fish       29.8 1.10 17 10     27.6     32.2
#>  soy        39.1 1.48 17 10     36.1     42.4
#>  skim       44.6 1.77 17  9     41.0     48.5
#> 
#> Results are averaged over the levels of: percent 
#> submodel: ~ source + percent 
#> Confidence level used: 0.95 
#> Intervals are back-transformed from the log scale 

# Type II ANOVA
joint_tests(pigs.rg, submodel = "type2")
#>  model term     df1 df2 F.ratio p.value
#>  source           2  17  28.291 <0.0001
#>  percent          3  17   7.827  0.0017
#>  source:percent   6  17   0.926  0.5011
#> 

## Changing levels of one factor
newrg <- pigs.rg
levels(newrg)$source <- 1:3
newrg
#>  source percent response   SE df n
#>       1       9     25.7 2.11 17 2
#>       2       9     34.4 2.31 17 3
#>       3       9     35.2 2.36 17 3
#>       1      12     30.9 2.07 17 3
#>       2      12     39.6 2.66 17 3
#>       3      12     43.2 2.90 17 3
#>       1      15     31.0 2.55 17 2
#>       2      15     39.2 2.63 17 3
#>       3      15     49.6 4.08 17 2
#>       1      18     32.3 2.17 17 3
#>       2      18     42.9 4.99 17 1
#>       3      18     59.8 6.95 17 1
#> 

## Unraveling a previously standardized covariate
zd = scale(fiber$diameter)
fibz.lm <- lm(strength ~ machine * zd, data = fiber)
(fibz.rg <- ref_grid(fibz.lm, at = list(zd = -2:2)))   ### 2*SD range
#>  machine zd prediction    SE df
#>  A       -2       30.7 2.020  9
#>  B       -2       34.2 2.470  9
#>  C       -2       31.1 1.410  9
#>  A       -1       35.4 1.280  9
#>  B       -1       37.9 1.580  9
#>  C       -1       34.8 0.803  9
#>  A        0       40.2 0.777  9
#>  B        0       41.6 0.858  9
#>  C        0       38.5 0.966  9
#>  A        1       45.0 0.979  9
#>  B        1       45.3 0.929  9
#>  C        1       42.3 1.690  9
#>  A        2       49.8 1.650  9
#>  B        2       49.0 1.690  9
#>  C        2       46.0 2.520  9
#> 
lev <- levels(fibz.rg)
levels(fibz.rg) <- list (
    machine = lev$machine,
    diameter = with(attributes(zd), 
                    `scaled:center` + `scaled:scale` * lev$zd) )
fibz.rg
#>  machine diameter prediction    SE df
#>  A           15.5       30.7 2.020  9
#>  B           15.5       34.2 2.470  9
#>  C           15.5       31.1 1.410  9
#>  A           19.8       35.4 1.280  9
#>  B           19.8       37.9 1.580  9
#>  C           19.8       34.8 0.803  9
#>  A           24.1       40.2 0.777  9
#>  B           24.1       41.6 0.858  9
#>  C           24.1       38.5 0.966  9
#>  A           28.5       45.0 0.979  9
#>  B           28.5       45.3 0.929  9
#>  C           28.5       42.3 1.690  9
#>  A           32.8       49.8 1.650  9
#>  B           32.8       49.0 1.690  9
#>  C           32.8       46.0 2.520  9
#> 

### Compactify results with a by variable
update(joint_tests(pigs.rg, by = "source"), by = NULL)
#>  model term source df1 df2 F.ratio p.value
#>  percent    fish     3  17   1.712  0.2023
#>  percent    soy      3  17   1.290  0.3097
#>  percent    skim     3  17   6.676  0.0035
#> 
```
