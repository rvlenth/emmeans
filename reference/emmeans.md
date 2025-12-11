# Estimated marginal means (Least-squares means)

Compute estimated marginal means (EMMs) for specified factors or factor
combinations in a linear model; and optionally, comparisons or contrasts
among them. EMMs are also known as least-squares means.

## Usage

``` r
emmeans(object, specs, by = NULL, fac.reduce = function(coefs) apply(coefs,
  2, mean), contr, options = get_emm_option("emmeans"), weights, offset, ...,
  tran)
```

## Arguments

- object:

  An object of class `emmGrid`; or a fitted model object that is
  supported, such as the result of a call to `lm` or `lmer`. Many
  fitted-model objects are supported; see
  [[`vignette("models", "emmeans")`](https://rvlenth.github.io/emmeans/articles/models.md)](https://rvlenth.github.io/emmeans/doc/models.md)
  for details.

- specs:

  A `character` vector specifying the names of the predictors over which
  EMMs are desired. `specs` may also be a `formula` or a `list`
  (optionally named) of valid `spec`s. Use of formulas is described in
  the Overview section below. Specifying `.` as the only factor name
  creates a list of specifications for all model terms.

  **Note:** We recommend *against* using two-sided formulas; see the
  note below for `contr`.

- by:

  A character vector specifying the names of predictors to condition on.

- fac.reduce:

  A function that combines the rows of a matrix into a single vector.
  This implements the “marginal averaging” aspect of EMMs. The default
  is the mean of the rows. Typically if it is overridden, it would be
  some kind of weighted mean of the rows. If `fac.reduce` is nonlinear,
  bizarre results are likely, and EMMs will not be interpretable. NOTE:
  If the `weights` argument is non-missing, `fac.reduce` is ignored.

- contr:

  A character value or `list` specifying contrasts to be added. See
  [`contrast`](https://rvlenth.github.io/emmeans/reference/contrast.md).
  **Note:** `contr` is ignored when `specs` is a formula. **Note 2:**:
  We recommend *against* using this argument; obtaining means and
  obtaining contrasts are two different things, and it is best to do
  them in separate steps, using the
  [`contrast`](https://rvlenth.github.io/emmeans/reference/contrast.md)
  function for the contrasts.

- options:

  If non-`NULL`, a named `list` of arguments to pass to
  [`update.emmGrid`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md),
  just after the object is constructed. (Options may also be included in
  `...`; see the ‘options’ section below.)

- weights:

  Character value, numeric vector, or numeric matrix specifying weights
  to use in averaging predictions. See “Weights” section below. Also, if
  `object` is not already a reference grid, `weights` (if it is
  character) is passed to `ref_grid` as `wt.nuis` in case nuisance
  factors are specified. We can override this by specifying `wt.nuis`
  explicitly. This more-or-less makes the weighting of nuisance factors
  consistent with that of primary factors.

- offset:

  Numeric vector or scalar. If specified, this adds an offset to the
  predictions, or overrides any offset in the model or its reference
  grid. If a vector of length differing from the number of rows in the
  result, it is subsetted or cyclically recycled.

- ...:

  When `object` is not already a `"emmGrid"` object, these arguments are
  passed to
  [`ref_grid`](https://rvlenth.github.io/emmeans/reference/ref_grid.md).
  Common examples are `at`, `cov.reduce`, `data`, `type`, `regrid`,
  `df`, `nesting`, and `vcov.`. Model-type-specific options (see
  [[`vignette("models", "emmeans")`](https://rvlenth.github.io/emmeans/articles/models.md)](https://rvlenth.github.io/emmeans/doc/models.md)),
  commonly `mode`, may be used here as well. In addition, if the model
  formula contains references to variables that are not predictors, you
  must provide a `params` argument with a list of their names. These
  arguments may also be used in lieu of `options`. See the ‘Options’
  section below.

- tran:

  Placeholder to prevent it from being included in `...`. If
  non-missing, it is added to \`options\`. See the ‘Options’ section.

## Value

When `specs` is a `character` vector or one-sided formula, an object of
class `"emmGrid"`. A number of methods are provided for further
analysis, including
[`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md),
[`confint.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md),
[`test.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md),
[`contrast.emmGrid`](https://rvlenth.github.io/emmeans/reference/contrast.md),
and
[`pairs.emmGrid`](https://rvlenth.github.io/emmeans/reference/contrast.md).
When `specs` is a `list` or a `formula` having a left-hand side, the
return value is an
[`emm_list`](https://rvlenth.github.io/emmeans/reference/emm_list-object.md)
object, which is simply a `list` of `emmGrid` objects.

## Details

Users should also consult the documentation for
[`ref_grid`](https://rvlenth.github.io/emmeans/reference/ref_grid.md),
because many important options for EMMs are implemented there, via the
`...` argument.

## Overview

Estimated marginal means or EMMs (sometimes called least-squares means)
are predictions from a linear model over a *reference grid*; or marginal
averages thereof. The
[`ref_grid`](https://rvlenth.github.io/emmeans/reference/ref_grid.md)
function identifies/creates the reference grid upon which `emmeans` is
based.

For those who prefer the terms “least-squares means” or “predicted
marginal means”, functions `lsmeans` and `pmmeans` are provided as
wrappers. See
[`wrappers`](https://rvlenth.github.io/emmeans/reference/wrappers.md).

If `specs` is a `formula`, it should be of the form `~ specs`,
`~ specs | by`, `contr ~ specs`, or `contr ~ specs | by`. The formula is
parsed and the variables therein are used as the arguments `specs`,
`by`, and `contr` as indicated. The left-hand side is optional (and we
don't recommend it), but if specified it should be the name of a
contrast family (e.g., `pairwise`). Operators like `*` or `:` are needed
in the formula to delineate names, but otherwise are ignored.

We now also allow using `.` in `specs`. If this is done, we run
[`joint_tests`](https://rvlenth.github.io/emmeans/reference/joint_tests.md)
on the side to determine all relevant model terms, then replace `specs`
with a corresponding list of specifications. This is a convenience, but
it can create a sizeable `emm_list` object and it is coded rather
inefficiently. While it is permissible to include contrasts via a
`contr` argument or formula left-hand-side, we recommend instead doing
this in a follow-up test with
[`contrast.emm_list`](https://rvlenth.github.io/emmeans/reference/emm_list-object.md).
*Caution*: In models with nested fixed effects, using `.` creates
results where nested factors interact with nesting factors; in those
cases, any contrasts you specify will go across nests, which is likely
not what is desired.

In the special case where the mean (or weighted mean) of all the
predictions is desired, specify `specs` as `~ 1` or `"1"`.

A number of standard contrast families are provided. They can be
identified as functions having names ending in `.emmc` – see the
documentation for
[`emmc-functions`](https://rvlenth.github.io/emmeans/reference/emmc-functions.md)
for details – including how to write your own `.emmc` function for
custom contrasts.

## Weights

If `weights` is a vector, its length must equal the number of
predictions to be averaged to obtain each EMM. If a matrix, each row of
the matrix is used in turn, wrapping back to the first row as needed.
When in doubt about what is being averaged (or how many), first call
`emmeans` with `weights = "show.levels"`.

If `weights` is a string, it should partially match one of the
following:

- `"equal"`:

  Use an equally weighted average.

- `"proportional"`:

  Weight in proportion to the frequencies (in the original data) of the
  factor combinations that are averaged over.

- `"outer"`:

  Weight in proportion to each individual factor's marginal frequencies.
  Thus, the weights for a combination of factors are the outer product
  of the one-factor margins

- `"cells"`:

  Weight according to the frequencies of the cells being averaged.

- `"flat"`:

  Give equal weight to all cells with data, and ignore empty cells.

- `"show.levels"`:

  This is a convenience feature for understanding what is being averaged
  over. Instead of a table of EMMs, this causes the function to return a
  table showing the levels that are averaged over, in the order that
  they appear.

Outer weights are like the 'expected' counts in a chi-square test of
independence, and will yield the same results as those obtained by
proportional averaging with one factor at a time. All except `"cells"`
uses the same set of weights for each mean. In a model where the
predicted values are the cell means, cell weights will yield the raw
averages of the data for the factors involved. Using `"flat"` is similar
to `"cells"`, except nonempty cells are weighted equally and empty cells
are ignored.

## Counterfactuals

Counterfactual reference grids (see the documentation for
[`ref_grid`](https://rvlenth.github.io/emmeans/reference/ref_grid.md))
contain pairs of imputed and actual factor levels, and are handled in a
special way. For starters, the `weights` argument is ignored and we
always use `"cells"` weights. Our understanding is that if factors
`A, B` are specified as counterfactuals, the marginal means for `A`
should still be the same as if `A` were the only counterfactual.
Accordingly, in computing these marginal means, we exclude all cases
where `B != actual_B`, because if `A` were the only counterfactual, `B`
will stay at its actual level. We also take special pains to "remember"
information about actual and imputed levels of counterfactuals so that
appropriate results are obtained when `emmeans` is applied to a previous
`emmeans` result.

## Offsets

Unlike in `ref_grid`, an offset need not be scalar. If not enough values
are supplied, they are cyclically recycled. For a vector of offsets, it
is important to understand that the ordering of results goes with the
first name in `specs` varying fastest. If there are any `by` factors,
those vary slower than all the primary ones, but the first `by` variable
varies the fastest within that hierarchy. See the examples.

## Options and `...`

Arguments that could go in `options` may instead be included in `...`,
typically, arguments such as `type`, `infer`, etc. that in essence are
passed to
[`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md).
Arguments in both places are overridden by the ones in `...`.

There is a danger that `...` arguments could partially match those used
by both `ref_grid` and `update.emmGrid`, creating a conflict. If these
occur, usually they can be resolved by providing complete (or at least
longer) argument names; or by isolating non-`ref_grid` arguments in
`options`; or by calling `ref_grid` separately and passing the result as
`object`. See a not-run example below.

Also, when `specs` is a two-sided formula, or `contr` is specified,
there is potential confusion concerning which `...` arguments apply to
the means, and which to the contrasts. When such confusion is possible,
we suggest doing things separately (a call to `emmeans` with no
contrasts, followed by a call to
[`contrast`](https://rvlenth.github.io/emmeans/reference/contrast.md)).
We treat `adjust` as a special case: it is applied to the `emmeans`
results *only* if there are no contrasts specified, otherwise it is
passed only to `contrast`.

## See also

[`ref_grid`](https://rvlenth.github.io/emmeans/reference/ref_grid.md),
[`contrast`](https://rvlenth.github.io/emmeans/reference/contrast.md),
[vignette("models",
"emmeans")](https://rvlenth.github.io/emmeans/doc/models.md)

## Examples

``` r
warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
emmeans (warp.lm,  ~ wool | tension)
#> tension = L:
#>  wool emmean   SE df lower.CL upper.CL
#>  A      44.6 3.65 48     37.2     51.9
#>  B      28.2 3.65 48     20.9     35.6
#> 
#> tension = M:
#>  wool emmean   SE df lower.CL upper.CL
#>  A      24.0 3.65 48     16.7     31.3
#>  B      28.8 3.65 48     21.4     36.1
#> 
#> tension = H:
#>  wool emmean   SE df lower.CL upper.CL
#>  A      24.6 3.65 48     17.2     31.9
#>  B      18.8 3.65 48     11.4     26.1
#> 
#> Confidence level used: 0.95 
# or equivalently emmeans(warp.lm, "wool", by = "tension")

# 'adjust' argument ignored in emmeans, passed to contrast part...
emmeans (warp.lm, poly ~ tension | wool, adjust = "sidak")
#> $emmeans
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
#> 
#> $contrasts
#> wool = A:
#>  contrast  estimate   SE df t.ratio p.value
#>  linear      -20.00 5.16 48  -3.878  0.0006
#>  quadratic    21.11 8.93 48   2.363  0.0439
#> 
#> wool = B:
#>  contrast  estimate   SE df t.ratio p.value
#>  linear       -9.44 5.16 48  -1.831  0.1412
#>  quadratic   -10.56 8.93 48  -1.182  0.4272
#> 
#> P value adjustment: sidak method for 2 tests 
#> 

# 'adjust' argument NOT ignored ...
emmeans (warp.lm, ~ tension | wool, adjust = "sidak")
#> wool = A:
#>  tension emmean   SE df lower.CL upper.CL
#>  L         44.6 3.65 48    35.53     53.6
#>  M         24.0 3.65 48    14.98     33.0
#>  H         24.6 3.65 48    15.53     33.6
#> 
#> wool = B:
#>  tension emmean   SE df lower.CL upper.CL
#>  L         28.2 3.65 48    19.20     37.2
#>  M         28.8 3.65 48    19.76     37.8
#>  H         18.8 3.65 48     9.76     27.8
#> 
#> Confidence level used: 0.95 
#> Conf-level adjustment: sidak method for 3 estimates 

# Get all sets of EMMs for this model
( allsets <- emmeans(warp.lm, ".") )
#> NOTE: Results may be misleading due to involvement in interactions
#> NOTE: Results may be misleading due to involvement in interactions
#> $`emmeans of wool`
#>  wool emmean   SE df lower.CL upper.CL
#>  A      31.0 2.11 48     26.8     35.3
#>  B      25.3 2.11 48     21.0     29.5
#> 
#> Results are averaged over the levels of: tension 
#> Confidence level used: 0.95 
#> 
#> $`emmeans of tension`
#>  tension emmean   SE df lower.CL upper.CL
#>  L         36.4 2.58 48     31.2     41.6
#>  M         26.4 2.58 48     21.2     31.6
#>  H         21.7 2.58 48     16.5     26.9
#> 
#> Results are averaged over the levels of: wool 
#> Confidence level used: 0.95 
#> 
#> $`emmeans of wool, tension`
#>  wool tension emmean   SE df lower.CL upper.CL
#>  A    L         44.6 3.65 48     37.2     51.9
#>  B    L         28.2 3.65 48     20.9     35.6
#>  A    M         24.0 3.65 48     16.7     31.3
#>  B    M         28.8 3.65 48     21.4     36.1
#>  A    H         24.6 3.65 48     17.2     31.9
#>  B    H         18.8 3.65 48     11.4     26.1
#> 
#> Confidence level used: 0.95 
#> 
contrast(allsets, "eff")    # all effects
#> $`contrasts of emmeans of wool`
#>  contrast estimate   SE df t.ratio p.value
#>  A effect     2.89 1.49 48   1.940  0.0582
#>  B effect    -2.89 1.49 48  -1.940  0.0582
#> 
#> Results are averaged over the levels of: tension 
#> P value adjustment: fdr method for 2 tests 
#> 
#> $`contrasts of emmeans of tension`
#>  contrast estimate   SE df t.ratio p.value
#>  L effect     8.24 2.11 48   3.914  0.0009
#>  M effect    -1.76 2.11 48  -0.836  0.4075
#>  H effect    -6.48 2.11 48  -3.078  0.0052
#> 
#> Results are averaged over the levels of: wool 
#> P value adjustment: fdr method for 3 tests 
#> 
#> $`contrasts of emmeans of wool, tension`
#>  contrast   estimate   SE df t.ratio p.value
#>  A L effect  16.4074 3.33 48   4.929 <0.0001
#>  B L effect   0.0741 3.33 48   0.022  0.9823
#>  A M effect  -4.1481 3.33 48  -1.246  0.4289
#>  B M effect   0.6296 3.33 48   0.189  0.9823
#>  A H effect  -3.5926 3.33 48  -1.079  0.4289
#>  B H effect  -9.3704 3.33 48  -2.815  0.0212
#> 
#> P value adjustment: fdr method for 6 tests 
#> 


if (FALSE) { # \dontrun{
  ### Offsets: Consider a silly example:
  emmeans(warp.lm, ~ tension | wool, offset = c(17, 23, 47)) @ grid
  # note that offsets are recycled so that each level of tension receives
  # the same offset for each wool.
  # But using the same offsets with ~ wool | tension will probably not
  # be what you want because the ordering of combinations is different.
} # }
```
