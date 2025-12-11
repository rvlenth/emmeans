# Support functions for model extensions

This documents some functions and methods that may be useful to package
developers wishing to add support for emmeans for their model objects.A
user or package developer may add emmeans support for a model class by
writing `recover_data` and `emm_basis` methods for that class. (Users in
need for a quick way to obtain results for a model that is not supported
may be better served by the
[`qdrg`](https://rvlenth.github.io/emmeans/reference/qdrg.md) function.)
There are several other exported functions that may be useful. See the
"xtending" vignette for more details.

## Usage

``` r
recover_data(object, ...)

# S3 method for class 'call'
recover_data(object, trms, na.action, data = NULL,
  params = "pi", frame, pwts, addl.vars, ...)

emm_basis(object, trms, xlev, grid, ...)

.emm_register(classes, pkgname, qdrg = FALSE)

.std.link.labels(fam, misc)

.combine.terms(...)

.aovlist.dffun(k, dfargs)

.cmpMM(X, weights = rep(1, nrow(X)), assign = attr(X$qr, "assign"))

.get.excl(levs, exc, inc)

.get.offset(terms, grid)

.my.vcov(object, vcov. = .statsvcov, ...)

.all.vars(expr, retain = c("\\$", "\\[\\[", "\\]\\]", "'", "\""),
  ...)

.diag(x, nrow, ncol)

.num.key(levs, key)

.emm_vignette(css = system.file("css", "clean-simple.css", package =
  "emmeans"), highlight = NULL, ...)

.hurdle.support(cmu, cshape, cp0, cmean, zmu, zshape, zp0)

.zi.support(zmu, zshape, zp0)
```

## Arguments

- object:

  An object of the same class as is supported by a new method.

- ...:

  Additional parameters that may be supported by the method.

- trms:

  The [`terms`](https://rdrr.io/r/stats/terms.html) component of
  `object` (typically with the response deleted, e.g. via
  [`delete.response`](https://rdrr.io/r/stats/delete.response.html))

- na.action:

  Integer vector of indices of observations to ignore; or `NULL` if none

- data:

  Data frame. Usually, this is `NULL`. However, if non-null, this is
  used in place of the reconstructed dataset. It must have all of the
  predictors used in the model, and any factor levels must match those
  used in fitting the model.

- params:

  Character vector giving the names of any variables in the model
  formula that are *not* predictors. For example, a spline model may
  involve a local variable `knots` that is not a predictor, but its
  value is needed to fit the model. Names of parameters not actually
  used are harmless, and the default value `"pi"` (the only numeric
  constant in base R) is provided in case the model involves it. An
  example involving splines may be found at
  <https://github.com/rvlenth/emmeans/issues/180>.

- frame:

  Optional `data.frame`. Many model objects contain the model frame used
  when fitting the model. In cases where there are no predictor
  transformations, this model frame has all the original predictor
  values and so is usable for recovering the data. Thus, if `frame` is
  non-missing and `data` is `NULL`, a check is made on `trms` and if
  there are no function calls, we use `data = frame`. This can be
  helpful because it provides a modicum of security against the
  possibility that the original data used when fitting the model has
  been altered or removed.

- pwts:

  Optional vector of prior weights. Typically, this may be obtained from
  the fitted `model` via `weights(model)`. If this is provided, it is
  used to set weights as long as it is non-`NULL` and the same length as
  the number of rows of the data.

- addl.vars:

  Character value or vector specifying additional predictors to include
  in the reference grid. These must be names of variables that exist, or
  you will get an error. This may be useful if you need to do additional
  computations later on that depend on these variables; e.g., bias
  adjustments for random slopes of variables not among the fixed
  predictors.

- xlev:

  Named list of factor levels (*excluding* ones coerced to factors in
  the model formula)

- grid:

  A `data.frame` (provided by `ref_grid`) containing the predictor
  settings needed in the reference grid

- classes:

  Character names of one or more classes to be registered. The package
  must contain the functions `recover_data.foo` and `emm_basis.foo` for
  each class `foo` listed in `classes`.

- pkgname:

  Character name of package providing the methods (usually should be the
  second argument of `.onLoad`)

- qdrg:

  Logical value. If `FALSE`, the `recover_data` and `emm_basis` methods
  are registered. If `TRUE`, the `qdrg` method for each class is
  registered instead.

- fam:

  Result of call to `family(object)`

- misc:

  A `list` intended for the `@misc` slot of an `emmGrid` object

- k, dfargs:

  Arguments to `.aovlist.dffun`, which is made available as a
  convenience to developers providing support similar to that provided
  for `aovlist` objects

- X, weights, assign:

  Arguments for `.cmpMM`, which compacts a model matrix `X` into a much
  smaller matrix that has the same row space. Specifically, it returns
  the R portion of its QR decomposition. If `X` is already of class
  `qr`, it is used directly. `weights` should be the weights used in the
  model fit, and `assign` is used for unravelling any pivoting done by
  [`qr`](https://rdrr.io/r/base/qr.html).

- levs, key:

  The `.num.key` function returns the numeric indices of the levels in
  `levs` to the set of all levels in `key`

- exc, inc:

  Arguments for `.get.excl` which is useful in writing `.emmc` functions
  for generating contrast coefficients, and supports arguments `exclude`
  or `include` for excluding or specifying which levels to use.

- terms:

  A `terms` component

- vcov.:

  Function or matrix that returns a suitable covariance matrix. The
  default is `.statsvcov` which is
  [`stats::vcov`](https://rdrr.io/r/stats/vcov.html). The `.my.vcov`
  function should be called in place of
  [`vcov`](https://rdrr.io/r/stats/vcov.html), and it supports the user
  being able to specify a different matrix or function via the optional
  `vcov.` argument.

- expr, retain:

  Arguments for `.all.vars`, which is an alternative to
  [`all.vars`](https://rdrr.io/r/base/allnames.html) that has special
  provisions for retaining the special characters in `retain`, thus
  allowing model specifications like `y ~ data$trt * df[["dose"]]`

- x, nrow, ncol:

  Arguments for `.diag`, which is an alternative to
  [`diag`](https://rdrr.io/r/base/diag.html) that lacks its idiosyncrasy
  of returning an identity matrix when `x` is of length 1.

- css, package, highlight:

  Arguments for `.emm_vignette`, which is a clean and simple alternative
  to such as `html_document` for use as the output style of a Markdown
  file. All the vignettes in the emmeans package use this output style.

- cmu, zmu:

  In `.hurdle.support` and `.zi.support`, these specify a vector of
  back-transformed estimates for the count and zero model, respectively

- cshape, zshape:

  Shape parameter for the count and zero model, respectively

- cp0, zp0:

  Function of `(mu, shape)` for computing Prob(Y = 0) for the count and
  zero model, respectively

- cmean:

  Function of `(mu, shape)` for computing the mean of the count model.
  Typically, this just returns `mu`

## Value

The `recover_data` method must return a
[`data.frame`](https://rdrr.io/r/base/data.frame.html) containing all
the variables that appear as predictors in the model, and attributes
`"call"`, `"terms"`, `"predictors"`, and `"responses"`.
(`recover_data.call` will provide these attributes.)

The `emm_basis` method should return a `list` with the following
elements:

- X:

  The matrix of linear functions over `grid`, having the same number of
  rows as `grid` and the number of columns equal to the length of
  `bhat`.

- bhat:

  The vector of regression coefficients for fixed effects. This should
  *include* any `NA`s that result from rank deficiencies.

- nbasis:

  A matrix whose columns form a basis for non-estimable functions of
  beta, or a 1x1 matrix of `NA` if there is no rank deficiency.

- V:

  The estimated covariance matrix of `bhat`.

- dffun:

  A function of `(k, dfargs)` that returns the degrees of freedom
  associated with `sum(k * bhat)`.

- dfargs:

  A `list` containing additional arguments needed for `dffun`

`.std.link.llabels` returns a modified version of `misc` with the
appropriate information included corresponding to the information in
`fam`

`combine.terms` returns a `terms` object resulting from combining all
the terms or formulas in `...`.

`.get.offset` returns the values, based on `grid`, of any `offset`
component in `terms`

`.hurdle.support` returns a matrix with 3 rows containing the estimated
mean responses and the differentials wrt `cmu` and `zmu`, resp.

`.zi.support` returns a matrix with 2 rows containing the estimated
probabilities of 0 and the differentials wrt `mu`. See the section on
hurdle and zero-inflated models.

## Note

Without an explicit `data` argument, `recover_data` returns the *current
version* of the dataset. If the dataset has changed since the model was
fitted, then this will not be the data used to fit the model. It is
especially important to know this in simulation studies where the data
are randomly generated or permuted, and in cases where several datasets
are processed in one step (e.g., using `dplyr`). In those cases, users
should be careful to provide the actual data used to fit the model in
the `data` argument.

## Details

To create a reference grid, the `ref_grid` function needs to reconstruct
the data used in fitting the model, and then obtain a matrix of linear
functions of the regression coefficients for a given grid of predictor
values. These tasks are performed by calls to `recover_data` and
`emm_basis` respectively. A vignette giving details and examples is
available via [vignette("xtending",
"emmeans")](https://rvlenth.github.io/emmeans/doc/xtending.md)

To extend emmeans's support to additional model types, one need only
write S3 methods for these two functions. The existing methods serve as
helpful guidance for writing new ones. Most of the work for
`recover_data` can be done by its method for class `"call"`, providing
the `terms` component and `na.action` data as additional arguments.
Writing an `emm_basis` method is more involved, but the existing methods
(e.g., `emmeans:::emm_basis.lm`) can serve as models. Certain
`recover_data` and `emm_basis` methods are exported from emmeans. (To
find out, do `methods("recover_data")`.) If your object is based on
another model-fitting object, it may be that all that is needed is to
call one of these exported methods and perhaps make modifications to the
results. Contact the developer if you need others of these exported.

If the model has a multivariate response, `bhat` needs to be “flattened”
into a single vector, and `X` and `V` must be constructed consistently.

In models where a non-full-rank result is possible (often, you can tell
by seeing if there is a `singular.ok` argument in the model-fitting
function),
[`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
and its relatives check the estimability of each prediction, using the
[`nonest.basis`](https://rvlenth.github.io/estimability/reference/nonest.basis.html)
function in the estimability package.

The models already supported are detailed in [the "models"
vignette](https://rvlenth.github.io/emmeans/doc/models.md). Some
packages may provide additional emmeans support for its object classes.

## Communication between methods

If the `recover_data` method generates information needed by
`emm_basis`, that information may be incorporated by creating a `"misc"`
attribute in the returned recovered data. That information is then
passed as the `misc` argument when `ref_grid` calls `emm_basis`.

## Optional hooks

Some models may need something other than standard linear estimates and
standard errors. If so, custom functions may be pointed to via the items
`misc$estHook`, `misc$vcovHook` and `misc$postGridHook`. If just the
name of the hook function is provided as a character string, then it is
retrieved using [`get`](https://rdrr.io/r/base/get.html).

The `estHook` function should have arguments `(object, do.se, tol, ...)`
where `object` is the `emmGrid` object, `do.se` is a logical flag for
whether to return the standard error, and `tol` is the tolerance for
assessing estimability. It should return a matrix with 3 columns: the
estimates, standard errors (`NA` when `do.se==FALSE`), and degrees of
freedom (`NA` for asymptotic). The number of rows should equal
`nrow(linfct(object)`. The `vcovHook` function should have arguments
`(object, tol, ...)` as described. It should return the covariance
matrix for the estimates. Finally, `postGridHook`, if present, is called
at the very end of `ref_grid`; it takes one argument, the constructed
`object`, and should return a suitably modified `emmGrid` object.

## Registering S3 methods for a model class

The `.emm_register` function is provided as a convenience to
conditionally register your S3 methods for a model class,
`recover_data.foo` and `emm_basis.foo`, where `foo` is the class name.
Your package should implement an `.onLoad` function and call
`.emm_register` if emmeans is installed. See the example.

## Support for Hurdle and Zero-inflated models

The functions `.hurdle.support` and `.zi.support` help facilitate
calculations needed to estimate the mean response (count model and zero
model combined) of these models. `.hurdle.support` returns a matrix of
three rows. The first is the estimated mean for a hurdle model, and the
2nd and 3rd rows are differentials for the count and zero models, which
needed for delta-method calculations. To use these, regard the `@linfct`
slot as comprising two sets of columns, for the count and zero models
respectively. To do the delta method calculations, multiply the rows of
the count part by its differentials times `link$mu.eta` evcaluated at
that part of the linear predictor. Do the same for the zero part, using
its differentials and `mu.eta`. If the resulting matrix is **A**, then
the covariance of the mean response is **AVA'** where **V**is the `@V`
slot of the object.

The function `zi.support` works the same way, only it is much simpler,
and is used to estimate the probability of 0 and its differential for
either part of a zero-inflated model or hurdle model.

See the code for `emm_basis.zeroinfl` and `emm_basis.hurdle` for how
these are used with models fitted by the pscl package.

## See also

[Vignette on extending
emmeans](https://rvlenth.github.io/emmeans/doc/xtending.md)

## Examples

``` r
if (FALSE) { # \dontrun{
#--- If your package provides recover_data and emm_grid methods for class 'mymod',
#--- put something like this in your package code -- say in zzz.R:
  .onLoad <- function(libname, pkgname) {
    if (requireNamespace("emmeans", quietly = TRUE))
      emmeans::.emm_register("mymod", pkgname)
  }
} # }
```
