# For developers: Extending \*\*emmeans\*\*

## Contents

This vignette explains how developers may incorporate **emmeans**
support in their packages. If you are a user looking for a quick way to
obtain results for an unsupported model, you are probably better off
trying to use the
[`qdrg()`](https://rvlenth.github.io/emmeans/reference/qdrg.md)
function.

1.  [Introduction](#intro)
2.  [Data example](#dataex)
3.  [Supporting `rlm` objects](#rlm)
4.  [Quick and dirty support](#qdrg)
5.  [Supporting `lqs` objects](#lqs)
6.  [Communication between methods](#communic)
7.  [Hook functions](#hooks)
8.  [Re-gridded basis](#regridded)
9.  [Exported methods from **emmeans**](#exported)
10. [Existing support for `rsm` objects](#rsm)
11. [Dispatching and restrictions](#dispatch)
12. [Exporting and registering your methods](#exporting)
13. [Conclusions](#concl)

[Index of all vignette
topics](https://rvlenth.github.io/emmeans/articles/vignette-topics.md)

## Introduction

Suppose you want to use **emmeans** for some type of model that it
doesn’t (yet) support. Or, suppose you have developed a new package with
a fancy model-fitting function, and you’d like it to work with
**emmeans**. What can you do? Well, there is hope because **emmeans** is
designed to be extended.

The first thing to do is to look at the help page for extending the
package:

``` r
help("extending-emmeans", package="emmeans")
```

It gives details about the fact that you need to write two S3 methods,
`recover_data` and `emm_basis`, for the class of object that your
model-fitting function returns. The `recover_data` method is needed to
recreate the dataset so that the reference grid can be identified. The
`emm_basis` method then determines the linear functions needed to
evaluate each point in the reference grid and to obtain associated
information—such as the variance-covariance matrix—needed to do
estimation and testing.

These methods must also be exported from your package so that they are
available to users. See the section on [exporting the
methods](#exporting) for details and suggestions.

This vignette presents an example where suitable methods are developed,
and discusses a few issues that arise.

[Back to Contents](#contents)

## Data example

The **MASS** package contains various functions that do robust or
outlier-resistant model fitting. We will cobble together some
**emmeans** support for these. But first, let’s create a suitable
dataset (a simulated two-factor experiment) for testing.

``` r
fake = expand.grid(rep = 1:5, A = c("a1","a2"), B = c("b1","b2","b3"))
fake$y = c(11.46,12.93,11.87,11.01,11.92,17.80,13.41,13.96,14.27,15.82,
           23.14,23.75,-2.09,28.43,23.01,24.11,25.51,24.11,23.95,30.37,
           17.75,18.28,17.82,18.52,16.33,20.58,20.55,20.77,21.21,20.10)
```

The `y` values were generated using predetermined means and
Cauchy-distributed errors. There are some serious outliers in these
data.

## Supporting `rlm`

The **MASS** package provides an `rlm` function that fits
robust-regression models using *M* estimation. We’ll fit a model using
the default settings for all tuning parameters:

``` r
library(MASS)
fake.rlm = rlm(y ~ A * B, data = fake)

library(emmeans)
emmeans(fake.rlm, ~ B | A)
```

``` ro
## A = a1:
##  B  emmean    SE df asymp.LCL asymp.UCL
##  b1   11.8 0.477 NA      10.9      12.8
##  b2   23.3 0.477 NA      22.4      24.2
##  b3   17.8 0.477 NA      16.9      18.7
## 
## A = a2:
##  B  emmean    SE df asymp.LCL asymp.UCL
##  b1   14.7 0.477 NA      13.7      15.6
##  b2   24.7 0.477 NA      23.8      25.6
##  b3   20.6 0.477 NA      19.7      21.6
## 
## Confidence level used: 0.95
```

The first lesson to learn about extending **emmeans** is that sometimes,
it already works! It works here because `rlm` objects inherit from `lm`,
which is supported by the **emmeans** package, and `rlm` objects aren’t
enough different to create any problems.

[Back to Contents](#contents)

## Quick and dirty support

Later, we will talk about how to fully support a model object in
**emmeans**. But it is often very easy to provide *partial* support via
the [`qdrg()`](https://rvlenth.github.io/emmeans/reference/qdrg.md)
(“quick and dirty reference grid”) function.
[`qdrg()`](https://rvlenth.github.io/emmeans/reference/qdrg.md) creates
a reference grid for a model, and that reference grid can subsequently
be used in a call to
[`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.md) or
many other functions (but not
[`emtrends()`](https://rvlenth.github.io/emmeans/reference/emtrends.md)).
In the past,
[`qdrg()`](https://rvlenth.github.io/emmeans/reference/qdrg.md) was just
a standalone function; but now, it can dispath S3 methods based on the
class of its `object` argument. Consider for example a model fitted in
the **robmixglm** package:

``` r
require(robmixglm, quietly = TRUE)
fit <- robmixglm(inverse(conc) ~ source + factor(percent), family = "gamma",
                 data = pigs, cores = 1)
```

Examining the model object `fmx`, it is not at all obvious how to get
the needed parameters for
[`qdrg()`](https://rvlenth.github.io/emmeans/reference/qdrg.md); but
they are buried in there. Here is an S3 method that supports this model
class:

``` r
qdrg.robmixglm <- function(object, data = eval(object$call$data), ...) {
    coef <- coef(object)
    idx <- seq_along(coef)
    qdrg(formula = formula(object), data = data, coef = coef,
         vcov = object$fit@vcov[idx, idx, drop = FALSE], ...)
}
```

In this method, the `vcov` matrix includes covariances for some extra
parameters besides the regression coefficients, so we have to use a
subset of that matrix. Note that
[`qdrg()`](https://rvlenth.github.io/emmeans/reference/qdrg.md)
functions as a standalone function as long as its `object` argument is
missing; and a typical `qdrg` method will just determine what arguments
to pass to this non-generic
[`qdrg()`](https://rvlenth.github.io/emmeans/reference/qdrg.md). Let’s
get the reference grid for a selection of `bp` values:

``` r
rg <- qdrg(object = fit, link = "log")
emmeans(rg, ~ ., type = "response")
```

``` ro
## $`emmeans of source`
##  source response    SE  df asymp.LCL asymp.UCL
##  fish       29.7 0.964 Inf      27.9      31.7
##  soy        39.0 1.290 Inf      36.6      41.7
##  skim       44.1 1.550 Inf      41.1      47.2
## 
## Results are averaged over the levels of: percent 
## Confidence level used: 0.95 
## Intervals are back-transformed from the inverse[log] scale 
## 
## $`emmeans of percent`
##  percent response   SE  df asymp.LCL asymp.UCL
##        9     31.2 1.13 Inf      29.0      33.4
##       12     37.4 1.27 Inf      35.0      40.0
##       15     38.7 1.50 Inf      35.9      41.7
##       18     42.1 1.98 Inf      38.4      46.2
## 
## Results are averaged over the levels of: source 
## Confidence level used: 0.95 
## Intervals are back-transformed from the inverse[log] scale
```

Here, the `link` argument got passed via `...` to
[`qdrg()`](https://rvlenth.github.io/emmeans/reference/qdrg.md); It was
needed since I don’t immediately see how to get the link function from
the object, though I imagine there is a way.

Be careful when using these methods; since `object` is not the first
argument of the generic
[`qdrg()`](https://rvlenth.github.io/emmeans/reference/qdrg.md)
function, it is best to specify `object =` in the call to this method,
though the code tries to work around this.

Package developers may provide minimal emmeans support by providing a
`qdrg` method like this. If you do this in your package, you should
[export and register the method](#exporting), and add
`emmeans (>=2.0.0)` and `estimability` to the package’s `Suggests` list.

[Back to Contents](#contents)

## Supporting `lqs` objects

The **MASS** resistant-regression functions `lqs`, `lmsreg`, and
`ltsreg` are another story, however. They create `lqs` objects that are
not extensions of any other class, and have other issues, including not
even having a `vcov` method. So for these, we really do need to write
new methods for `lqs` objects. First, let’s fit a model.

``` r
fake.lts = ltsreg(y ~ A * B, data = fake)
```

### The `recover_data` method

It is usually an easy matter to write a `recover_data` method. Look at
the one for `lm` objects:

``` r
emmeans:::recover_data.lm
```

``` ro
## function (object, frame = object$model, ...) 
## {
##     fcall = object$call
##     recover_data(fcall, delete.response(terms(object)), object$na.action, 
##         frame = frame, pwts = weights(object), ...)
## }
## <bytecode: 0x563ff86fe9d8>
## <environment: namespace:emmeans>
```

Note that all it does is obtain the `call` component and call the method
for class `call`, with additional arguments for its `terms` component
and `na.action`. It happens that we can access these attributes in
exactly the same way as for `lm` objects; so:

``` r
recover_data.lqs = emmeans:::recover_data.lm
```

Let’s test it:

``` r
rec.fake = recover_data(fake.lts)
head(rec.fake)
```

``` ro
##    A  B
## 1 a1 b1
## 2 a1 b1
## 3 a1 b1
## 4 a1 b1
## 5 a1 b1
## 6 a2 b1
```

Our recovered data excludes the response variable `y` (owing to the
`delete.response` call), and this is fine.

#### Special arguments

By the way, there are two special arguments `data` and `params` that may
be handed to `recover_data` via `ref_grid` or `emmeans` or a related
function; and you may need to provide for if you don’t use the
`recover_data.call` function. The `data` argument is needed to cover a
desperate situation that occurs with certain kinds of models where the
underlying data information is not saved with the object—e.g., models
that are fitted by iteratively modifying the data. In those cases, the
only way to recover the data is to for the user to give it explicitly,
and `recover_data` just adds a few needed attributes to it.

The `params` argument is needed when the model formula refers to
variables besides predictors. For example, a model may include a spline
term, and the knots are saved in the user’s environment as a vector and
referred to in the call to fit the model. In trying to recover the data,
we try to construct a data frame containing all the variables present on
the right-hand side of the model, but if some of those are scalars or of
different lengths than the number of observations, an error occurs. So
you need to exclude any names in `params` when reconstructing the data.

Many model objects contain the model frame as a slot; for example, a
model fitted with `lm(..., model = TRUE)` has a member `$model`
containing the model frame. This can be useful for recovering the data,
provided none of the predictors are transformed (when predictors are
transformed, the original predictor values are not in the model frame so
it’s harder to recover them). Therefore, when the model frame is
available in the model object, it should be provided in the `frame`
argument of
[`recover_data.call()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md);
then when `data = NULL`, a check is made on `trms`, and if it has no
function calls, then `data` is set to `frame`. Of course, in the rarer
case where the original data are available in the model object, specify
that as `data`.

#### Error handling

If you check for any error conditions in `recover_data`, simply have it
return a character string with the desired message, rather than invoking
`stop`. This provides a cleaner exit. The reason is that whenever
`recover_data` throws an error, an informative message suggesting that
`data` or `params` be provided is displayed. But a character return
value is tested for and throws a different error with your string as the
message.

### The `emm_basis` method

The `emm_basis` method has four required arguments:

``` r
args(emmeans:::emm_basis.lm)
```

``` ro
## function (object, trms, xlev, grid, ...) 
## NULL
```

These are, respectively, the model object, its `terms` component (at
least for the right-hand side of the model), a `list` of levels of the
factors, and the grid of predictor combinations that specify the
reference grid.

The function must obtain six things and return them in a named `list`.
They are the matrix `X` of linear functions for each point in the
reference grid, the regression coefficients `bhat`; the
variance-covariance matrix `V`; a matrix `nbasis` for non-estimable
functions; a function `dffun(k,dfargs)` for computing degrees of freedom
for the linear function `sum(k*bhat)`; and a list `dfargs` of arguments
to pass to `dffun`. Optionally, the returned list may include a
`model.matrix` element (the model matrix for the data or a compact
version thereof obtained via
[`.cmpMM()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)),
which, if included, enables the `submodel` option.

To write your own `emm_basis` function, examining some of the existing
methods can help; but the best resource is the `predict` method for the
object in question, looking carefully to see what it does to predict
values for a new set of predictors (e.g., `newdata` in `predict.lm`).
Following this advice, let’s take a look at it:

``` r
MASS:::predict.lqs
```

``` ro
## function (object, newdata, na.action = na.pass, ...) 
## {
##     if (missing(newdata)) 
##         return(fitted(object))
##     Terms <- delete.response(terms(object))
##     m <- model.frame(Terms, newdata, na.action = na.action, xlev = object$xlevels)
##     if (!is.null(cl <- attr(Terms, "dataClasses"))) 
##         .checkMFClasses(cl, m)
##     X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
##     drop(X %*% object$coefficients)
## }
## <bytecode: 0x564002cbd5b0>
## <environment: namespace:MASS>
```

###### 

Based on this, here is a listing of an `emm_basis` method for `lqs`
objects:

``` r
emm_basis.lqs = function(object, trms, xlev, grid, ...) { 
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)
    X = model.matrix(trms, m, contrasts.arg = object$contrasts) 
    bhat = coef(object) 
    Xmat = model.matrix(trms, data=object$model)                      # 5
    V = rev(object$scale)[1]^2 * solve(t(Xmat) %*% Xmat)
    nbasis = matrix(NA) 
    dfargs = list(df = nrow(Xmat) - ncol(Xmat))
    dffun = function(k, dfargs) dfargs$df
    list(X = X, bhat = bhat, nbasis = nbasis, V = V,                  #10
         dffun = dffun, dfargs = dfargs)
}
```

Before explaining it, let’s verify that it works:

``` r
emmeans(fake.lts, ~ B | A)
```

``` ro
## A = a1:
##  B  emmean    SE df lower.CL upper.CL
##  b1   11.9 0.238 24     11.5     12.4
##  b2   23.2 0.238 24     22.7     23.7
##  b3   17.8 0.238 24     17.4     18.3
## 
## A = a2:
##  B  emmean    SE df lower.CL upper.CL
##  b1   14.3 0.238 24     13.8     14.8
##  b2   24.0 0.238 24     23.5     24.5
##  b3   20.8 0.238 24     20.3     21.3
## 
## Confidence level used: 0.95
```

Hooray! Note the results are comparable to those we had for `fake.rlm`,
albeit the standard errors are quite a bit smaller. (In fact, the SEs
could be misleading; a better method for estimating covariances should
probably be implemented, but that is beyond the scope of this vignette.)

[Back to Contents](#contents)

### Dissecting `emm_basis.lqs`

Let’s go through the listing of this method, line-by-line:

- Lines 2–3: Construct the linear functions, `X`. This is a pretty
  standard two-step process: First obtain a model frame, `m`, for the
  grid of predictors, then pass it as data to `model.matrix` to create
  the associated design matrix. As promised, this code is essentially
  identical to what you find in `predict.lqs`.

- Line 4: Obtain the coefficients, `bhat`. Most model objects have a
  `coef` method.

- Lines 5–6: Obtain the covariance matrix, `V`, of `bhat`. In many
  models, this can be obtained using the object’s `vcov` method. But not
  in this case. Instead, I cobbled one together using the inverse of the
  **X’X** matrix as in ordinary regression, and the variance estimate
  found in the last element of the `scale` element of the object. This
  probably under-estimates the variances and distorts the covariances,
  because robust estimators have some efficiency loss.

- Line 7: Compute the basis for non-estimable functions. This applies
  only when there is a possibility of rank deficiency in the model. But
  `lqs` methods don’t allow rank deficiencies, so it we have fitted such
  a model, we can be sure that all linear functions are estimable; we
  signal that by setting `nbasis` equal to a 1 x 1 matrix of `NA`. If
  rank deficiency were possible, the **estimability** package (which is
  required by **emmeans**) provides a `nonest.basis` function that makes
  this fairly painless—I would have coded
  `nbasis = estimability::nonest.basis(Xmat)`.

  There some subtleties you need to know regarding estimability. Suppose
  the model is rank-deficient, so that the design matrix **X** has *p*
  columns but rank *r* \< *p*. In that case, `bhat` should be of length
  *p* (not *r*), and there should be *p* - *r* elements equal to `NA`,
  corresponding to columns of **X** that were excluded from the fit.
  Also, `X` should have all *p* columns. In other words, do not alter or
  throw-out columns of `X` or their corresponding elements of
  `bhat`—even those with `NA` coefficients—as they are essential for
  assessing estimability. `V` should be *r* x *r*, however—the
  covariance matrix for the non-excluded predictors.

- Lines 8–9: Obtain `dffun` and `dfargs`. This is a little awkward
  because it is designed to allow support for mixed models, where
  approximate methods may be used to obtain degrees of freedom. The
  function `dffun` is expected to have two arguments: `k`, the vector of
  coefficients of `bhat`, and `dfargs`, a list containing any additional
  arguments. In this case (and in many other models), the degrees of
  freedom are the same regardless of `k`. We put the required degrees of
  freedom in `dfargs` and write `dffun` so that it simply returns that
  value. (Note: If asymptotic tests and CIs are desired, return `Inf`
  degrees of freedom.)

- Line 10: Return these results in a named list.

[Back to Contents](#contents)

## Communication between methods

If you need to pass information obtained in
[`recover_data()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
to the
[`emm_basis()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
method, simply incorporate it as `attr(data, "misc")` where `data` is
the dataset returned by
[`recover_data()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md).
Subsequently, that attribute is available in `emm_grid()` by adding a
`misc` argument.

## Hook functions

Most linear models supported by **emmeans** have straightforward
structure: Regression coefficients, their covariance matrix, and a set
of linear functions that define the reference grid. However, a few are
more complex. An example is the `clm` class in the **ordinal** package,
which allows a scale model in addition to the location model. When a
scale model is used, the scale parameters are included in the model
matrix, regression coefficients, and covariance matrix, and we can’t
just use the usual matrix operations to obtain estimates and standard
errors. To facilitate using custom routines for these tasks, the
`emm_basis.clm` function function provided in **emmeans** includes, in
its `misc` part, the names (as character constants) of two “hook”
functions: `misc$estHook` has the name of the function to call when
computing estimates, standard errors, and degrees of freedom (for the
`summary` method); and `misc$vcovHook` has the name of the function to
call to obtain the covariance matrix of the grid values (used by the
`vcov` method). These functions are called in lieu of the usual built-in
routines for these purposes, and return the appropriately sized
matrices.

In addition, you may want to apply some form of special post-processing
after the reference grid is constructed. To provide for this, give the
name of your function to post-process the object in `misc$postGridHook`.
Again, `clm` objects (as well as `polr` in the **MASS** package) serve
as an example. They allow a `mode` specification that in two cases,
calls for post-processing. The `"cum.prob"` mode uses the `regrid`
function to transform the linear predictor to the cumulative-probability
scale. And the `"prob"` mode performs this, as well as applying the
contrasts necessary to convert the cumulative probabilities into the
class probabilities.

## Re-gridded basis

Sometimes your `emm_basis` method may essentially create a re-gridded
basis, where `X` and `bhat` are not actually a model matrix and
regression coefficients, but instead, `X` is the identity, `bhat`
comprises the predictions at each grid point, and `V` is the covariance
matrix of those predictions. In those cases, we recommend also setting
`misc$regrid.flag = TRUE`. Currently, this flag is used only for
checking whether the `nuisance` argument can be used in
[`ref_grid()`](https://rvlenth.github.io/emmeans/reference/ref_grid.md),
and it is not absolutely necessary because we also check to see if `X`
is the identity. But it provides a more efficient and reliable check.
The code for nuisance factors relies on the structure of model matrices
where columns are associated with model terms. So it is not possible to
process nuisance factors with a re-gridded basis.

[Back to Contents](#contents)

## Exported methods from **emmeans**

For package developers’ convenience, **emmeans** exports some of its S3
methods for `recover_data` and/or `emm_basis`—use
`methods("recover_data")` and `methods("emm_basis")` to discover which
ones. It may be that all you need is to invoke one of those methods and
perhaps make some small changes—especially if your model-fitting
algorithm makes heavy use of an existing model type supported by
**emmeans**.

A few additional functions are exported because they may be useful to
developers. They are as follows:

- `emmeans::.all.vars(expr, retain)` Some users of your package may
  include `$` or `[[]]` operators in their model formulas. If you need
  to get the variable names,
  [`base::all.vars`](https://rdrr.io/r/base/allnames.html) will probably
  not give you what you need. For example, if
  `form = ~ data$x + data[[5]]`, then `base::all.vars(form)` returns the
  names `"data"` and `"x"`, whereas `emmeans::.all.vars(form)` returns
  the names `"data$x"` and `"data[[5]]"`. The `retain` argument may be
  used to specify regular expressions for patterns to retain as parts of
  variable names.

- `emmeans::.diag(x, nrow, ncol)` The base `diag` function has a booby
  trap whereby, for example, `diag(57.6)` returns a 57 x 57 identity
  matrix rather than a 1 x 1 matrix with 57.6 as its only element. But
  `emmeans::.diag(57.6)` will return the latter. The function works
  identically to `diag` except for its tail run around the
  identity-matrix trap.

- `emmeans::.aovlist.dffun(k, dfargs)` This function is exported because
  it is needed for computing degrees of freedom for models fitted using
  `aov`, but it may be useful for other cases where Satterthwaite
  degrees-of-freedom calculations are needed. It requires the `dfargs`
  slot to contain analogous contents.

- `emmeans::.get.offset(terms, grid)` If `terms` is a model formula
  containing an `offset` call, this is will compute that offset in the
  context of `grid` (a `data.frame`).

- `emmeans::.my.vcov(object, ...)` In a call to `ref_grid`, `emmeans`,
  etc., the user may use `vcov.` to specify an alternative function or
  matrix to use as the covariance matrix of the fixed-effects
  coefficients. This function supports that feature. Calling `.my.vcov`
  in place of the `vcov` method will substitute the user’s `vcov.` when
  it is specified.

- `emmeans::.std.link.labels(fam, misc)` This is useful in `emm_basis`
  methods for generalized linear models. Call it with `fam` equal to the
  `family` object for your model, and `misc` either an existing list, or
  just [`list()`](https://rdrr.io/r/base/list.html) if none. It returns
  a new `misc` list containing the link function and, in some cases,
  extra features that are used for certain types of link functions
  (e.g., for a log link, the setups for returning ratio comparisons with
  `type = "response"`).

- `emmeans::.num.key(levs, key)` Returns integer indices of elements of
  `key` in `levs` when `key` is a character vector; or just returns
  integer values if already integer. Also throws an error if levels are
  mismatched or indices exceed legal range. This is useful in custom
  contrast functions (`.emmc` functions).

- `emmeans::.get.excl(levs, exclude, include)` This is support for the
  `exclude` and `include` arguments of contrast functions. It checks
  legality and returns an integer vector of `exclude` indices in `levs`,
  given specified integer or character arguments `exclude` and
  `include`. In your `.emmc` function, `exclude` should default to
  `integer(0)` and `include` should have no default.

- `emmeans::.cmpMM(X, weights, assign)` creates a compact version of the
  model matrix `X` (or, preferably, its QR decomposition). This is
  useful if we want an
  [`emm_basis()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
  method to return a `model.matrix` element. The returned result is just
  the R portion of the QR decomposition of `diag(sqrt(weights)) %*% X`,
  with the `assign` attribute added. If `X` is a `qr` object, we assume
  the weights are already incorporated, as is true of the `qr` slot of a
  `lm` object.

[Back to Contents](#contents)

## Existing support for `rsm` objects

As a nontrivial example of how an existing package supports **emmeans**,
we show the support offered by the **rsm** package. Its `rsm` function
returns an `rsm` object which is an extension of the `lm` class. Part of
that extension has to do with `coded.data` structures whereby, as is
typical in response-surface analysis, models are fitted to variables
that have been linearly transformed (coded) so that the scope of each
predictor is represented by plus or minus 1 on the coded scale.

Without any extra support in **rsm**, `emmeans` will work just fine with
`rsm` objects; but if the data are coded, it becomes awkward to present
results in terms of the original predictors on their original, uncoded
scale. The `emmeans`-related methods in **rsm** provide a `mode`
argument that may be used to specify whether we want to work with coded
or uncoded data. The possible values for `mode` are `"asis"` (ignore any
codings, if present), `"coded"` (use the coded scale), and `"decoded"`
(use the decoded scale). The first two are actually the same in that no
decoding is done; but it seems clearer to provide separate options
because they represent two different situations.

### The `recover_data` method

Note that coding is a *predictor* transformation, not a response
transformation (we could have that, too, as it’s already supported by
the **emmeans** infrastructure). So, to handle the `"decode"` mode, we
will need to actually decode the predictors used to construct he
reference grid. That means we need to make `recover_data` a lot fancier!
Here it is:

``` r
recover_data.rsm = function(object, data, mode = c("asis", "coded", "decoded"), ...) {
    mode = match.arg(mode)
    cod = rsm::codings(object)
    fcall = object$call
    if(is.null(data))                                                 # 5
        data = emmeans::recover_data(fcall, 
                   delete.response(terms(object)), object$na.action, 
                   weights = weights(object), ...)
    if (!is.null(cod) && (mode == "decoded")) {
        pred = cpred = attr(data, "predictors")
        trms = attr(data, "terms")                                    #10
        data = rsm::decode.data(rsm::as.coded.data(data, formulas = cod))
        for (form in cod) {
            vn = all.vars(form)
            if (!is.na(idx <- grep(vn[1], pred))) { 
                pred[idx] = vn[2]                                     #15
                cpred = setdiff(cpred, vn[1])
            }
        }
        attr(data, "predictors") = pred
        new.trms = update(trms, reformulate(c("1", cpred)))           #20
        attr(new.trms, "orig") = trms
        attr(data, "terms") = new.trms
        attr(data, "misc") = cod
    }
    data
}
```

Lines 2–7 ensure that `mode` is legal, retrieves the codings from the
object, and obtain the results we would get from `recover_data` had it
been an `lm` object. If `mode` is not `"decoded"`, *or* if no codings
were used, that’s all we need. Otherwise, we need to return the decoded
data. However, it isn’t quite that simple, because the model equation is
still defined on the coded scale. Rather than to try to translate the
model coefficients and covariance matrix to the decoded scale, we
elected to remember what we will need to do later to put things back on
the coded scale. In lines 9–10, we retrieve the attributes of the
recovered data that provide the predictor names and `terms` object on
the coded scale. In line 11, we replace the recovered data with the
decoded data.

By the way, the codings comprise a list of formulas with the coded name
on the left and the original variable name on the right. It is possible
that only some of the predictors are coded (for example, blocking
factors will not be). In the `for` loop in lines 12–18, the coded
predictor names are replaced with their decoded names. For technical
reasons to be discussed later, we also remove these coded predictor
names from a copy, `cpred`, of the list of all predictors in the coded
model. In line 19, the `"predictors"` attribute of `data` is replaced
with the modified version.

Now, there is a nasty technicality. The `ref_grid` function in
**emmeans** has a few lines of code after `recover_data` is called that
determine if any terms in the model convert covariates to factors or
vice versa; and this code uses the model formula. That formula involves
variables on the coded scale, and those variables are no longer present
in the data, so an error will occur if it tries to access them. Luckily,
if we simply take those terms out of the formula, it won’t hurt because
those coded predictors would not have been converted in that way. So in
line 20, we update `trms` with a simpler model with the coded variables
excluded (the intercept is explicitly included to ensure there will be a
right-hand side even is `cpred` is empty). We save that as the `terms`
attribute, and the original terms as a new `"orig"` attribute to be
retrieved later. The `data` object, modified or not, is returned. If
data have been decoded, `ref_grid` will construct its grid using decoded
variables.

In line 23, we save the codings as the `"misc"` attribute, to be
accessed later by
[`emm_basis()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md).

### The `emm_basis` method

Now comes the `emm_basis` method that will be called after the grid is
defined. It is listed below:

``` r
emm_basis.rsm = function(object, trms, xlev, grid, 
                         mode = c("asis", "coded", "decoded"), misc, ...) {
    mode = match.arg(mode)
    cod = misc
    if(!is.null(cod) && mode == "decoded") {                          # 5
        grid = rsm::coded.data(grid, formulas = cod)
        trms = attr(trms, "orig")
    }
    
    m = model.frame(trms, grid, na.action = na.pass, xlev = xlev)     #10
    X = model.matrix(trms, m, contrasts.arg = object$contrasts)
    bhat = as.numeric(object$coefficients) 
    V = emmeans::.my.vcov(object, ...)
    
    if (sum(is.na(bhat)) > 0)                                         #15
        nbasis = estimability::nonest.basis(object$qr)
    else
        nbasis = estimability::all.estble
    dfargs = list(df = object$df.residual)
    dffun = function(k, dfargs) dfargs$df                             #20

    list(X = X, bhat = bhat, nbasis = nbasis, V = V, 
         dffun = dffun, dfargs = dfargs, misc = list())
}
```

This is much simpler. The coding formulas are obtained from `misc` (line
4) so that we don’t have to re-obtain them from the object. All we have
to do is determine if decoding was done (line 5); and, if so, convert
the grid back to the coded scale (line 6) and recover the original
`terms` attribute (line 7). The rest is borrowed directly from the
`emm_basis.lm` method in **emmeans**. Note that line 13 uses one of the
exported functions we described in the preceding section. Lines 15–18
use functions from the **estimability** package to handle the
possibility that the model is rank-deficient.

### A demonstration

Here’s a demonstration of this **rsm** support. The standard example for
`rsm` fits a second-order model `CR.rs2` to a dataset organized in two
blocks and with two coded predictors.

``` r
library("rsm")
example("rsm")   ### (output is not shown) ###
```

First, let’s look at some results on the coded scale—which are the same
as for an ordinary `lm` object.

``` r
emmeans(CR.rs2, ~ x1 * x2, mode = "coded", 
        at = list(x1 = c(-1, 0, 1), x2 = c(-2, 2)))
```

``` ro
##  x1 x2 emmean    SE df lower.CL upper.CL
##  -1 -2   75.0 0.298  7     74.3     75.7
##   0 -2   77.0 0.240  7     76.4     77.5
##   1 -2   76.4 0.298  7     75.6     77.1
##  -1  2   76.8 0.298  7     76.1     77.5
##   0  2   79.3 0.240  7     78.7     79.9
##   1  2   79.2 0.298  7     78.5     79.9
## 
## Results are averaged over the levels of: Block 
## Confidence level used: 0.95
```

Now, the coded variables `x1` and `x2` are derived from these coding
formulas for predictors `Time` and `Temp`:

``` r
codings(CR.rs1)
```

``` ro
## $x1
## x1 ~ (Time - 85)/5
## 
## $x2
## x2 ~ (Temp - 175)/5
```

Thus, for example, a coded value of `x1 = 1` corresponds to a time of
85 + 1 x 5 = 90. Here are some results working with decoded predictors.
Note that the `at` list must now be given in terms of `Time` and `Temp`:

``` r
emmeans(CR.rs2, ~ Time * Temp, mode = "decoded", 
        at = list(Time = c(80, 85, 90), Temp = c(165, 185)))
```

``` ro
##  Time Temp emmean    SE df lower.CL upper.CL
##    80  165   75.0 0.298  7     74.3     75.7
##    85  165   77.0 0.240  7     76.4     77.5
##    90  165   76.4 0.298  7     75.6     77.1
##    80  185   76.8 0.298  7     76.1     77.5
##    85  185   79.3 0.240  7     78.7     79.9
##    90  185   79.2 0.298  7     78.5     79.9
## 
## Results are averaged over the levels of: Block 
## Confidence level used: 0.95
```

Since the supplied settings are the same on the decoded scale as were
used on the coded scale, the EMMs are identical to those in the previous
output.

## Dispatching and restrictions

The **emmeans** package has internal support for a number of model
classes. When
[`recover_data()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
and
[`emm_basis()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
are dispatched, a search is made for external methods for a given class;
and if found, those methods are used instead of the internal ones.
However, certain restrictions apply when you aim to override an existing
internal method:

1.  The class name being extended must appear in the first or second
    position in the results of `class(object)`. That is, you may have a
    base class for which you provide
    [`recover_data()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
    and
    [`emm_basis()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
    methods, and those will also work for *direct* descendants thereof;
    but any class in third place or later in the inheritance is ignored.
2.  Certain classes vital to the correct operation of the package, e.g.,
    `"lm"`, `"glm"`, etc., may not be overridden.

If there are no existing internal methods for the class(es) you provide
methods for, there are no restrictions on them.

## Exporting and registering your methods

To make your methods available to users of your package, the methods
must be exported. R and CRAN are evolving in a way that having S3
methods in the registry is increasingly important; so it is a good idea
to provide for that. The problem is not all of your package users will
have **emmeans** installed.

Thus, registering the methods must be done conditionally. We provide a
courtesy function
[`.emm_register()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
to make this simple. Suppose that your package offers two model classes
`foo` and `bar`, and it includes the corresponding functions
`recover_data.foo`, `recover_data.bar`, `emm_basis.foo`, and
`emm_basis.bar`. Then to register these methods, add or modify the
`.onLoad` function in your package (traditionally saved in the source
file `zzz.R`):

``` r
.onLoad <- function(libname, pkgname) {
    if (requireNamespace("emmeans", quietly = TRUE))
        emmeans::.emm_register(c("foo", "bar"), pkgname)
}
```

You should also add `emmeans (>= 1.4)` and `estimability` (which is
required by **emmeans**) to the `Suggests` field of your `DESCRIPTION`
file.

When registering a `qdrg` method, do the same as shown above, but add
the argument `qdrg = TRUE` to the
[`.emm_register()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
call, and in `Suggests`, use `emmeans (>= 1.12)` as well as
`estimability`.

[Back to Contents](#contents)

## Conclusions

It is relatively simple to write appropriate methods that work with
**emmeans** for model objects it does not support. I hope this vignette
is helpful for understanding how. Furthermore, if you are the developer
of a package that fits linear models, I encourage you to include
`recover_data` and `emm_basis` methods for those classes of objects, so
that users have access to **emmeans** support.

[Back to Contents](#contents)

[Index of all vignette
topics](https://rvlenth.github.io/emmeans/articles/vignette-topics.md)
