# The `emmGrid` class

The `emmGrid` class encapsulates linear functions of regression
parameters, defined over a grid of predictors. This includes reference
grids and grids of marginal means thereof (aka estimated marginal
means). Objects of class \`emmGrid\` may be used independently of the
underlying model object. Instances are created primarily by
[`ref_grid`](https://rvlenth.github.io/emmeans/reference/ref_grid.md)
and [`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md),
and several related functions.

## Slots

- `model.info`:

  list. Contains the elements `call` (the call that produced the model),
  `terms` (its `terms` object), and `xlev` (factor-level information)

- `roles`:

  list. Contains at least the elements `predictors`, `responses`, and
  `multresp`. Each is a character vector of names of these variables.

- `grid`:

  data.frame. Contains the combinations of the variables that define the
  reference grid. In addition, there is an auxiliary column named
  `".wgt."` holding the observed frequencies or weights for each factor
  combination (excluding covariates). If the model has one or more
  [`offset()`](https://rdrr.io/r/stats/offset.html) calls, there is an
  another auxiliary column named `".offset."`. Auxiliary columns are not
  considered part of the reference grid. (However, any variables
  included in `offset` calls *are* in the reference grid.)

- `levels`:

  list. Each entry is a character vector with the distinct levels of
  each variable in the reference grid. Note that `grid` is obtained by
  applying the function
  [`expand.grid`](https://rdrr.io/r/base/expand.grid.html) to this list

- `matlevs`:

  list. Like `levels` but has the levels of any matrices in the original
  dataset. Matrix columns are always concatenated and treated as a
  single variable for purposes of the reference grid

- `linfct`:

  matrix. Each row consists of the linear function of the regression
  coefficients for predicting its corresponding element of the reference
  grid. The rows of this matrix go in one-to-one correspondence with the
  rows of `grid`, and the columns with elements of `bhat`.

- `bhat`:

  numeric. The regression coefficients. If there is a multivariate
  response, the matrix of coefficients is flattened to a single vector,
  and `linfct` and `V` redefined appropriately. Important: `bhat` must
  *include* any `NA` values produced as a result of collinearity in the
  predictors. These are taken care of later in the estimability check.

- `nbasis`:

  matrix. The basis for the non-estimable functions of the regression
  coefficients. Every EMM will correspond to a linear combination of
  rows of `linfct`, and that result must be orthogonal to all the
  columns of `nbasis` in order to be estimable. If everything is
  estimable, `nbasis` should be a 1 x 1 matrix of `NA`.

- `V`:

  matrix. The symmetric variance-covariance matrix of `bhat`

- `dffun`:

  function having two arguments. `dffun(k, dfargs)` should return the
  degrees of freedom for the linear function `sum(k*bhat)`, or `NA` if
  unavailable

- `dfargs`:

  list. Used to hold any additional information needed by `dffun`.

- `misc`:

  list. Additional information used by methods. These include at least
  the following: `estName` (the label for the estimates of linear
  functions), and the default values of `infer`, `level`, and `adjust`
  to be used in the
  [`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
  method. Elements in this slot may be modified if desired using the
  [`update.emmGrid`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md)
  method.

- `post.beta`:

  matrix. A sample from the posterior distribution of the regression
  coefficients, if MCMC methods were used; or a 1 x 1 matrix of `NA`
  otherwise. When it is non-trivial, the
  [`as.mcmc.emmGrid`](https://rvlenth.github.io/emmeans/reference/mcmc-support.md)
  method returns `post.beta %*% t(linfct)`, which is a sample from the
  posterior distribution of the EMMs.

## Methods

All methods for these objects are S3 methods except for `show`. They
include
[`[.emmGrid`](https://rvlenth.github.io/emmeans/reference/rbind.emmGrid.md),
[`as.glht.emmGrid`](https://rvlenth.github.io/emmeans/reference/glht-support.md),
[`as.mcmc.emmGrid`](https://rvlenth.github.io/emmeans/reference/mcmc-support.md),
[`as.mcmc.list.emmGrid`](https://rvlenth.github.io/emmeans/reference/mcmc-support.md)
(see coda),
[`cld.emmGrid`](https://rvlenth.github.io/emmeans/reference/CLD.emmGrid.md)
(see multcomp),
[`coef.emmGrid`](https://rvlenth.github.io/emmeans/reference/contrast.md),
[`confint.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md),
[`contrast.emmGrid`](https://rvlenth.github.io/emmeans/reference/contrast.md),
[`pairs.emmGrid`](https://rvlenth.github.io/emmeans/reference/contrast.md),
[`plot.emmGrid`](https://rvlenth.github.io/emmeans/reference/plot.md),
[`predict.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md),
[`print.emmGrid`](https://rvlenth.github.io/emmeans/reference/emmGrid-methods.md),
[`rbind.emmGrid`](https://rvlenth.github.io/emmeans/reference/rbind.emmGrid.md),
`show.emmGrid`,
[`str.emmGrid`](https://rvlenth.github.io/emmeans/reference/emmGrid-methods.md),
[`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md),
[`test.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md),
[`update.emmGrid`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md),
[`vcov.emmGrid`](https://rvlenth.github.io/emmeans/reference/emmGrid-methods.md),
and
[`xtable.emmGrid`](https://rvlenth.github.io/emmeans/reference/xtable.emmGrid.md)
