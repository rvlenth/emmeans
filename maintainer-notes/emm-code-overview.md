---
title: "Code overview"
author: "Russ Lenth"
date: "2025-10-29"
output: emmeans::.emm_vignette
---

Here we have a bit about the main functions and how they work. This is very, very far from
a line-by=line description of what's going on in the code, but it is hoped that
this helps explain the fundamentals.

## Contents
  * [`ref_grid`](#ref-grid)
  * [`emmeans`](#emmeans)
  * [Nested effects](#nesting)
  * [`contrast`](#contrast)
  * [`.find.by.rows`](#byrows)
  * [`summary`](#summary)
  * [`confint` and `test`](#confint)
  * [`joint_tests`](#joint)
  * [`.est.se.df`](#est)
  * [`emtrends`](#emtrends)
  * [`regrid`](#regrid)
  * [Bayesian models, `summary.hpd`](#bayes)
  * [Bias adjustment](#bias)
  * [`mvcontrast`](#mvcontrast)
  * [Satterthwaite method](#satt)
  * [Estimability](#estble)


## `ref_grid` {#ref-grid}
This is the foundation of the package. It looks at the model and extracts the necessary 
info to create an `emmGrid` object. The slots are

  * `bhat` -- $b$: the fixed-effects coefficients:
  * `V` -- $V$: the variance-covariance matrix of $b$
  * `linfct` -- $L$: the linear functions of $b$, i.e., the linear predictor for the 
     points in the reference grid is $Lb$
  * `grid`: a data frame with the factor/covariate combinations corresponding to each row of $L$.
    It also usually has a column named `.wgt.` (weights for each point, used when weighted means
    are requested), and sometimes an `.offset.` column (if so, it is added to $Lb$)
  * `levels`: a named list of levels. Basically, `@grid = do.call(expand.grid, @levels)`
  * `dffun` and `dfargs`: a function and named list of arguments for computing the 
    degrees of freedom associated with a linear prediction $x'b$
  * `nbasis` -- $N$: either a $1\times1$ matrix of `NA` or a matrix whose columns span the null space 
    of the coefficient matrix underlying $b$. It is used for assessing estimability of
    $x'b$; if $N'x \ne 0$, then $x'b$ is not uniquely estimable.
  * `post.beta`: for Bayesian models, this is a posterior sample of $b$ instances
  * `model.info`, `roles`, `matlevs`: lists with other model information
  * `misc`: a list of extra information used by package functions. See documentation for `update.emmGrid`

`ref_grid()` calls a `recover_data()` method for the model to recover the data
used to fit the model, looks also at the model's `terms`, and figures out what
variables are involved, which are factors, and what are the levels for each
predictor. It then creates `@grid` and passes this to the `emm_basis()` method
for the model, which returns `bhat`, `V`, and `X -> linfct`, `dffun`, `dfargs`,
and maybe some elements of `misc` such as the link function. Then it figures out 
if there was a response transformation. If the response is multivariate, then `X`
is really just the matrix for each response, while `bhat` is a stretched-out matrix;
so we figure out new variable name(s) and levels for the multivariate response,
and expand `X` via a kronecker product, and expand `grid` accordingly as well.

## `emmeans` {#emmeans}
`emmeans` creates a new `emmGrid` object corresponding to marginal means of the
reference grid. (If provided the model instead, it calls `ref_grid()`.) It does
this by averaging the rows of `linfct` in the same way. To keep track of
everything properly, we first ensure that `grid` (and `linfct`) is sorted in
standard order (first factor varying fastest, last the slowest). (It depends on
the grid being regular; and if it isn't, an error is thrown.) Then we create the
index vector $1, 2, ..., n$ of row numbers, and puts this into an `array` with
dimensions equal to the lengths of `levels`; then by identifying which
dimensions need to be averaged over, that also identifies the row indexes of
`linfct` that we need to average over, possibly with weights. The resturned
`object` is basically the same as the input one, except for `grid` and `linfct`.

### Nested fixed effects (and `.nested_emm`) {#nesting}
A nuance here is that if the model has nested fixed effects, we store the
details in `model.info`. Depending on how the marginal means desired in
`emmeans()` relate to the nesting pattern, `emmeans()` may need to be called
with each level of the nesting factor(s). Meanwhile, because we need a regular
grid in `emmeans()`, the reference grid is created with all predictors crossed,
and a logical vector `misc$display` is created that is `TRUE` for the rows that
exist in a nest, and `FALSE` if not.

## `contrast` {#contrast}
This function creates a new `emmGrid` object from the supplied one, replacing
`linfct` ($L$) by $ML$, where $M$ is a matrix of contrast coefficients. There
are a bunch of standard contrast families that are implemented as `.emmc`
functions (e.g., `pairwise.emmc()`) whose arguments are the factor levels and
perhaps some other arguments (that all have to have defaults). The function returns 
a data frame such that each column has the coefficients of a linear function to
apply to the rows of $L$ (or some subset thereof, as determined by `by`). When 
we have nested effects, there is a `.nested_contrast` function that does this for
each nest. The `.emmc` function also returns the names to assign the contrast (and
this modifies `levels` and grid` appropriately), and perhaps a default `adjust`
method that is saved in `misc`.

### `.find.by.rows` {#byrows}
When there is a non-trivial `by` specification in `contrast` (and some other functions),
we need to identify which rows of the input `emmGrid` object correspond to each
combination of levels of the `by` factors, and that is the purpose of the `.find.by.rows()` 
function. It returns a list of integer vectors for the respective `by` groups.

## `summary` {#summary}
This function does the actual estimation and returns an object of class
`emm_summary`, which extends `data.frame`. Its `print` method formats the
results and displays them in much the way that a data frame is, but often with
extra annotations (from `misc$mesg`). Its argument `infer` is a logical vector
of length 2 that decides whether to also show confidence intervals or tests/P
values, respectively (if of length only 1, that value is used for both). The
default for `infer` is in `misc$infer` which is initialized to `(FALSE, FALSE)`
by `ref_grid()`, `c(TRUE, FALSE)` by `emmeans()`, and `c(FALSE, TRUE)` by `contrast()`.

`summary` calls `.est.se.df()` (see below) and if `type = "response"`, the
estimates are back-transformed via `link$linkinv(est)` and the SEs are
multiplied by `link$mu.eta(est)` (this is the delta method). If `adjust` is
other than `"none"`, we determine whether it is a "legal" adjustment (for
example, the Tukey adjustment is not legal unless we have exactly *one* family
of pairwise comparisons, and there are bookkeeping provisions to ascertain
that); then for confidence intervals, critical values are adjusted, and for
tests, P values are adjusted.

The `summary` function also allows for a nonzero null hypothesis and/or one-sided tests,
and implements the various significance and equivalence tests in the documentation.

### `confint` {#confint}
This just calls `summary` with `infer = c(TRUE, FALSE)`.

### `test` {#test}
This calls `summary` with `infer = c(FALSE, TRUE)`. But it also has a `joint`
argument that if `TRUE`, computes a joint test of all rows of `linfct`. Note
that $\mbox{cov}(Lb) = LVL'$ but we need to take extra care to account for any
non-estimable or linearly-dependent rows. Assuming we've done that, the Wald
statistic is $(Lb - \mu_0)'(LVL')^{-1}(Lb - \mu_0)$ and the numerator degrees of
freedom is the rank of $L$. We kind of have to punt for the denominator d.f.
since the `dffun` slot doesn't account for multivariate cases. Taking care of 
non-estimability entails a projection from the **estimability** package, and
linear dependencies are handled via the `qr()` function.

### `joint_tests` {#joint}
This function goes through `linfct` with respect to the model terms and breaks
it down into independent pieces relating to contrasts of each term, and calls 
`test(..., joint = TRUE)` for each piece. Then it assembles the results into a 
`summary_emm` object. By default, any terms that have zero d.f. are omitted.

Note that, `joint_tests` returns no tests of covariates if just called
with a default reference grid or `emmeans` result. If we call this with a model
object, then it calls `ref_grid` with a different default for covariates,
reducing them to an interval around their mean instead of just their mean. This
works correctly so long as the covariate effects are linear. We need more than
the two values if nonlinear.

### `.est.se.df` {#est}
Is a primary workhorse for `summary`. Per its name, it puts together the linear 
predictions $Lb$, standard errors $\sqrt{\mbox{diag}(LVL')}$, and degrees of freedom
via `dffun` and `dfargs`, and determines the link function, if present, and bookkeeping
factors for multiplicity adjustments.

## `emtrends` {#emtrends}
This function requires a model object, and returns a special kind of reference
grid based on difference ratios of the covariate `var`. It uses `ref_grid` to do
most of the work, and there are special hooks there to hand it an expanded set
of `var` values, which it then uses for the difference quotients. It allows for
higher-order polynomial effects, using Newton's divided-difference formulas.

## `regrid` {#regrid}
This function completely overhauls the reference grid, basically divorcing it
from $b$ and $L$. If called with `transform = "response"`, it replaces `bhat`
with `h($Lb$)` (where `h` is `link$linkinv`), `linfct` with the identity matrix,
and `V` with $DLVL'D$ where $D$ is the diagonal matrix of `link$mu.eta` values.
The `transform` argument also allows it to be transformed to other transformed
scales using the same ideas in reverse, after first transforming to the response
scale. With `transform = "none"`, we do like `transform = "response"` but using
`h()` as the identity function. And with `transform = "pass"` we change nothing
unless `N.sim` is non-missing, in which case the `post.beta` slot is added with
simulated $b$ values according to the `sim` argument.

## Bayesian models (or simulated reference grids); `summary.hpd`{#bayes}
When `post.beta` is non-`NA`,`emmeans`, `contrast`, `emtrends`, and `regrid`
apply whatever they did to `bhat` to each row of `post.beta`, and return the
object with `post.beta` replaced by that result. When `summary` is called and 
`post.beta` is not `NA`, it diverts to `summary.hpd` unless called with 
`frequentist = TRUE`; in that case, the reference grid is usually initialized with
`bhat` as the average of the rows and `V = cov(post.beta)`.

In some ways, Bayesian models are easier to support, because all we *really* need
is the posterior sample of estimates, and with back-transformations and such
we just compute thenappropriate posterior sample of back-transformed estimates, 
and summarize with `summary.hpd`. We don't need the delta method.

## Bias adjustment {#bias}
This is explained in the documentation for `summary`. It is implemented by doing
a fixup to `link` where we replace `link$linkinv` by $h(\eta)+\frac12h''(\eta)$,
the latter term being estimated using a difference quotient on `link$mu.eta`.

## `mvcontrast` {#mvcontrast}
This is like contrast except we identify one factor to treat as multivariate
and then we use the specified contrast method on that multivariate vector
and test the result using Hotelling's $T^2$.

## Satterthwaite method {#satt}
First, I have to say that Iowa has special ownership of this since Satterthwaite 
earned his PhD in our department (well, actually Math at the time) in 1943 by
developing this method. It's about the only thing he ever did though, as he
suffered from schizophrenia and just wasn't able to function well most of the rest
of his life. Anyway, the idea is that we have a variance estimator $W$ and we like
the idea that it would be proportional to a $\chi^2$ r.v. So we find the d.f. based on
the first two moments. Now if $W/c \sim \chi^2_\nu$, then $E(W)=c\nu$ and $\mbox{var}(W) = 2c^2\nu$,
implying that $2\cdot[E(W)]^2/\mbox{var(W)} = (2c^2\nu^2)/(2c^2\nu) = \nu$.
Accordingly, we estimate the degrees of freedom using  $\hat\nu = 2W^2 / \hat{\mbox{var}}(W)$.

In the support functions for `nlme::gls` objects, it is feasible to do this by
obtaining the jacobian of the variance matrix of the random effects and using
that to estimate the variance of the variance estimate in question. This is the
function `gls_grad`. For `lme4::mermod` and `nlme::lme` objects, we don't have
that variance matrix, so instead we have an approximate method
(`"appx-satterthwaite"`) that perturbs the response values slightly and refits
the model. It only takes a few of these perturbations to do a respectable job of
estimating the required variances; this is done by the `gradV.kludge` function.

# Estimability {#estble}
An important part of the package is that we assess estimability of our predictions. 
This is important because some models allow rank-deficient model matrices; if we have a
rank-deficient model, there are infinitely many possible solutions for the fixed
effects, but we have only one of them. A prediction is defined as estimable if it is
the same, no matter which solution we used. Equivalently, given the model matrix $X$,
$x'\beta$ is estimable if, and only if, $x$ is in the row space of $X$.

The **estimability** packagen provides the functions needed to assess estimability. 
We do this by creating a basis $N$ for the null space of $X$, i.e., $XN=0$, and so
$x'\beta$ is estimable iff $x'N=0$. We store the matrix $N$ in the `@nbasis` slot
of an `emmGrid` object.

For more details, see the package documentation for **estimability** and its 
[vignette](https://cran.r-project.org/web/packages/estimability/vignettes/add-est-check.html).