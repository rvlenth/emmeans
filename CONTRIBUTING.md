# NA

Welcome! We presume you’re reading this because you want to contribute
to **emmeans** (pronounced “E-M-means”). Open source software relies on
the user community lending their skills and time to improving the
ecosystem, so we welcome your interest.

If your interest is solely to determine how to have **emmeans** support
your linear modelling package, please see our guide [“Extending
emmeans”](https://rvlenth.github.io/emmeans/articles/xtending.html).

## Contents

- [Code of Conduct](#coc)
- [Noncoding way to contribute](#noncoding)
- [How to contribute](#how)
- [Setting up your development environment](#dev-env)
- [Coding style](#code-style)
- [Documentation and vignette style](#doc-style)
- [Unit tests](#testing)
- [Generative AI policy](#ai-policy)
- [Requests for new features](#feature-requests)
- [Code overview](#code-overview)
  - [The `ref_grid` function](#ref-grid)
  - [The `emmeans` function](#emmeans)
    - [Nested fixed effects](#nesting)
  - [The `contrast` function](#contrast)
    - [The `.find.by.rows` function](#byrows)
  - [The `summary` function](#summary)
    - [The `confint` and `test` functions](#confint)
    - [The `joint_tests` function](#joint)
    - [The `.est.se.df` function](#est)
  - [The `emtrends` function](#emtrends)
  - [The `regrid` function](#regrid)
  - [Bayesian models and the `summary.hpd` function](#bayes)
  - [Bias adjustment](#bias)
  - [The `mvcontrast` function](#mvcontrast)
  - [Satterthwaite method](#satt)
  - [Estimability](#estble)

## Code of Conduct

All contributers must follow the Contributor Covenant [contributer code
of
conduct](https://www.contributor-covenant.org/version/2/1/code_of_conduct/)
in order to contribute to **emmeans**.

## Noncoding ways to contribute

1.  Make improvements to the documentation. Perhaps the help files need
    expansion or you have a suggestion for a vignette. At the time of
    publishing this, we do have an [open
    request](https://github.com/rvlenth/emmeans/issues/521) to create a
    cheatsheet for emmeans. It would be a nice helper, but we do not
    have time to work on this.

2.  Answer an [issue](https://github.com/rvlenth/emmeans/issues).
    **emmeans** does not tend to have a large number of open issues, but
    this is generally good advice for supporting packages.

## How to contribute

1.  Read through the [code overview](#code-overview).

2.  File an issue before proceeding with any bug fix or feature
    addition. This is less about asking permission and more for
    coordination.

3.  Once changes are agreed upon, fork the repository and make a new
    branch for the planned change.

4.  After implementing a change, please run existing unit tests
    ([`testthat::test_check()`](https://testthat.r-lib.org/reference/test_package.html))
    to ensure existing functionality still works. Please update the
    NEWS.md file with your change(s). Please do not update the version
    or date in the DESCRIPTION file, but *do* add your name to the end
    of the authors field with the label “ctb” indicating you are a
    contributor to emmeans.

5.  Make a pull request to the main branch with a reasonably descriptive
    message describing the changes made and linking to the issue filed.

## Setting up your development environment

Download the latest version of emmeans from Github:

    remotes::install_github("rvlenth/emmeans", dependencies = TRUE)

Alternatively, **devtools** has the same function if you prefer that
package.

Make sure you are using a version of R and **roxygen2** that is
supported by **emmeans**; see the `DESCRIPTION` file for this
information.

## Coding style

This package has a few conventions that depart from the default
expectations established by the
[lintr](https://lintr.r-lib.org/index.html) package:

- We prefer that new functions use camel case,
  e.g. [`joint_tests()`](https://rvlenth.github.io/emmeans/reference/joint_tests.md),
  but some older function names use dotted case, e.g
  [`make.meanint()`](https://rvlenth.github.io/emmeans/reference/joint_tests.md).

- Default indentation is 4 spaces.

- Acceptable assignment operators: “=”, “\<-”, and “\<\<-”

- Curly braces: the opening brace is the last character of a line but it
  does not need its own line; the closing brace needs its own line.

- `else {` always starts on a new line ()

Run `lintr::lint_dir(path = c("R", "tests")` to check your style. A
.lintr file is included in the **emmeans** root directory for
**emmeans** custom style choices. For existing **emmeans** R/Rmd files,
you will see lint errors.

## Documentation and vignette style

Like many R packages, **emmeans** uses
[**roxygen2**](https://roxygen2.r-lib.org/) for function documentation.
Please read up on that package prior to making changes to function
argument descriptions and other aspects of function help files.

**emmeans** has a large number of vignettes, and there is likely a need
for more! If there is a topic you would like to see addressed, file an
issue. We may not be able to work on it, but some other interested
**emmeans** user (perhaps you!) might. All vignettes are written with
Rmarkdown (with the .Rmd file extension) and use this YAML header and
opening chunk:

    ---
    title: "{A Brief Title Describing the Vignette}"
    author: "emmeans package, Version `r packageVersion('emmeans')`"
    output: emmeans::.emm_vignette
    vignette: >
      %\VignetteIndexEntry{{A brief description}}
      %\VignetteEngine{knitr::rmarkdown}
      %\VignetteEncoding{UTF-8}
    ---

```` R
```{r, echo = FALSE, results = "hide", message = FALSE}
require("emmeans")
knitr::opts_chunk$set(fig.width = 4.5, class.output = "ro")
```
````

There may be additional information in the opening R chunk, depending on
the purpose of the vignette (e.g. other packages are loaded). A common
addition is `options(show.signif.stars = FALSE, width = 100)` indicating
that the emmeans documentation does not typically show significance
stars.

Please create a section at the top indexing all the level 2 (`##`)
headers in the vignette with their html anchors so each section has a
direct web link.

There is a custom CSS file for **emmeans** that
[`emmeans::.emm_vignette`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
will use when building the documentation.

Please see the [vignettes
directory](https://github.com/rvlenth/emmeans/tree/main/vignettes) for
specific examples.

## Testing

There is a need for additional unit tests of the package code, as some
parts of **emmeans** have low coverage. If you do want to contribute a
new feature, then please expect to write new unit tests to verify their
behavior and outputs under expected and extreme cases.

## Generative AI policy

Generative AI (“genAI”), referring to large language models, chatbots,
agents and similar, is a useful tool for generating code, but there are
serious concerns regarding quality of code that is generated by genAI,
energy usage required for genAI tools, and copyright violations to
construct large language models (LLMs).

- We have found it sufficiently useful that AI-assisted bug checking is
  welcome.

- At this time, fully AI-generated code (including unit tests) are not
  acceptable due to unpredictable quality.

- Perhaps LLMs (at least the cheap ones most of us can access) will
  improve in quality and we can visit this policy, but at this time, we
  want a human hand engaged in improving **emmeans**.

## Feature extension requests

For us to incorporate a new function into the package, we need to be
convinced that you have addressed the following points:

1.  The contribution is useful to more than a very narrow group of
    users. Can you describe how broadly your contribution might be used?

2.  We thoroughly understand its underpinnings. This is especially
    important, because any questions about the use of the function, or
    any bug reports, come to us, the maintainers of the package. Can you
    give literature reference(s) for the methods you implemented?

3.  How easily could it be inadvertently misused? For example, is it
    appropriate for use with only certain models of a given class, but
    not others?

4.  What about non-estimability issues? Can this arise? How do we handle
    this?

5.  Generally, what else can go wrong, even with intentional use? For
    example, is it robust to irregularities that might exist in the data
    or model?

## Code Overview

Here we have a bit about the main functions and how they work. This is
very, very far from a line-by-line description of what’s going on in the
code, but it is hoped that this helps explain the fundamentals.

### The `ref_grid` function

This is the foundation of the package. It looks at the model and
extracts the necessary info to create an `emmGrid` object. The slots are

- `bhat` – \\b\\: the fixed-effects coefficients:
- `V` – \\V\\: the variance-covariance matrix of \\b\\
- `linfct` – \\L\\: the linear functions of \\b\\, i.e., the linear
  predictor for the points in the reference grid is \\Lb\\
- `grid`: a data frame with the factor/covariate combinations
  corresponding to each row of \\L\\. It also usually has a column named
  `.wgt.` (weights for each point, used when weighted means are
  requested), and sometimes an `.offset.` column (if so, it is added to
  \\Lb\\)
- `levels`: a named list of levels. Basically,
  `@grid = do.call(expand.grid, @levels)`
- `dffun` and `dfargs`: a function and named list of arguments for
  computing the degrees of freedom associated with a linear prediction
  \\x'b\\
- `nbasis` – \\N\\: either a \\1\times1\\ matrix of `NA` or a matrix
  whose columns span the null space of the coefficient matrix underlying
  \\b\\. It is used for assessing estimability of \\x'b\\; if \\N'x \ne
  0\\, then \\x'b\\ is not uniquely estimable.
- `post.beta`: for Bayesian models, this is a posterior sample of \\b\\
  instances
- `model.info`, `roles`, `matlevs`: lists with other model information
- `misc`: a list of extra information used by package functions. See
  documentation for `update.emmGrid`

[`ref_grid()`](https://rvlenth.github.io/emmeans/reference/ref_grid.md)
calls a
[`recover_data()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
method for the model to recover the data used to fit the model, looks
also at the model’s `terms`, and figures out what variables are
involved, which are factors, and what are the levels for each predictor.
It then creates `@grid` and passes this to the
[`emm_basis()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
method for the model, which returns `bhat`, `V`, and `X -> linfct`,
`dffun`, `dfargs`, and maybe some elements of `misc` such as the link
function. Then it figures out if there was a response transformation. If
the response is multivariate, then `X` is really just the matrix for
each response, while `bhat` is a stretched-out matrix; so we figure out
new variable name(s) and levels for the multivariate response, and
expand `X` via a kronecker product, and expand `grid` accordingly as
well.

### The `emmeans` function

`emmeans` creates a new `emmGrid` object corresponding to marginal means
of the reference grid. (If provided the model instead, it calls
[`ref_grid()`](https://rvlenth.github.io/emmeans/reference/ref_grid.md).)
It does this by averaging the rows of `linfct` in the same way. To keep
track of everything properly, we first ensure that `grid` (and `linfct`)
is sorted in standard order (first factor varying fastest, last the
slowest). (It depends on the grid being regular; and if it isn’t, an
error is thrown.) Then we create the index vector \\1, 2, ..., n\\ of
row numbers, and puts this into an `array` with dimensions equal to the
lengths of `levels`; then by identifying which dimensions need to be
averaged over, that also identifies the row indexes of `linfct` that we
need to average over, possibly with weights. The returned `object` is
basically the same as the input one, except for `grid` and `linfct`.

#### Nested fixed effects (and the `.nested_emm` function)

A nuance here is that if the model has nested fixed effects, we store
the details in `model.info`. Depending on how the marginal means desired
in [`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.md)
relate to the nesting pattern,
[`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.md)
may need to be called with each level of the nesting factor(s).
Meanwhile, because we need a regular grid in
[`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.md),
the reference grid is created with all predictors crossed, and a logical
vector `misc$display` is created that is `TRUE` for the rows that exist
in a nest, and `FALSE` if not.

### The `contrast` function

This function creates a new `emmGrid` object from the supplied one,
replacing `linfct` (\\L\\) by \\ML\\, where \\M\\ is a matrix of
contrast coefficients. There are a bunch of standard contrast families
that are implemented as `.emmc` functions (e.g.,
[`pairwise.emmc()`](https://rvlenth.github.io/emmeans/reference/emmc-functions.md))
whose arguments are the factor levels and perhaps some other arguments
(that all have to have defaults). The function returns a data frame such
that each column has the coefficients of a linear function to apply to
the rows of \\L\\ (or some subset thereof, as determined by `by`). When
we have nested effects, there is a `.nested_contrast` function that does
this for each nest. The `.emmc` function also returns the names to
assign the contrast (and this modifies `levels` and
grid`appropriately), and perhaps a default`adjust`method that is saved in`misc\`.

#### The `.find.by.rows` function

When there is a non-trivial `by` specification in `contrast` (and some
other functions), we need to identify which rows of the input `emmGrid`
object correspond to each combination of levels of the `by` factors, and
that is the purpose of the `.find.by.rows()` function. It returns a list
of integer vectors for the respective `by` groups.

### The `summary` function

This function does the actual estimation and returns an object of class
`emm_summary`, which extends `data.frame`. Its `print` method formats
the results and displays them in much the way that a data frame is, but
often with extra annotations (from `misc$mesg`). Its argument `infer` is
a logical vector of length 2 that decides whether to also show
confidence intervals or tests/P values, respectively (if of length only
1, that value is used for both). The default for `infer` is in
`misc$infer` which is initialized to `(FALSE, FALSE)` by
[`ref_grid()`](https://rvlenth.github.io/emmeans/reference/ref_grid.md),
`c(TRUE, FALSE)` by
[`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.md),
and `c(FALSE, TRUE)` by
[`contrast()`](https://rvlenth.github.io/emmeans/reference/contrast.md).

`summary` calls `.est.se.df()` (see below) and if `type = "response"`,
the estimates are back-transformed via `link$linkinv(est)` and the SEs
are multiplied by `link$mu.eta(est)` (this is the delta method). If
`adjust` is other than `"none"`, we determine whether it is a “legal”
adjustment (for example, the Tukey adjustment is not legal unless we
have exactly *one* family of pairwise comparisons, and there are
bookkeeping provisions to ascertain that); then for confidence
intervals, critical values are adjusted, and for tests, P values are
adjusted.

The `summary` function also allows for a nonzero null hypothesis and/or
one-sided tests, and implements the various significance and equivalence
tests in the documentation.

#### The `confint` and `test` functions

`confint` just calls `summary` with `infer = c(TRUE, FALSE)`.

`test` calls `summary` with `infer = c(FALSE, TRUE)`. But it also has a
`joint` argument that if `TRUE`, computes a joint test of all rows of
`linfct`. Note that \\\mbox{cov}(Lb) = LVL'\\ but we need to take extra
care to account for any non-estimable or linearly-dependent rows.
Assuming we’ve done that, the Wald statistic is \\(Lb -
\mu_0)'(LVL')^{-1}(Lb - \mu_0)\\ and the numerator degrees of freedom is
the rank of \\L\\. We kind of have to punt for the denominator d.f.
since the `dffun` slot doesn’t account for multivariate cases. Taking
care of non-estimability entails a projection from the **estimability**
package, and linear dependencies are handled via the
[`qr()`](https://rdrr.io/r/base/qr.html) function.

#### The `joint_tests` function

This function goes through `linfct` with respect to the model terms and
breaks it down into independent pieces relating to contrasts of each
term, and calls `test(..., joint = TRUE)` for each piece. Then it
assembles the results into a `summary_emm` object. By default, any terms
that have zero d.f. are omitted.

Note that, `joint_tests` returns no tests of covariates if just called
with a default reference grid or `emmeans` result. If we call this with
a model object, then it calls `ref_grid` with a different default for
covariates, reducing them to an interval around their mean instead of
just their mean. This works correctly so long as the covariate effects
are linear. We need more than the two values if nonlinear.

#### The `.est.se.df` function

Is a primary workhorse for `summary`. Per its name, it puts together the
linear predictions \\Lb\\, standard errors \\\sqrt{\mbox{diag}(LVL')}\\,
and degrees of freedom via `dffun` and `dfargs`, and determines the link
function, if present, and bookkeeping factors for multiplicity
adjustments.

### The `emtrends` function

This function requires a model object, and returns a special kind of
reference grid based on difference ratios of the covariate `var`. It
uses `ref_grid` to do most of the work, and there are special hooks
there to hand it an expanded set of `var` values, which it then uses for
the difference quotients. It allows for higher-order polynomial effects,
using Newton’s divided-difference formulas.

### The `regrid` function

This function completely overhauls the reference grid, basically
divorcing it from \\b\\ and \\L\\. If called with
`transform = "response"`, it replaces `bhat` with \\h(Lb)\\ (where \\h\\
is `link$linkinv`), `linfct` with the identity matrix, and `V` with
\\DLVL'D\\ where \\D\\ is the diagonal matrix of `link$mu.eta` values.
The `transform` argument also allows it to be transformed to other
transformed scales using the same ideas in reverse, after first
transforming to the response scale. With `transform = "none"`, we do
like `transform = "response"` but using \\h\\ as the identity function.
And with `transform = "pass"` we change nothing unless `N.sim` is
non-missing, in which case the `post.beta` slot is added with simulated
\\b\\ values according to the `sim` argument.

### Bayesian models (or simulated reference grids); `summary.hpd` function

When `post.beta` is non-`NA`,`emmeans`, `contrast`, `emtrends`, and
`regrid` apply whatever they did to `bhat` to each row of `post.beta`,
and return the object with `post.beta` replaced by that result. When
`summary` is called and `post.beta` is not `NA`, it diverts to
`summary.hpd` unless called with `frequentist = TRUE`; in that case, the
reference grid is usually initialized with `bhat` as the average of the
rows and `V = cov(post.beta)`.

In some ways, Bayesian models are easier to support, because all we
*really* need is the posterior sample of estimates, and with
back-transformations and such we just compute the appropriate posterior
sample of back-transformed estimates, and summarize with `summary.hpd`.
We don’t need the delta method.

### Bias adjustment

This is explained in the documentation for `summary`. It is implemented
by doing a fixup to `link` where we replace `link$linkinv` by
\\h(\eta)+\frac12h''(\eta)\\, the latter term being estimated using a
difference quotient on `link$mu.eta`.

### The `mvcontrast` function

This is like contrast except we identify one factor to treat as
multivariate and then we use the specified contrast method on that
multivariate vector and test the result using Hotelling’s \\T^2\\.

### Satterthwaite method

The idea is that we have a variance estimator \\W\\ and we posit that it
would be proportional to a \\\chi^2\\ r.v. So we find the d.f. based on
the first two moments. Now if \\W/c \sim \chi^2\_\nu\\, then
\\E(W)=c\nu\\ and \\\mbox{var}(W) = 2c^2\nu\\, implying that
\\2\cdot\[E(W)\]^2/\mbox{var(W)} = (2c^2\nu^2)/(2c^2\nu) = \nu\\.
Accordingly, we estimate the degrees of freedom using \\\hat\nu = 2W^2 /
\hat{\mbox{var}}(W)\\.

In the support functions for
[`nlme::gls`](https://rdrr.io/pkg/nlme/man/gls.html) objects, it is
feasible to do this by obtaining the jacobian of the variance matrix of
the random effects and using that to estimate the variance of the
variance estimate in question. This is the function `gls_grad`. For
`lme4::mermod` and [`nlme::lme`](https://rdrr.io/pkg/nlme/man/lme.html)
objects, we don’t have that variance matrix, so instead we have an
approximate method (`"appx-satterthwaite"`) that perturbs the response
values slightly and refits the model. It only takes a few of these
perturbations to do a respectable job of estimating the required
variances; this is done by the `gradV.kludge` function.

### Estimability

An important part of the package is that we assess estimability of our
predictions. This is important because some models allow rank-deficient
model matrices; if we have a rank-deficient model, there are infinitely
many possible solutions for the fixed effects, but we have only one of
them. A prediction is defined as estimable if it is the same, no matter
which solution we used. Equivalently, given the model matrix \\X\\,
\\x'\beta\\ is estimable if, and only if, \\x\\ is in the row space of
\\X\\.

The **estimability** package provides the functions needed to assess
estimability. We do this by creating a basis \\N\\ for the null space of
\\X\\, i.e., \\XN=0\\, and so \\x'\beta\\ is estimable iff \\x'N=0\\. We
store the matrix \\N\\ in the `@nbasis` slot of an `emmGrid` object.

For more details, see the package documentation for **estimability** and
its
[vignette](https://cran.r-project.org/web/packages/estimability/vignettes/add-est-check.html).
