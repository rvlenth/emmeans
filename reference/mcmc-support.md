# Support for MCMC-based estimation

When a model is fitted using Markov chain Monte Carlo (MCMC) methods,
its reference grid contains a `post.beta` slot. These functions
transform those posterior samples to posterior samples of EMMs or
related contrasts. They can then be summarized or plotted using, e.g.,
functions in the coda package.

## Usage

``` r
# S3 method for class 'emmGrid'
as.mcmc(x, names = TRUE, sep.chains = TRUE, likelihood,
  NE.include = FALSE, ...)

# S3 method for class 'emm_list'
as.mcmc(x, which = 1, ...)

# S3 method for class 'emmGrid'
as.mcmc.list(x, names = TRUE, ...)

# S3 method for class 'emm_list'
as.mcmc.list(x, which = 1, ...)
```

## Arguments

- x:

  An object of class `emmGrid`

- names:

  Logical scalar or vector specifying whether variable names are
  appended to levels in the column labels for the `as.mcmc` or
  `as.mcmc.list` result – e.g., column names of `treat A` and `treat B`
  versus just `A` and `B`. When there is more than one variable
  involved, the elements of `names` are used cyclically.

- sep.chains:

  Logical value. If `TRUE`, and there is more than one MCMC chain
  available, an
  [`mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html) object is
  returned by `as.mcmc`, with separate EMMs posteriors in each chain.

- likelihood:

  Character value or function. If given, simulations are made from the
  corresponding posterior predictive distribution. If not given, we
  obtain the posterior distribution of the parameters in `object`. See
  Prediction section below.

- NE.include:

  Logical value. If `TRUE`, non-estimable columns are kept but returned
  as columns of `NA` values (this may create errors or warnings in
  subsequent analyses using, say, coda). If `FALSE`, non-estimable
  columns are dropped, and a warning is issued. (If all are
  non-estimable, an error is thrown.)

- ...:

  arguments passed to other methods

- which:

  item in the `emm_list` to use

## Value

An object of class [`mcmc`](https://rdrr.io/pkg/coda/man/mcmc.html) or
[`mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html).

## Details

When the object's `post.beta` slot is non-trivial, `as.mcmc` will return
an [`mcmc`](https://rdrr.io/pkg/coda/man/mcmc.html) or
[`mcmc.list`](https://rdrr.io/pkg/coda/man/mcmc.list.html) object that
can be summarized or plotted using methods in the coda package. In these
functions, `post.beta` is transformed by post-multiplying it by
`t(linfct)`, creating a sample from the posterior distribution of LS
means. In `as.mcmc`, if `sep.chains` is `TRUE` and there is in fact more
than one chain, an `mcmc.list` is returned with each chain's results.
The `as.mcmc.list` method is guaranteed to return an `mcmc.list`, even
if it comprises just one chain.

## Prediction

When `likelihood` is specified, it is used to simulate values from the
posterior predictive distribution corresponding to the given likelihood
and the posterior distribution of parameter values. Denote the
likelihood function as \\f(y\|\theta,\phi)\\, where \\y\\ is a response,
\\\theta\\ is the parameter estimated in `object`, and \\\phi\\
comprises zero or more additional parameters to be specified. If
`likelihood` is a function, that function should take as its first
argument a vector of \\\theta\\ values (each corresponding to one row of
`object@grid`). Any \\\phi\\ values should be specified as additional
named function arguments, and passed to `likelihood` via `...`. This
function should simulate values of \\y\\.

A few standard likelihoods are available by specifying `likelihood` as a
character value. They are:

- `"normal"`:

  The normal distribution with mean \\\theta\\ and standard deviation
  specified by additional argument `sigma`

- `"binomial"`:

  The binomial distribution with success probability \\theta\\, and
  number of trials specified by `trials`

- `"poisson"`:

  The Poisson distribution with mean \\theta\\ (no additional
  parameters)

- `"gamma"`:

  The gamma distribution with scale parameter \\\theta\\ and shape
  parameter specified by `shape`

## Examples

``` r
if(requireNamespace("coda")) 
    emm_example("as.mcmc-coda")
#> 
#> --- Running code from 'system.file("extexamples", "as.mcmc-coda.R", package = "emmeans")'
#> 
#> > cbpp.rg <- do.call(emmobj, readRDS(system.file("extdata", 
#> +     "cbpplist", package = "emmeans")))
#> 
#> > pred.incidence <- coda::as.mcmc(regrid(cbpp.rg), likelihood = "binomial", 
#> +     trials = 20)
#> 
#> > summary(pred.incidence)
#> 
#> Iterations = 1:250
#> Thinning interval = 1 
#> Number of chains = 2 
#> Sample size per chain = 250 
#> 
#> 1. Empirical mean and standard deviation for each variable,
#>    plus standard error of the mean:
#> 
#>                                 Mean    SD Naive SE Time-series SE
#> size 15.0357142857143 period 1 3.860 1.912  0.08551        0.10055
#> size 15.0357142857143 period 2 1.652 1.391  0.06221        0.06134
#> size 15.0357142857143 period 3 1.598 1.250  0.05592        0.05594
#> size 15.0357142857143 period 4 1.058 1.062  0.04752        0.05217
#> 
#> 2. Quantiles for each variable:
#> 
#>                                2.5% 25% 50% 75% 97.5%
#> size 15.0357142857143 period 1    1   3   4   5     8
#> size 15.0357142857143 period 2    0   1   1   2     5
#> size 15.0357142857143 period 3    0   1   1   2     5
#> size 15.0357142857143 period 4    0   0   1   2     4
#> 
#> 
    # Use emm_example("as.mcmc-coda", list = TRUE) # to see just the code
    
```
