# Summarize an emmGrid from a Bayesian model

This function computes point estimates and HPD (or quantile) intervals
for each factor combination in `object@emmGrid`. While this function may
be called independently, it is called automatically by the S3 method
[`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
when the object is based on a Bayesian model. (Note: the `level`
argument, or its default, is passed as `prob`).

## Usage

``` r
hpd.summary(object, prob, by, type, point.est = median,
  ci.method = get_emm_option("post.ci.method"), delta,
  bias.adjust = get_emm_option("back.bias.adj"), sigma, ...)
```

## Arguments

- object:

  an `emmGrid` object having a non-missing `post.beta` slot

- prob:

  numeric probability content for HPD intervals (note: when not
  specified, the current `level` option is used; see
  [`emm_options`](https://rvlenth.github.io/emmeans/reference/emm_options.md))

- by:

  factors to use as `by` variables

- type:

  prediction type as in
  [`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)

- point.est:

  function to use to compute the point estimates from the posterior
  sample for each grid point

- ci.method:

  character value matching `"HPD"` (default) or `"quantile"` (but
  actually not case-sensitive). The default is to use HPD intervals
  (which are the shortest possible intervals). Alternatively, choosing
  `"quantile"` uses the quantiles of the posterior having equal tail
  probabilities.

- delta:

  Numeric equivalence threshold (on the linear predictor scale
  regardless of `type`). See the section below on equivalence testing.

- bias.adjust:

  Logical value for whether to adjust for bias in back-transforming
  (`type = "response"`). This requires a value of `sigma` to exist in
  the object or be specified.

- sigma:

  Error SD assumed for bias correction (when `type = "response"`. If not
  specified, `object@misc$sigma` is used, and a warning if it is not
  found or invalid. *Note:* `sigma` may be a vector, as long as it
  conforms to the number of observations in the posterior sample.

- ...:

  required but not used

## Value

an object of class `summary_emm`

## Note

HPD intervals require the coda package to be installed on your system;
otherwise an error is thrown. (So one way to sidestep that error is to
specify `ci.method = "quantile"`).

## Equivalence testing note

If `delta` is positive, two columns labeled `p.equiv` and `odds.eq` are
appended to the summary. `p.equiv` is the fraction of posterior
estimates having absolute values less than `delta`. The `odds.eq` column
is just `p.equiv` converted to an odds ratio; so it is the posterior
odds of equivalence.

A high value of `p.equiv` is evidence in favor of equivalence. It can be
used to obtain something equivalent (in spirit) to the frequentist
Schuirmann (TOST) procedure, whereby we would conclude equivalence at
significance level \\\alpha\\ if the \\(1 - 2\alpha)\\ confidence
interval falls entirely in the interval \\\[-\delta, \delta\]\\.
Similarly in the Bayesian context, an equally strong argument for
equivalence is obtained if `p.equiv` exceeds \\1 - 2\alpha\\.

A closely related quantity is the ROPE (region of practical
equivalence), obtainable via
`bayestestR::rope(object, range = c(-delta, delta))`. Its value is
approximately `100 * p.equiv / 0.95` if the default `ci = 0.95` is used.
See also bayestestR's [issue
\#567](https://github.com/easystats/bayestestR/issues/567).

Finally, a Bayes factor for equivalence is obtainable by dividing
`odds.eq` by the prior odds of equivalence, assessed or elicited
separately.

## See also

summary.emmGrid

## Examples

``` r
if(require("coda")) 
    emm_example("hpd.summary-coda")
#> Loading required package: coda
#> 
#> --- Running code from 'system.file("extexamples", "hpd.summary-coda.R", package = "emmeans")'
#> 
#> > cbpp.rg <- do.call(emmobj, readRDS(system.file("extdata", 
#> +     "cbpplist", package = "emmeans")))
#> 
#> > cbpp.emm <- emmeans(cbpp.rg, "period")
#> 
#> > hpd.summary(cbpp.emm)
#>  period emmean lower.HPD upper.HPD
#>  1       -1.43     -1.96    -0.894
#>  2       -2.39     -3.13    -1.823
#>  3       -2.52     -3.19    -1.862
#>  4       -2.97     -3.88    -1.998
#> 
#> Point estimate displayed: median 
#> Results are given on the logit (not the response) scale. 
#> HPD interval probability: 0.95 
#> 
#> > summary(pairs(cbpp.emm), type = "response", delta = log(2))
#>  contrast          odds.ratio lower.HPD upper.HPD p.equiv odds.eq
#>  period1 / period2       2.60     1.264      4.34   0.210  0.2658
#>  period1 / period3       2.90     1.266      5.42   0.114  0.1287
#>  period1 / period4       4.64     1.595      9.76   0.026  0.0267
#>  period2 / period3       1.14     0.486      2.36   0.918 11.1951
#>  period2 / period4       1.78     0.580      4.35   0.604  1.5253
#>  period3 / period4       1.55     0.437      4.01   0.700  2.3333
#> 
#> Point estimate displayed: median 
#> 'p.equiv' and 'odds.eq' based on posterior P(|lin. pred.| < 0.6931) 
#> Results are back-transformed from the log odds ratio scale 
#> HPD interval probability: 0.95 
#> 
    # Use emm_example("hpd.summary-coda", list = TRUE) # to see just the code
```
