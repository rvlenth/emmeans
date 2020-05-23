R package **emmeans**: Estimated marginal means
====

##### *Note: **emmeans** is a continuation of the package **lsmeans**. The latter will eventually be retired.*

## Features
Estimated marginal means (EMMs, previously known as least-squares means in the
context of traditional regression models) are derived by using a model to make
predictions over a regular grid of predictor combinations (called a *reference
grid*). These predictions may possibly be averaged (typically with equal
weights) over one or more of the predictors. Such marginally-averaged
predictions are useful for describing the results of fitting a model,
particularly in presenting the effects of factors. The **emmeans** package can
easily produce these results, as well as various graphs of them
(interaction-style plots and side-by-side intervals).


  * Estimation and testing of pairwise comparisons of EMMs, and several other
    types of contrasts, are provided. There is also a `cld` method for display of
    grouping  symbols.
    
  * Two-way support of the `glht` function in the **multcomp** package.
  
  * For models where continuous predictors interact with factors, the package's
    `emtrends` function works in terms of a reference grid of predicted slopes of
    trend lines for each factor combination.
    
  * Vignettes are provided on various aspects of EMMs and using the package. 
    See the [CRAN page](https://CRAN.R-project.org/package=emmeans)


## Model support


  * The package incorporates support for many types of models, including 
    standard models fitted using `lm`, `glm`, and relatives, 
    various mixed models, GEEs, survival models, count models,
    ordinal responses, zero-inflated models, and others. Provisions for
    some models include special modes for accessing different types of 
    predictions; for example, with zero-inflated models, one may opt for
    the estimated response including zeros, just the linear predictor, 
    or the zero model.
    For details, see
    [`vignette("models", package = "emmeans")`](https://CRAN.R-project.org/package=emmeans/vignettes/models.html)
    
  * Various Bayesian models (**carBayes**, **MCMCglmm**, **MCMCpack**) are
    supported by way of creating a posterior sample of least-squares means or
    contrasts thereof, which may then be examined using tools such as in the
    **coda** package.
    
  * Package developers may provide **emmeans** support for their models by
    writing `recover_data` and `emm_basis` methods. See `vignette("extending",
    package = "emmeans")`
    

## Versions and installation


  * **CRAN** The latest CRAN version may be found at [https://CRAN.R-project.org/package=emmeans](https://CRAN.R-project.org/package=emmeans).
    Also at that site, formatted versions of this package's vignettes 
    may be viewed.

  * **Github** To install the latest development version from Github, 
    install the newest version (definitely 2.0 or higher) of the **devtools** 
    package; then run
    
```r
remotes::install_github("rvlenth/emmeans", dependencies = TRUE, build_opts = "")

### To install without vignettes (faster):
remotes::install_github("rvlenth/emmeans")
```
*Note:* If you are a Windows user, you should also first download and
      install the latest version of
      [`Rtools`](https://cran.r-project.org/bin/windows/Rtools/).

For the latest release notes on this development version, see the 
[NEWS file](https://github.com/rvlenth/emmeans/blob/master/NEWS.md)
