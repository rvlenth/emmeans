R package **emmeans**: Estimated marginal means (least-squares means)
====

#### Note
**emmeans** is a continuation of the package **lsmeans**. The latter will eventually be retired.

## Features
Estimated marginal means (EMMs, previously known as least-squares means in the context of traditional regression models) are derived by using a model to make predictions over a regular grid of pridictor combinations (called a *reference grid*). These predictions may possibly be averaged (typically with equal weights) over one or more of the predictors. Such marginally-averaged predictions are useful for describing the results of fitting a model, particularly in presenting the effects of factors. The **emmeans** package can easily produce these results, as well as various graphs of them (interaction-style plots and side-by-side intervals).
* Estimation and testing of pairwise comparisons of EMMs, and several other types of contrasts, are provided. There is also a `cld` method for display of grouping symbols.
* Two-way support of the `glht` function in the **multcomp** package.
* For models where continuous predictors interact with factors, the package's `emtrends` function works in terms of a reference grid of predicted slopes of trend lines for each factor combination.
* Incorporates support for many types of models, including those in **stats** package (`lm`, `glm`, `aov`, `aovlist`), linear and genearalized linear mixed models (e.g., **nlme**, **lme4**, **afex**), ordinal-response models (e.g., **ordinal**, **MASS**), survival analysis (e.g., **survival**, **coxme**), generalized estimating equations (**gee**, **geepack**), and others. See `help("models", package = "emmeans")`
* Various Bayesian models (**carBayes**, **MCMCglmm**, **MCMCpack**) are supported by way of creating a posterior sample of least-squares means or contrasts thereof, which may then be examined using tools such as in the **coda** package.
* Package developers may provide **emmeans** support for their models by providing `recover_data` and `emm_basis` methods. See `vignette("extending", package = "emmeans")`


* To install the latest development version from Github, have the newest **devtools** package installed, then run
```
devtools::install_github("rvlenth/emmeans", dependencies = TRUE,
                        build_vignettes = TRUE)
```
For latest release notes on this development version, see the [NEWS file](https://github.com/rvlenth/emmeans/blob/master/inst/NEWS)
