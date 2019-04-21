## NEWS for the emmeans package

emmeans 1.3.4
-------------

  * Un-naming a lot of arguments in `do.call(paste, ...)` and `do.call(order, ...)`,
    to prevent problems with factor names like `method` that are argument names
    for these functions (#94)
  * Fix to a logic error in `summary.emmGrid()` whereby transformations of class
    `list` were ignored.
  * Enhancement to `update.emmGrid(..., levels = levs)` whereby we can easily
    relabel the reference grid and ensure that the `grid` and `roles` slots
    stay consistent. Added vignette example.
  * Clarified ordering rules used by `emmeans()`. We now ensure that the
    original order of the reference grid is preserved. Previously, the grid 
    was re-ordered if any numeric or character levels occurred out of order, 
    per `order()`
  * Curbing use of "statistical significance" language. This includes
    additional vignette material and plans to deprecate `CLD()` due to its 
    misleading display of pairwise-comparison tests.
  * Bug fix for `betareg` objects, where the wrong `terms` component was 
    sometimes used.
  * Correction to logic error that affected multiplicity adjustments when
    `by` variables are present (#98).
  * Addition of `pwpp()` function to plot *P* values of comparisons
  * Improvement to `summary(..., adjust = "scheffe")`. We now actually
    compute and use the rank of the matrix of linear functions to obtain
    the *F* numerator d.f., rather than trying to guess the likely correct 
    value.
  * Removal of vignette on transitioning from **lsmeans** -- 
    it's been a long enough time now.
  

emmeans 1.3.3
-------------

  * Fix to unintended consequence of #71 that caused incorrect ordering 
    of `contrast()` results if they are later used by `emmeans()`.
    This was first noticed with ordinal models in `prob` mode (#83).
  * Improved checking of conformability of parameters -- for models
    with rank deficiency not handled same way as lm()'s NA convention
  * Added basic support for `sommer::mmer`, `MuMIn::averaging`, and
    `mice::mira` objects
  * Fix in `nnet::multinom` support when there are 2 outcomes (#19)
  * Added Satterthwaite d.f. to `gls` objects
  * `famSize` now correct when `exclude` or `include` is used in 
    a contrast function (see #68)
  * Stronger warnings of possible bias with `aovList` objects, in part
    due to the popularity of `afex::aov_ez()` which uses these models.
  * Updates to FAQs vignette



emmeans 1.3.2
-------------

  * I decided to enable "optimal digits" display by default. In summaries,
    we try to show enough---but not too much---precision in estimates and
    confidence intervals. If you don't like this and want to revert
    to the old (exaggerated precision) behavior, do 
    `emm_options(opt.digits = FALSE)`
  * Added `include` argument to most `.emmc` functions (#67)
  * Now allow character values for `ref`, `exclude`, and `include` in
    `.emmc` functions (#68)
  * Better handling of matrix predictors (#66)
  * Fixed over-zealous choice to not pass `...` arguments in `emmeans()`
    when two-sided formulas are present
  * Fix to `clm` support when model is rank-deficient
  * Fix to `regrid(..., transform = "log")` error when there are
    existing non-estimable cases (issue #65)
  * Improvements to `brmsfit` support (#43)
  * Added support for `mgcv::gam` and `mgcv::gamm` models
  * `.my.vcov()` now passes `...` to clients
  * Removed **glmmADMB** support. This package appears to be dormant
  * Fixed ordering bug for nested models (#71)
  * Support for `manova` object no longer requires `data` keyword (#72)
  * Added support for multivariate response in `aovlist` models (#73)
  * Documentation clarification (#76)
  * Fix to `CLD` fatal error when `sort = TRUE` (#77)
  * Fix to issue with weights and incomplete cases with `lme` objects (#75)
  * Nested fixed-effects yielded NonEsts when two factors are nested 
    in the same factor(s) (#79)


emmeans 1.3.1
-------------

  * `"mvt"` adjustment ignored `by` grouping
  * `contrast()` mis-labeled estimates when levels varied among `by` groups
    (most prominently this happened in `CLD(..., details = TRUE)`)
  * Changed `aovlist` support so it re-fits the model when non-sum-to-zero
    contrasts were used
  * `print.summary_emm()` now cleans up numeric columns with `zapsmall()`
  * More robust handling of `nesting` in `ref_grid()` and `update()`,
    and addition of `covnest` argument for whether to include covariates
    when auto-detecting nesting
  * Revision of some vignettes
  * Fixed bug in `hpd.summary()` and handoff to it from `summary()`
  * Fixed bug where `ref_grid()` ignored `mult.levs`
  * Fixes in emmeans where it passes `...` where it shouldn't
  * `CLD()` now works for MCMC models (uses frequentist summary)
  * Addition of `opt.digits` option


emmeans 1.3.0
-------------

  * Deprecated functions like `ref.grid()` put to final rest, and we no 
    longer support packages that provide `recover.data` or `lsm.basis` methods
  * Courtesy exports `.recover_data()` and `.emm_basis()` to provide
    access for extension developers to all available methods
  * Streamlining of a stored example in `inst/extdata`
  * Fix to `.all.vars()` that could cause errors when response variable
    has a function call with character constants.
  * Relabeling of differences as ratios when appropriate in `regrid()`
    (so results match `summary()` labeling with `type = "response"`).
  * `plot.emmGrid(..., comparisons = TRUE, type = "response")`
    produced incorrect comparison arrows; now fixed


emmeans 1.2.4
-------------

  * Support for model formulas such as `df$y ~ df$treat + df[["cov"]]`. This had
    failed previously for two obscure reasons, but now works correctly.
  * New `simplify.names` option for above types of models
  * `emm_options()` with no arguments now returns all options in force,
    including the defaults. This makes it more consistent with `options()`
  * Bug fix for `emtrends()`; produced incorrect results in models with offsets. 
  * Separated the help pages for `update.emmGrid()` and `emm_options()`
  * New `qdrg()` function (quick and dirty reference grid) for help with
    unsupported model objects


emmeans 1.2.3
-------------

  * S3 methods involving packages **multcomp** and **coda** are now
    dynamically registered, not merely exported as functions.
    This passes checks when S3 methods are required to be registered.
  * `cld()` has been deprecated in favor of `CLD()`. This had been a
    headache. **multcomp** is the wrong place for the generic to be; 
    it is too fancy a dance to export `cld` with or without having
    **multcomp** installed.
  * Added vignette caution regarding interdependent covariates
  * Improved **glmmADMB** support to recover contrasts correctly
  

emmeans 1.2.2
-------------

  * Removed ggplot2, multcomp, and coda to Suggests -- thus vastly
    reducing dependencies
  * Added a FAQ to the FAQs vignette
  * Modified advice in `xtending.Rmd` vignette on how to export methods
  * Fixes to `revpairwise.emmc` and `cld` regarding comparing only 1 EMM
  * `cld.emm_list` now returns results only for `object[[ which[1] ]]`,
    along with a warning message.
  * Deprecated `emmeans` specs like `cld ~ group`, a vestige of **lsmeans**
    as it did not work correctly (and was already undocumented)


emmeans 1.2.1
-------------

  * Moved **brms** to `Suggests` (dozens and dozens fewer dependencies)


emmeans 1.2
-----------

  * Index of vignette topics added
  * New, improved (to my taste) vignette formats
  * Fixed df bug in regrid (#29)
  * Fixed annotation bug for nested models (#30)
  * Better documentation for `lme` models in "models" vignette
  * Additional fixes for arguments passed to `.emmc` functions (#22)
  * Support added for logical predictors (who knew we could have those? not me)
  * Replaced tex/pdf "Extending" vignette with Rmd/html
  * Overhauled the faulty logic for df methods in emm_basis.merMod
  * Added Henrik to contributors list (long-standing oversight)
  * Added `exclude` argument to most `.emmc` functions: allows
    user to omit certain levels when computing contrasts
  * New `hpd.summary()` function for Bayesian models to show HPD intervals
    rather than frequentist summary. Note: `summary()` automatically
    reroutes to it. Also `plot()` and `emmip()` play along.
  * Rudimentary support for **brms** package
  * *Ad hoc* Satterthwaite method for `nlme::lme` models


emmeans 1.1.3
-------------

  * Formatting corrections in documentation
  * Fixed bug for survival models where `Surv()` was interpreted
    as a response transformation.
  * Fixed bug (issue #19) in multinom support
  * Fixed bug (issue #22) in optional arguments with 
    interaction contrasts
  * Fixed bug (issue #23) in weighting with character predictors
  * Clarifying message when `cld()` is applied to an `emm_list` (issue #24)
  * Added `offset` argument to `ref_grid()` (scalar offset only) and to
    `emmeans()` (vector offset allowed) -- (issue #18)
  * New optional argument for `[.summary_emm` to choose whether to 
    retain its class or coerce to a `data.frame` (relates to issue #14)
  * Added `reverse` option for `trt.vs.ctrl` and relatives (#27)
  

emmeans 1.1.2
-------------
 
  * Changed the way `terms` is accessed with `lme` objects to make
    it more robust
  * `emmeans:::convert_scripts()` renames output file more simply
  * Added `[` method for class `summary_emm`
  * Added `simple` argument for `contrast` - essentially the complement of `by`
  * Improved estimability handling in `joint_tests()`
  * Made `ref_grid()` accept `ylevs` list of length > 1; 
    also slight argument change: `mult.name` -> `mult.names`
  * Various bug fixes, bullet-proofing
  * Fixes to make Markdown files render better


emmeans 1.1
-----------

  * Fixed a bug in `emmeans()` wherein `weights` was ignored 
    when `specs` is a `list`
  * Coerce `data` argument, if supplied to a data.frame 
    (`recover_data()` doesn't like tibbles...)
  * Added `as.data.frame` method for `emmGrid` objects, making it
    often possible to pass it directly to other functions as a `data` 
    argument.
  * Fixed bug in `contrast()` where `by` was ignored for 
    interaction contrasts
  * Fixed bug in `as.glht()` where it choked on `df = Inf`
  * Fixed bug occurring when a model call has no `data` or `subset`
  * New `joint_tests()` function tests all [interaction] contrasts


emmeans 1.0
-----------

  * Added preliminary support for `gamlss` objects (but doesn't support
    smoothing). Additional argument is `what = c("mu", "sigma", "nu", "tau")`
    It seems to be flaky when the model of interest is just `~ 1`.
  * Improved support for models with fancy variable names 
    (containing spaces and such)
  * Fixed a bug whereby `emmeans()` might pass `data` to `contrast()`
  * Added some missing documentation for `summary.emmGrid()`
  * Repaired handling of `emm_options(summary = ...)` to work as
    advertised. 
  * Changed many object names in examples and vignettes from xxx.emmGrid
    to xxx.emm (result of overdoing the renaming the object class itself)
  * Changed `emmGrid()` function to `emm()` as had been intended
    as alternative to `mcp()` in `multcomp::glht()` (result of ditto).
  * Fixed error in exporting `cld.emm_list()`
  * Fixed a bug whereby all CIs were computed using the first estimate's 
    degrees of freedom.
  * Now using `Inf` to display d.f. for asymptotic (z) tests. (`NA` will
    still work too but `Inf` is a better choice for consistency and meaning.)
  * Bug fix in nesting-detection code when model has only an intercept
    

emmeans 0.9.1
-------------

  * Documentation corrections (broken links, misspellings, mistakes)
  * More sophisticated check for randomized data in `recover_data()`
    now throws an error when it finds recovered data not reproducible
  * Added support for gam::gam objects
  * Fixes to `vcov()` calls to comply with recent R-devel changes
  

emmeans 0.9
-----------

This is the initial major version that replaces the **lsmeans** package.
Changes shown below are changes made to the last real release of **lsmeans**
(version 2.27-2). **lsmeans** versions greater than that are transitional
to that package being retired.

  * We now emphasize the terminology "estimated marginal means" rather 
    than "least-squares means"
  * The flagship functions are now `emmeans()`, `emtrends()`, `emmip()`, etc.
    But `lsmeans()`, `lstrends()`, etc. as well as `pmmeans()` etc. are
    mapped to their corresponding `emxxxx()` functions.
  * In addition, we are trying to avoid names that could get confused as
    S3 methods. So, `ref.grid -> ref_grid`, `lsm.options -> emm_options`, etc.
  * Classes `ref.grid` and `lsmobj` are gone.
    Both are replaced by class `emmGrid`. An `as.emmGrid()` function
    is provided to convert old objects to class `emmGrid`.
  * I decided to revert back to "kenward-roger" as the default degrees-of-freedom
    method for `lmerMod models`. Also added options `disable.lmerTest`
    and `lmerTest.limit`, similar to those for **pbkrtest**.
  * Documentation and NAMESPACE are now "ROxygenated"
  * Additional `neuralgia` and `pigs` datasets
  * Dispatching of `emmmeans()` methods is now top-down rather than convoluted
    intermingling of S3 methods
  * Improved display of back-transformed contrasts when log or logit 
    transformation was used: We change any ` - `s in labels to ` / `s
    to emphasize that thnese results are ratios.
  * A message is now displayed when nesting is auto-detected in `ref_grid`.
    (Can be disabled via `emm_options()`)
  * Options were added for several messages that users may want to suppress,
    e.g., ones about interactions and nesting.
  * Greatly overhauled help page for models. It is now a vignette, 
    with a quick reference chart linked to details, and is
    organized by similarities instead of packages.
  * Support for 'mer' objects (lme4.0 package) removed.
  * A large number of smaller interlinked vignettes replaces the one
    big one on using the package. Several vignettes are linked in the
    help pages.
  * Graphics methods `plot()` and `emmip()` are now **ggplot2**-based.
    Old **lattice**-based functionality is still available too,
    and there is a `graphics.engine` option to choose the default.
  * Non-exported utilities convert_workspace() and convert_scripts() to
    help with transition
  * Moved `Suggests` pkgs to `Enhances` when not needed for 
    building/testing

### NOTE: **emmeans** is a continuation of the **lsmeans** package. 
New developments will take place in **emmeans**, and **lsmeans** 
will remain static and eventually will be archived.
    

