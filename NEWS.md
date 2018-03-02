emmeans 1.1.2.9001
------------------

  * Formatting corrections in documentation
  

emmeans 1.1.2
------------------
 
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
    

