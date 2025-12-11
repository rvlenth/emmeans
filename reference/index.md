# Package index

## Overview

- [`emmeans-package`](https://rvlenth.github.io/emmeans/reference/emmeans-package.md)
  : Estimated marginal means (aka Least-squares means)

## Estimation of marginal means and trends

- [`ref_grid()`](https://rvlenth.github.io/emmeans/reference/ref_grid.md)
  : Create a reference grid from a fitted model
- [`emmeans()`](https://rvlenth.github.io/emmeans/reference/emmeans.md)
  : Estimated marginal means (Least-squares means)
- [`emtrends()`](https://rvlenth.github.io/emmeans/reference/emtrends.md)
  : Estimated marginal means of linear trends
- [`lsmeans()`](https://rvlenth.github.io/emmeans/reference/wrappers.md)
  [`lstrends()`](https://rvlenth.github.io/emmeans/reference/wrappers.md)
  [`lsmip()`](https://rvlenth.github.io/emmeans/reference/wrappers.md)
  [`lsm()`](https://rvlenth.github.io/emmeans/reference/wrappers.md)
  [`lsmobj()`](https://rvlenth.github.io/emmeans/reference/wrappers.md)
  [`lsm.options()`](https://rvlenth.github.io/emmeans/reference/wrappers.md)
  [`get.lsm.option()`](https://rvlenth.github.io/emmeans/reference/wrappers.md)
  : Wrappers for alternative naming of EMMs
- [`qdrg()`](https://rvlenth.github.io/emmeans/reference/qdrg.md) :
  Quick and dirty reference grid

## Follow-up analyses

- [`contrast()`](https://rvlenth.github.io/emmeans/reference/contrast.md)
  [`pairs(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/contrast.md)
  [`coef(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/contrast.md)
  [`weights(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/contrast.md)
  : Contrasts and linear functions of EMMs
- [`eff_size()`](https://rvlenth.github.io/emmeans/reference/eff_size.md)
  : Calculate Cohen effect sizes and confidence bounds thereof
- [`joint_tests()`](https://rvlenth.github.io/emmeans/reference/joint_tests.md)
  [`make.meanint()`](https://rvlenth.github.io/emmeans/reference/joint_tests.md)
  [`meanint()`](https://rvlenth.github.io/emmeans/reference/joint_tests.md)
  [`make.symmint()`](https://rvlenth.github.io/emmeans/reference/joint_tests.md)
  [`symmint()`](https://rvlenth.github.io/emmeans/reference/joint_tests.md)
  : Compute joint tests of the terms in a model
- [`mvcontrast()`](https://rvlenth.github.io/emmeans/reference/mvcontrast.md)
  : Multivariate contrasts
- [`pairwise.emmc()`](https://rvlenth.github.io/emmeans/reference/emmc-functions.md)
  [`revpairwise.emmc()`](https://rvlenth.github.io/emmeans/reference/emmc-functions.md)
  [`tukey.emmc()`](https://rvlenth.github.io/emmeans/reference/emmc-functions.md)
  [`poly.emmc()`](https://rvlenth.github.io/emmeans/reference/emmc-functions.md)
  [`opoly.emmc()`](https://rvlenth.github.io/emmeans/reference/emmc-functions.md)
  [`trt.vs.ctrl.emmc()`](https://rvlenth.github.io/emmeans/reference/emmc-functions.md)
  [`trt.vs.ctrl1.emmc()`](https://rvlenth.github.io/emmeans/reference/emmc-functions.md)
  [`trt.vs.ctrlk.emmc()`](https://rvlenth.github.io/emmeans/reference/emmc-functions.md)
  [`dunnett.emmc()`](https://rvlenth.github.io/emmeans/reference/emmc-functions.md)
  [`eff.emmc()`](https://rvlenth.github.io/emmeans/reference/emmc-functions.md)
  [`del.eff.emmc()`](https://rvlenth.github.io/emmeans/reference/emmc-functions.md)
  [`consec.emmc()`](https://rvlenth.github.io/emmeans/reference/emmc-functions.md)
  [`mean_chg.emmc()`](https://rvlenth.github.io/emmeans/reference/emmc-functions.md)
  [`helmert.emmc()`](https://rvlenth.github.io/emmeans/reference/emmc-functions.md)
  [`nrmlz.emmc()`](https://rvlenth.github.io/emmeans/reference/emmc-functions.md)
  [`wtcon.emmc()`](https://rvlenth.github.io/emmeans/reference/emmc-functions.md)
  [`identity.emmc()`](https://rvlenth.github.io/emmeans/reference/emmc-functions.md)
  : Contrast families

## Summaries

- [`summary(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
  [`confint(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
  [`test()`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
  [`predict(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
  [`as.data.frame(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
  [`` `[`( ``*`<summary_emm>`*`)`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
  :

  Summaries, predictions, intervals, and tests for `emmGrid` objects

- [`hpd.summary()`](https://rvlenth.github.io/emmeans/reference/hpd.summary.md)
  : Summarize an emmGrid from a Bayesian model

- [`untidy`](https://rvlenth.github.io/emmeans/reference/untidy.md) :
  Dare to be un-"tidy"!

## Modifying objects

- [`comb_facs()`](https://rvlenth.github.io/emmeans/reference/manip-factors.md)
  [`split_fac()`](https://rvlenth.github.io/emmeans/reference/manip-factors.md)
  [`add_grouping()`](https://rvlenth.github.io/emmeans/reference/manip-factors.md)
  [`add_submodels()`](https://rvlenth.github.io/emmeans/reference/manip-factors.md)
  [`permute_levels()`](https://rvlenth.github.io/emmeans/reference/manip-factors.md)
  : Manipulate factors in a reference grid

- [`rbind(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/rbind.emmGrid.md)
  [`` `+`( ``*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/rbind.emmGrid.md)
  [`` `[`( ``*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/rbind.emmGrid.md)
  [`head(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/rbind.emmGrid.md)
  [`tail(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/rbind.emmGrid.md)
  [`subset(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/rbind.emmGrid.md)
  [`rbind(`*`<emm_list>`*`)`](https://rvlenth.github.io/emmeans/reference/rbind.emmGrid.md)
  [`rbind(`*`<summary_emm>`*`)`](https://rvlenth.github.io/emmeans/reference/rbind.emmGrid.md)
  [`force_regular()`](https://rvlenth.github.io/emmeans/reference/rbind.emmGrid.md)
  :

  Combine or subset `emmGrid` objects

- [`regrid()`](https://rvlenth.github.io/emmeans/reference/regrid.md) :
  Reconstruct a reference grid with a new transformation or simulations

- [`mvregrid()`](https://rvlenth.github.io/emmeans/reference/mvregrid.md)
  : Multivariate regridding

- [`update(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md)
  [`` `levels<-`( ``*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md)
  [`update(`*`<summary_emm>`*`)`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md)
  :

  Update an `emmGrid` object

## Displays

- [`cld(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/CLD.emmGrid.md)
  [`cld(`*`<emm_list>`*`)`](https://rvlenth.github.io/emmeans/reference/CLD.emmGrid.md)
  : Compact letter displays

- [`emmip()`](https://rvlenth.github.io/emmeans/reference/emmip.md)
  [`emmip_ggplot()`](https://rvlenth.github.io/emmeans/reference/emmip.md)
  [`emmip_lattice()`](https://rvlenth.github.io/emmeans/reference/emmip.md)
  : Interaction-style plots for estimated marginal means

- [`plot(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/plot.md)
  [`plot(`*`<summary_emm>`*`)`](https://rvlenth.github.io/emmeans/reference/plot.md)
  :

  Plot an `emmGrid` or `summary_emm` object

- [`pwpm()`](https://rvlenth.github.io/emmeans/reference/pwpm.md) :
  Pairwise P-value matrix (plus other statistics)

- [`pwpp()`](https://rvlenth.github.io/emmeans/reference/pwpp.md) :
  Pairwise P-value plot

## Data sets

- [`auto.noise`](https://rvlenth.github.io/emmeans/reference/auto.noise.md)
  : Auto Pollution Filter Noise
- [`feedlot`](https://rvlenth.github.io/emmeans/reference/feedlot.md) :
  Feedlot data
- [`fiber`](https://rvlenth.github.io/emmeans/reference/fiber.md) :
  Fiber data
- [`MOats`](https://rvlenth.github.io/emmeans/reference/MOats.md) : Oats
  data in multivariate form
- [`neuralgia`](https://rvlenth.github.io/emmeans/reference/neuralgia.md)
  : Neuralgia data
- [`nutrition`](https://rvlenth.github.io/emmeans/reference/nutrition.md)
  : Nutrition data
- [`oranges`](https://rvlenth.github.io/emmeans/reference/oranges.md) :
  Sales of oranges
- [`pigs`](https://rvlenth.github.io/emmeans/reference/pigs.md) :
  Effects of dietary protein on free plasma leucine concentration in
  pigs
- [`ubds`](https://rvlenth.github.io/emmeans/reference/ubds.md) :
  Unbalanced dataset

## Structures

- [`emmGrid-class`](https://rvlenth.github.io/emmeans/reference/emmGrid-class.md)
  :

  The `emmGrid` class

- [`emmobj()`](https://rvlenth.github.io/emmeans/reference/emmobj.md) :

  Construct an `emmGrid` object from scratch

- [`contrast(`*`<emm_list>`*`)`](https://rvlenth.github.io/emmeans/reference/emm_list-object.md)
  [`pairs(`*`<emm_list>`*`)`](https://rvlenth.github.io/emmeans/reference/emm_list-object.md)
  [`test(`*`<emm_list>`*`)`](https://rvlenth.github.io/emmeans/reference/emm_list-object.md)
  [`confint(`*`<emm_list>`*`)`](https://rvlenth.github.io/emmeans/reference/emm_list-object.md)
  [`plot(`*`<emm_list>`*`)`](https://rvlenth.github.io/emmeans/reference/emm_list-object.md)
  [`coef(`*`<emm_list>`*`)`](https://rvlenth.github.io/emmeans/reference/emm_list-object.md)
  [`linfct(`*`<emm_list>`*`)`](https://rvlenth.github.io/emmeans/reference/emm_list-object.md)
  [`str(`*`<emm_list>`*`)`](https://rvlenth.github.io/emmeans/reference/emm_list-object.md)
  [`summary(`*`<emm_list>`*`)`](https://rvlenth.github.io/emmeans/reference/emm_list-object.md)
  [`print(`*`<emm_list>`*`)`](https://rvlenth.github.io/emmeans/reference/emm_list-object.md)
  [`as.data.frame(`*`<emm_list>`*`)`](https://rvlenth.github.io/emmeans/reference/emm_list-object.md)
  [`as.data.frame(`*`<summary_eml>`*`)`](https://rvlenth.github.io/emmeans/reference/emm_list-object.md)
  :

  The `emm_list` class

## Utilities

- [`emm()`](https://rvlenth.github.io/emmeans/reference/glht-support.md)
  [`as.glht()`](https://rvlenth.github.io/emmeans/reference/glht-support.md)
  :

  Support for
  [`multcomp::glht`](https://rdrr.io/pkg/multcomp/man/glht.html)

- [`emm_example()`](https://rvlenth.github.io/emmeans/reference/emm_example.md)
  : Run or list additional examples

- [`emm_options()`](https://rvlenth.github.io/emmeans/reference/emm_options.md)
  [`get_emm_option()`](https://rvlenth.github.io/emmeans/reference/emm_options.md)
  [`with_emm_options()`](https://rvlenth.github.io/emmeans/reference/emm_options.md)
  [`emm_defaults`](https://rvlenth.github.io/emmeans/reference/emm_options.md)
  : Set or change emmeans options

- [`str(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/emmGrid-methods.md)
  [`print(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/emmGrid-methods.md)
  [`vcov(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/emmGrid-methods.md)
  [`linfct()`](https://rvlenth.github.io/emmeans/reference/emmGrid-methods.md)
  :

  Miscellaneous methods for `emmGrid` objects

- [`make.tran()`](https://rvlenth.github.io/emmeans/reference/make.tran.md)
  [`inverse()`](https://rvlenth.github.io/emmeans/reference/make.tran.md)
  : Response-transformation extensions

## Convert to outside formats

- [`as.list(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/as.emmGrid.md)
  [`as.emm_list()`](https://rvlenth.github.io/emmeans/reference/as.emmGrid.md)
  [`as.emmGrid()`](https://rvlenth.github.io/emmeans/reference/as.emmGrid.md)
  :

  Convert to and from `emmGrid` objects

- [`as.mcmc(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/mcmc-support.md)
  [`as.mcmc(`*`<emm_list>`*`)`](https://rvlenth.github.io/emmeans/reference/mcmc-support.md)
  [`as.mcmc.list(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/mcmc-support.md)
  [`as.mcmc.list(`*`<emm_list>`*`)`](https://rvlenth.github.io/emmeans/reference/mcmc-support.md)
  : Support for MCMC-based estimation

- [`xtable(`*`<emmGrid>`*`)`](https://rvlenth.github.io/emmeans/reference/xtable.emmGrid.md)
  [`xtable(`*`<summary_emm>`*`)`](https://rvlenth.github.io/emmeans/reference/xtable.emmGrid.md)
  [`print(`*`<xtable_emm>`*`)`](https://rvlenth.github.io/emmeans/reference/xtable.emmGrid.md)
  :

  Using `xtable` for EMMs

## Models and extending

- [`models`](https://rvlenth.github.io/emmeans/reference/models.md) :

  Models supported in emmeans

- [`recover_data()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
  [`emm_basis()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
  [`.emm_register()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
  [`.std.link.labels()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
  [`.combine.terms()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
  [`.aovlist.dffun()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
  [`.cmpMM()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
  [`.get.excl()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
  [`.get.offset()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
  [`.my.vcov()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
  [`.all.vars()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
  [`.diag()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
  [`.num.key()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
  [`.emm_vignette()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
  [`.hurdle.support()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
  [`.zi.support()`](https://rvlenth.github.io/emmeans/reference/extending-emmeans.md)
  : Support functions for model extensions
