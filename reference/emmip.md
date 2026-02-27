# Interaction-style plots for estimated marginal means

Creates an interaction plot of EMMs based on a fitted model and a simple
formula specification.

## Usage

``` r
emmip(object, formula, ...)

# Default S3 method
emmip(object, formula, type, CIs = FALSE, PIs = FALSE,
  style, engine = get_emm_option("graphics.engine"), plotit = TRUE,
  nesting.order = FALSE, abbr.len = c(0, 0), ...)

emmip_ggplot(emms, style = "factor", dodge = 0.2, xlab = labs$xlab,
  ylab = labs$ylab, tlab = labs$tlab, facetlab = "label_both", scale,
  dotarg = list(shape = "circle"), linearg = list(linetype = "solid"),
  CIarg = list(alpha = 0.4), PIarg = list(alpha = 0.25), col = NULL, ...)

emmip_lattice(emms, style = "factor", xlab = labs$xlab, ylab = labs$ylab,
  tlab = labs$tlab, pch = c(1, 2, 6, 7, 9, 10, 15:20), lty = 1,
  col = NULL, ...)
```

## Arguments

- object:

  An object of class `emmGrid`, or a fitted model of a class supported
  by the emmeans package

- formula:

  Formula of the form `trace.factors ~ x.factors | by.factors`. The EMMs
  are plotted against `x.factor` for each level of `trace.factors`.
  `by.factors` is optional, but if present, it determines separate
  panels. Each element of this formula may be a single factor in the
  model, or a combination of factors using the `*` operator.

- ...:

  Additional arguments passed to
  [`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md)
  (when `object` is not already an `emmGrid` object), `predict.emmGrid`,
  `emmip_ggplot`, or `emmip_lattice`.

- type:

  As in
  [`predict.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md),
  this determines whether we want to inverse-transform the predictions
  (`type = "response"`) or not (any other choice). The default is
  `"link"`, unless the `"predict.type"` option is in force; see
  [`emm_options`](https://rvlenth.github.io/emmeans/reference/emm_options.md).
  In addition, the user may specify `type = "scale"` to create a
  transformed scale for the vertical axis based on `object`'s response
  transformation or link function.

- CIs:

  Logical value. If `TRUE`, confidence intervals (or HPD intervals for
  Bayesian models) are added to the plot (works only with
  `engine = "ggplot"`).

- PIs:

  Logical value. If `TRUE`, prediction intervals are added to the plot
  (works only with `engine = "ggplot"`). This is allowed only if the
  underlying model family is `"gaussian"`. If both `CIs` and `PIs` are
  `TRUE`, the prediction intervals will be somewhat longer, lighter, and
  thinner than the confidence intervals. Additional parameters to
  [`predict.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
  (e.g., `sigma`) may be passed via `...`. For Bayesian models, PIs
  require `frequentist = TRUE` and a value for `sigma`.

- style:

  Optional character value. This has an effect only when the horizontal
  variable is a single numeric variable. If `style` is unspecified or
  `"numeric"`, the horizontal scale will be numeric and curves are
  plotted using lines (and no symbols). With `style = "factor"`, the
  horizontal variable is treated as the levels of a factor (equally
  spaced along the horizontal scale), and curves are plotted using lines
  and symbols. When the horizontal variable is character or factor, or a
  combination of more than one predictor, `"factor"` style is always
  used.

- engine:

  Character value matching `"ggplot"` (default), `"lattice"`, or
  `"none"`. The graphics engine to be used to produce the plot. These
  require, respectively, the ggplot2 or lattice package to be installed.
  Specifying `"none"` is equivalent to setting `plotit = FALSE`.

- plotit:

  Logical value. If `TRUE`, a graphical object is returned; if `FALSE`,
  a data.frame is returned containing all the values used to construct
  the plot.

- nesting.order:

  Logical value. If `TRUE`, factors that are nested are presented in
  order according to their nesting factors, even if those nesting
  factors are not present in `formula`. If `FALSE`, only the variables
  in `formula` are used to order the variables.

- abbr.len:

  Integer vector of length 1 or 2 used to specify that factor levels be
  abbreviated to the specified lengths, using
  [`abbreviate`](https://rdrr.io/r/base/abbreviate.html). Any values
  less than 1 are ignored. `abbr.len[1]` applies to a factor plotted
  along the horizontal axis, and `abbr.len[2]` (if present) applies to
  trace factors (shown in a legend) or by factors (shown in facet
  labels).

- emms:

  A `data.frame` created by calling `emmip` with `plotit = FALSE`.
  Certain variables and attributes are expected to exist in this data
  frame; see the section detailing the rendering functions.

- dodge:

  Numerical amount passed to
  [`ggplot2::position_dodge`](https://ggplot2.tidyverse.org/reference/position_dodge.html)
  by which points and intervals are offset so they do not collide.

- xlab, ylab, tlab:

  Character labels for the horizontal axis, vertical axis, and traces
  (the different curves), respectively. The `emmip` function generates
  these automatically and provides therm via the `labs` attribute, but
  the user may override these if desired.

- facetlab:

  Labeller for facets (when by variables are in play). Use
  `"label_value"` to show just the factor levels, or `"label_both"` to
  show both the factor names and factor levels; `"label_context"`
  decides which based on how many `by` factors there are. See the
  documentation for
  [`ggplot2::labellers`](https://ggplot2.tidyverse.org/reference/labellers.html).

- scale:

  If not missing, an object of class
  [`scales::trans`](https://scales.r-lib.org/reference/new_transform.html)
  specifying a (usually) nonlinear scaling for the vertical axis. For
  example, `scales = scales::log_trans()` specifies a logarithmic scale.
  For fine-tuning purposes, additional arguments to
  [`ggplot2::scale_y_continuous`](https://ggplot2.tidyverse.org/reference/scale_continuous.html)
  may be included in `...` .

- dotarg:

  `list` of arguments passed to `geom_point` to customize appearance of
  points. If this includes `shape` and it is of length greater than 1,
  it is interpreted as a shape palette instead.

- linearg:

  `list` of arguments passed to `geom_line` to customize appearance of
  lines. If this includes `linetype` and it is of length greater than 1,
  it is interpreted as a linetype.

- CIarg, PIarg:

  `list`s of arguments passed to `geom_linerange` to customize
  appearance of intervals. (Note: the `linetype` aesthetic defaults to
  `"solid"` under the hood)

- col:

  With `emmip_ggplot`, this adds `color = col` (not `colour`) to all of
  the `*arg` lists. This is intended for setting a common color for
  everything, such as a black-and-white plot. With `emmip_lattice`,
  `col` specifies the colors to use for each group, recycled as needed.
  If not specified, the default trellis colors are used.

- pch:

  (Lattice only) The plotting characters to use for each group (i.e.,
  levels of `trace.factors`). They are recycled as needed.

- lty:

  (Lattice only) The line types to use for each group. Recycled as
  needed.

## Value

If `plotit = FALSE`, a `data.frame` (actually, a `summary_emm` object)
with the table of EMMs that would be plotted. The variables plotted are
named `xvar` and `yvar`, and the trace factor is named `tvar`. This data
frame has an added `"labs"` attribute containing the labels `xlab`,
`ylab`, and `tlab` for these respective variables. The confidence limits
are also included, renamed `LCL` and `UCL`.

If `plotit = TRUE`, the function returns an object of class `"ggplot"`
or a `"trellis"`, depending on `engine`.

## Note

Conceptually, this function is equivalent to
[`interaction.plot`](https://rdrr.io/r/stats/interaction.plot.html)
where the summarization function is thought to return the EMMs.

## Details

If `object` is a fitted model,
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md) is
called with an appropriate specification to obtain estimated marginal
means for each combination of the factors present in `formula` (in
addition, any arguments in `...` that match `at`, `trend`, `cov.reduce`,
or `fac.reduce` are passed to `emmeans`). Otherwise, if `object` is an
`emmGrid` object, its first element is used, and it must contain one
estimate for each combination of the factors present in `formula`.

## Rendering functions

The functions `emmip_ggplot` and `emmip_lattice` are called when
`plotit == TRUE` to render the plots; but they may also be called later
on an object saved via `plotit = FALSE` (or `engine = "none"`). The
functions require that `emms` contains variables `xvar`, `yvar`, and
`tvar`, and attributes `"labs"` and `"vars"`. Confidence intervals are
plotted if variables `LCL` and `UCL` exist; and prediction intervals are
plotted if `LPL` and `UPL` exist. Finally, it must contain the variables
named in `attr(emms, "vars")`.

In `emmip_ggplot`, colors, linetypes, and shapes are all assigned to
groups (according to `tvar`) unless overridden. So, for example, one may
have different symbols for each group by simply specifying
`dotarg = list()`.

## See also

[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md),
[`interaction.plot`](https://rdrr.io/r/stats/interaction.plot.html).

## Examples

``` r
#--- Three-factor example
noise.lm = lm(noise ~ size * type * side, data = auto.noise)

# Separate interaction plots of size by type, for each side
emmip(noise.lm, type ~ size | side)


# Same using the "lattice" engine
emmip(noise.lm, type ~ size | side, engine = "lattice")


# One interaction plot, using combinations of size and side as the x factor
# ... with added confidence intervals and some formatting changes
emmip(noise.lm, type ~ side * size, CIs = TRUE,
    CIarg = list(linewidth = 1.5, alpha = 1, color = "orange"),
    dotarg = list(size = 2, shape = "square", color = "black"))

    
# Same using legacy theme
with_emm_options(gg.theme = 1,
    emmip(noise.lm, type ~ side * size, CIs = TRUE,
        CIarg = list(linewidth = 1.5, alpha = 1, color = "orange"),
        dotarg = list(size = 2, shape = "square", color = "black")))


# Create a black-and-white version of above with different linetypes
# (Let the linetypes and symbols default to the palette)
emmip(noise.lm, type ~ side * size, CIs = TRUE, col = "black",
      linearg = list(), dotarg = list(size = 4), CIarg = list(alpha = 1)) +
    ggplot2::theme_bw()


# One interaction plot using combinations of type and side as the trace factor
emmip(noise.lm, type * side ~ size)


# Individual traces in panels
emmip(noise.lm, ~ size | type * side)


# Example for the 'style' argument
fib.lm = lm(strength ~ machine * sqrt(diameter), data = fiber)
fib.rg = ref_grid(fib.lm, at = list(diameter = c(12, 14, 15, 25, 36)))
emmip(fib.rg, machine ~ diameter)   # curves (because diameter is numeric)

emmip(fib.rg, machine ~ diameter, style = "factor")  # points and lines


# For an example using extra ggplot2 code, see 'vignette("messy-data")',
# in the section on nested models.

### Options with transformations or link functions
neuralgia.glm <- glm(Pain ~ Treatment * Sex + Age, family = binomial(), 
                     data = neuralgia) 

# On link scale:
emmip(neuralgia.glm, Treatment ~ Sex)


# On response scale:
emmip(neuralgia.glm, Treatment ~ Sex, type = "response")


# With transformed axis scale and custom scale divisions
emmip(neuralgia.glm, Treatment ~ Sex, type = "scale",
    breaks = seq(0.10, 0.90, by = 0.10))
```
