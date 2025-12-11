# Contrasts and linear functions of EMMs

These methods provide for follow-up analyses of `emmGrid` objects:
Contrasts, pairwise comparisons, tests, and confidence intervals. They
may also be used to compute arbitrary linear functions of predictions or
EMMs.

## Usage

``` r
contrast(object, ...)

# S3 method for class 'emmGrid'
contrast(object, method = "eff", interaction = FALSE, by,
  offset = NULL, scale = NULL, name = "contrast",
  options = get_emm_option("contrast"), type, adjust, simple,
  combine = FALSE, ratios = TRUE, parens, enhance.levels = TRUE, wts,
  ...)

# S3 method for class 'emmGrid'
pairs(x, reverse = FALSE, ...)

# S3 method for class 'emmGrid'
coef(object, ...)

# S3 method for class 'emmGrid'
weights(object, ...)
```

## Arguments

- object:

  An object of class `emmGrid`

- ...:

  Additional arguments passed to other methods

- method:

  Character value giving the root name of a contrast method (e.g.
  `"pairwise"` – see
  [emmc-functions](https://rvlenth.github.io/emmeans/reference/emmc-functions.md)).
  Alternatively, a function of the same form, or a named `list` of
  coefficients (for a contrast or linear function) that must each
  conform to the number of results in each `by` group. In a multi-factor
  situation, the factor levels are combined and treated like a single
  factor.

- interaction:

  Character vector, logical value, or list. If this is specified,
  `method` is ignored. See the “Interaction contrasts” section below for
  details.

- by:

  Character names of variable(s) to be used for “by” groups. The
  contrasts or joint tests will be evaluated separately for each
  combination of these variables. If `object` was created with by
  groups, those are used unless overridden. Use `by = NULL` to use no by
  groups at all.

- offset, scale:

  Numeric vectors of the same length as each `by` group. The `scale`
  values, if supplied, multiply their respective linear estimates, and
  any `offset` values are added. Scalar values are also allowed. (These
  arguments are ignored when `interaction` is specified.)

- name:

  Character name to use to override the default label for contrasts used
  in table headings or subsequent contrasts of the returned object.

- options:

  If non-`NULL`, a named `list` of arguments to pass to
  [`update.emmGrid`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md),
  just after the object is constructed.

- type:

  Character: prediction type (e.g., `"response"`) – added to `options`

- adjust:

  Character: adjustment method (e.g., `"bonferroni"`) – added to
  `options`

- simple:

  Character vector or list: Specify the factor(s) *not* in `by`, or a
  list thereof. See the section below on simple contrasts.

- combine:

  Logical value that determines what is returned when `simple` is a
  list. See the section on simple contrasts.

- ratios:

  Logical value determining how log and logit transforms are handled.
  These transformations are exceptional cases in that there is a valid
  way to back-transform contrasts: differences of logs are logs of
  ratios, and differences of logits are odds ratios. If `ratios = TRUE`
  and summarized with `type = "response"`, `contrast` results are
  back-transformed to ratios whenever we have true contrasts
  (coefficients sum to zero). For other transformations, there is no
  natural way to back-transform contrasts, so even when summarized with
  `type = "response"`, contrasts are computed and displayed on the
  linear-predictor scale. Similarly, if `ratios = FALSE`, log and logit
  transforms are treated in the same way as any other transformation.

- parens:

  character or `NULL`. If a character value, the labels for levels being
  contrasted are parenthesized if they match the regular expression in
  `parens[1]` (via [`grep`](https://rdrr.io/r/base/grep.html)). The
  default is `emm_option("parens")`. Optionally, `parens` may contain
  second and third elements specifying what to use for left and right
  parentheses (default `"("` and `")"`). Specify `parens = NULL` or
  `parens = "a^"` (which won't match anything) to disable all
  parenthesization.

- enhance.levels:

  character or logical. If character, the levels of the named factors
  that are contrasted are enhanced by appending the name of the factor;
  e.g., if a factor named `"trt"` has levels `A` and `B`, a `trt`
  comparison is labeled `trtA - trtB`. If `enhance.levels` is logical,
  then if `TRUE` (the default), only factors with numeric levels are
  enhanced; and of course if `FALSE`, no levels are enhanced. The levels
  of `by` variables are not enhanced, and any names of factors that
  don't exist are silently ignored. To enhance the labels beyond what is
  done here, change them directly via
  [`levels<-`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md).

- wts:

  The `wts` argument for some contrast methods. You should omit this
  argument unless you want unequal `wts`. Otherwise we recommend
  specifying `wts = NA` which instructs that `wts` be obtained from
  `object`, *separately* for each `by` group. If numerical `wts` are
  specified, they must conform to the number of levels in each `by`
  group, and those same weights are used in each group.

- x:

  An `emmGrid` object

- reverse:

  Logical value - determines whether to use `"pairwise"` (if `TRUE`) or
  `"revpairwise"` (if `FALSE`).

## Value

`contrast` and `pairs` return an object of class `emmGrid`. Its grid
will correspond to the levels of the contrasts and any `by` variables.
The exception is that an
[`emm_list`](https://rvlenth.github.io/emmeans/reference/emm_list-object.md)
object is returned if `simple` is a list and `combine` is `FALSE`.

`coef` returns a `data.frame` containing the "parent" object's grid,
along with columns named `c.1, c.2, ...` containing the contrast
coefficients used to produce the linear functions embodied in the
object. [`coef()`](https://rdrr.io/r/stats/coef.html) only returns
coefficients if `object` is the result of a call to `contrast()`, and
the parent object is the object that was handed to `contrast`. This is
most useful for understanding interaction contrasts.

`weights` returns the weights stored for each row of `object`, or a
vector of 1s if no weights are saved.

## Note

When `object` has a nesting structure (this can be seen via
`str(object)`), then any grouping factors involved are forced into
service as `by` variables, and the contrasts are thus computed
separately in each nest. This in turn may lead to an irregular grid in
the returned `emmGrid` object, which may not be valid for subsequent
`emmeans` calls.

## Pairs method

The call `pairs(object)` is equivalent to
`contrast(object, method = "pairwise")`; and
`pairs(object, reverse = TRUE)` is the same as
`contrast(object, method = "revpairwise")`.

## Interaction contrasts

When `interaction` is specified, interaction contrasts are computed.
Specifically contrasts are generated for each factor separately, one at
a time; and these contrasts are applied to the object (the first time
around) or to the previous result (subsequently). (Any factors specified
in `by` are skipped.) The final result comprises contrasts of contrasts,
or, equivalently, products of contrasts for the factors involved. Any
named elements of `interaction` are assigned to contrast methods; others
are assigned in order of appearance in `object@levels`. The contrast
factors in the resulting `emmGrid` object are ordered the same as in
`interaction`.

`interaction` may be a character vector or list of valid contrast
methods (as documented for the `method` argument). If the vector or list
is shorter than the number needed, it is recycled. Alternatively, if the
user specifies `contrast = TRUE`, the contrast specified in `method` is
used for all factors involved.

## Simple contrasts

`simple` is essentially the complement of `by`: When `simple` is a
character vector, `by` is set to all the factors in the grid *except*
those in `simple`. If `simple` is a list, each element is used in turn
as `simple`, and assembled in an `"emm_list"`. To generate *all* simple
main effects, use `simple = "each"` (this works unless there actually is
a factor named `"each"`). Note that a non-missing `simple` will cause
`by` to be ignored.

Ordinarily, when `simple` is a list or `"each"`, the return value is an
[`emm_list`](https://rvlenth.github.io/emmeans/reference/emm_list-object.md)
object with each entry in correspondence with the entries of `simple`.
However, with `combine = TRUE`, the elements are all combined into one
family of contrasts in a single
[`emmGrid`](https://rvlenth.github.io/emmeans/reference/emmGrid-class.md)
object using
[`rbind.emmGrid`](https://rvlenth.github.io/emmeans/reference/rbind.emmGrid.md)..
In that case, the `adjust` argument sets the adjustment method for the
combined set of contrasts.

## Examples

``` r
warp.lm <- lm(breaks ~ wool*tension, data = warpbreaks)
(warp.emm <- emmeans(warp.lm, ~ tension | wool))
#> wool = A:
#>  tension emmean   SE df lower.CL upper.CL
#>  L         44.6 3.65 48     37.2     51.9
#>  M         24.0 3.65 48     16.7     31.3
#>  H         24.6 3.65 48     17.2     31.9
#> 
#> wool = B:
#>  tension emmean   SE df lower.CL upper.CL
#>  L         28.2 3.65 48     20.9     35.6
#>  M         28.8 3.65 48     21.4     36.1
#>  H         18.8 3.65 48     11.4     26.1
#> 
#> Confidence level used: 0.95 

contrast(warp.emm, "poly")    # inherits 'by = "wool"' from warp.emm
#> wool = A:
#>  contrast  estimate   SE df t.ratio p.value
#>  linear      -20.00 5.16 48  -3.878  0.0003
#>  quadratic    21.11 8.93 48   2.363  0.0222
#> 
#> wool = B:
#>  contrast  estimate   SE df t.ratio p.value
#>  linear       -9.44 5.16 48  -1.831  0.0733
#>  quadratic   -10.56 8.93 48  -1.182  0.2432
#> 

### Custom contrast coefs (we already have wool as 'by' thus 3 means to contrast)
contrast(warp.emm, list(mid.vs.ends = c(-1,2,-1)/2, lo.vs.hi = c(1,0,-1)))
#> wool = A:
#>  contrast    estimate   SE df t.ratio p.value
#>  mid.vs.ends   -10.56 4.47 48  -2.363  0.0222
#>  lo.vs.hi       20.00 5.16 48   3.878  0.0003
#> 
#> wool = B:
#>  contrast    estimate   SE df t.ratio p.value
#>  mid.vs.ends     5.28 4.47 48   1.182  0.2432
#>  lo.vs.hi        9.44 5.16 48   1.831  0.0733
#> 

pairs(warp.emm)
#> wool = A:
#>  contrast estimate   SE df t.ratio p.value
#>  L - M      20.556 5.16 48   3.986  0.0007
#>  L - H      20.000 5.16 48   3.878  0.0009
#>  M - H      -0.556 5.16 48  -0.108  0.9936
#> 
#> wool = B:
#>  contrast estimate   SE df t.ratio p.value
#>  L - M      -0.556 5.16 48  -0.108  0.9936
#>  L - H       9.444 5.16 48   1.831  0.1704
#>  M - H      10.000 5.16 48   1.939  0.1389
#> 
#> P value adjustment: tukey method for comparing a family of 3 estimates 

# Effects (dev from mean) of the 6 factor combs, with enhanced levels:
contrast(warp.emm, "eff", by = NULL, 
    enhance.levels = c("wool", "tension"))  
#>  contrast              estimate   SE df t.ratio p.value
#>  tensionL woolA effect  16.4074 3.33 48   4.929 <0.0001
#>  tensionM woolA effect  -4.1481 3.33 48  -1.246  0.4289
#>  tensionH woolA effect  -3.5926 3.33 48  -1.079  0.4289
#>  tensionL woolB effect   0.0741 3.33 48   0.022  0.9823
#>  tensionM woolB effect   0.6296 3.33 48   0.189  0.9823
#>  tensionH woolB effect  -9.3704 3.33 48  -2.815  0.0212
#> 
#> P value adjustment: fdr method for 6 tests 
    
pairs(warp.emm, simple = "wool") # same as pairs(warp.emm, by = "tension")
#> tension = L:
#>  contrast estimate   SE df t.ratio p.value
#>  A - B       16.33 5.16 48   3.167  0.0027
#> 
#> tension = M:
#>  contrast estimate   SE df t.ratio p.value
#>  A - B       -4.78 5.16 48  -0.926  0.3589
#> 
#> tension = H:
#>  contrast estimate   SE df t.ratio p.value
#>  A - B        5.78 5.16 48   1.120  0.2682
#> 

# Do all "simple" comparisons, combined into one family
pairs(warp.emm, simple = "each", combine = TRUE)
#>  wool tension contrast estimate   SE df t.ratio p.value
#>  A    .       L - M      20.556 5.16 48   3.986  0.0021
#>  A    .       L - H      20.000 5.16 48   3.878  0.0029
#>  A    .       M - H      -0.556 5.16 48  -0.108  1.0000
#>  B    .       L - M      -0.556 5.16 48  -0.108  1.0000
#>  B    .       L - H       9.444 5.16 48   1.831  0.6594
#>  B    .       M - H      10.000 5.16 48   1.939  0.5255
#>  .    L       A - B      16.333 5.16 48   3.167  0.0241
#>  .    M       A - B      -4.778 5.16 48  -0.926  1.0000
#>  .    H       A - B       5.778 5.16 48   1.120  1.0000
#> 
#> P value adjustment: bonferroni method for 9 tests 

if (FALSE) { # \dontrun{

## Note that the following are NOT the same:
contrast(warp.emm, simple = c("wool", "tension"))
contrast(warp.emm, simple = list("wool", "tension"))
## The first generates contrasts for combinations of wool and tension
##   (same as by = NULL)
## The second generates contrasts for wool by tension, and for 
##   tension by wool, respectively.
} # }

# An interaction contrast for tension:wool
tw.emm <- contrast(warp.emm, interaction = c(tension = "poly", wool = "consec"), 
                   by = NULL)
tw.emm          # see the estimates
#>  tension_poly wool_consec estimate    SE df t.ratio p.value
#>  linear       B - A           10.6  7.29 48   1.447  0.1543
#>  quadratic    B - A          -31.7 12.60 48  -2.507  0.0156
#> 
coef(tw.emm)    # see the contrast coefficients
#>   tension wool c.1 c.2
#> 1       L    A   1  -1
#> 2       M    A   0   2
#> 3       H    A  -1  -1
#> 4       L    B  -1   1
#> 5       M    B   0  -2
#> 6       H    B   1   1

# Use of scale and offset
#   an unusual use of the famous stack-loss data...
mod <- lm(Water.Temp ~ poly(stack.loss, degree = 2), data = stackloss)
(emm <- emmeans(mod, "stack.loss", at = list(stack.loss = 10 * (1:4))))
#>  stack.loss emmean    SE df lower.CL upper.CL
#>          10   18.8 0.463 18     17.9     19.8
#>          20   22.3 0.564 18     21.1     23.5
#>          30   24.9 0.646 18     23.5     26.3
#>          40   26.7 0.958 18     24.6     28.7
#> 
#> Confidence level used: 0.95 
# Convert results from Celsius to Fahrenheit:
confint(contrast(emm, "identity", scale = 9/5, offset = 32))
#>  contrast     estimate    SE df lower.CL upper.CL
#>  stack.loss10     65.9 0.833 18     64.1     67.6
#>  stack.loss20     72.1 1.020 18     70.0     74.3
#>  stack.loss30     76.8 1.160 18     74.4     79.3
#>  stack.loss40     80.0 1.720 18     76.4     83.6
#> 
#> Confidence level used: 0.95 
```
