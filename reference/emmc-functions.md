# Contrast families

Functions with an extension of `.emmc` provide for named contrast
families. One of the standard ones documented here may be used, or the
user may write such a function.

## Usage

``` r
pairwise.emmc(levs, exclude = integer(0), include, ...)

revpairwise.emmc(levs, exclude = integer(0), include, ...)

tukey.emmc(levs, reverse = FALSE, ...)

poly.emmc(levs, max.degree = min(6, k - 1), ...)

opoly.emmc(levs, max.degree = min(6, k - 1), scores, exclude = integer(0),
  include, ...)

trt.vs.ctrl.emmc(levs, ref = 1, reverse = FALSE, exclude = integer(0),
  include, ...)

trt.vs.ctrl1.emmc(levs, ref = 1, ...)

trt.vs.ctrlk.emmc(levs, ref = length(levs), ...)

dunnett.emmc(levs, ref = 1, ...)

eff.emmc(levs, exclude = integer(0), include, wts = rep(1, length(levs)),
  ...)

del.eff.emmc(levs, exclude = integer(0), include, wts = rep(1,
  length(levs)), ...)

consec.emmc(levs, reverse = FALSE, exclude = integer(0), include, ...)

mean_chg.emmc(levs, reverse = FALSE, exclude = integer(0), include, ...)

helmert.emmc(levs, exclude = integer(0), include, ...)

nrmlz.emmc(levs, family, ...)

wtcon.emmc(levs, wts = rep(1, length(levs)), cmtype = "GrandMean", ...)

identity.emmc(levs, exclude = integer(0), include, ...)
```

## Arguments

- levs:

  Vector of factor levels

- exclude:

  integer vector of indices, or character vector of levels to exclude
  from consideration. These levels will receive weight 0 in all
  contrasts. Character levels must exactly match elements of `levs`.

- include:

  integer or character vector of levels to include (the complement of
  `exclude`). An error will result if the user specifies both `exclude`
  and `include`.

- ...:

  Additional arguments, passed to related methods as appropriate

- reverse:

  Logical value to determine the direction of comparisons

- max.degree:

  Integer specifying the maximum degree of polynomial contrasts

- scores:

  Set of values of length `length(levs)` over which orthogonal
  polynomials are computed. The default scores are the consecutive
  integers `seq_along(levs)`. (If `exclude` or `include` are used, the
  default scores are subsetted accordingly; and if `scores` is
  specified, its length must be the same as that of the subsetted
  `levs`).

- ref:

  Integer(s) or character(s) specifying which level(s) to use as the
  reference. Character values must exactly match elements of `levs`
  (including any enhancements – see examples)

- wts:

  Optional weights to use with `eff.emmc` and `del.eff.emmc` contrasts.
  These default to equal weights. If `exclude` or `include` are
  specified, `wts` may be either the same length as `levs` or the length
  of the included levels. In the former case, weights for any excluded
  levels are set to zero. `wts` has no impact on the results unless
  there are at least three levels included in the contrast.

- family:

  name of contrast family to use

- cmtype:

  the `type` argument passed to
  [`contrMat`](https://rdrr.io/pkg/multcomp/man/contrMat.html)

## Value

A data.frame, each column containing contrast coefficients for levs. The
"desc" attribute is used to label the results in emmeans, and the
"adjust" attribute gives the default adjustment method for multiplicity.

## Details

Each standard contrast family has a default multiple-testing adjustment
as noted below. These adjustments are often only approximate; for a more
exacting adjustment, use the interfaces provided to `glht` in the
multcomp package.

`pairwise.emmc`, `revpairwise.emmc`, and `tukey.emmc` generate contrasts
for all pairwise comparisons among estimated marginal means at the
levels in levs. The distinction is in which direction they are
subtracted. For factor levels A, B, C, D, `pairwise.emmc` generates the
comparisons A-B, A-C, A-D, B-C, B-D, and C-D, whereas `revpairwise.emmc`
generates B-A, C-A, C-B, D-A, D-B, and D-C. `tukey.emmc` invokes
`pairwise.emmc` or `revpairwise.emmc` depending on `reverse`. The
default multiplicity adjustment method is `"tukey"`, which is only
approximate when the standard errors differ.

`poly.emmc` and `opoly.emmc` generate orthogonal polynomial contrasts.
`poly.emmc` uses equally-spaced factor levels; coefficients are derived
from the [`poly`](https://rdrr.io/r/stats/poly.html) function, but an
*ad hoc* algorithm is used to scale them to integer coefficients that
are (usually) the same as in published tables of orthogonal polynomial
contrasts. On the other hand, `opoly.emmc`'s coefficients are always
normalized (sum of squares equals 1), but allows the user to choose
alternate reference points in `scores`, as in the
[`contr.poly`](https://rdrr.io/r/stats/contrast.html) function. In both
cases, the default multiplicity adjustment method is `"none"`.

`trt.vs.ctrl.emmc` and its relatives generate contrasts for comparing
one level (or the average over specified levels) with each of the other
levels. The argument `ref` should be the index(es) (not the labels) of
the reference level(s). `trt.vs.ctrl1.emmc` is the same as
`trt.vs.ctrl.emmc` with a reference value of 1, and `trt.vs.ctrlk.emmc`
is the same as `trt.vs.ctrl` with a reference value of `length(levs)`.
`dunnett.emmc` is the same as `trt.vs.ctrl`. The default multiplicity
adjustment method is `"dunnettx"`, a close approximation to the Dunnett
adjustment. *Note* in all of these functions, it is illegal to have any
overlap between the `ref` levels and the `exclude` levels. If any is
found, an error is thrown.

`consec.emmc` and `mean_chg.emmc` are useful for contrasting treatments
that occur in sequence. For a factor with levels A, B, C, D,
`consec.emmc` generates the comparisons B-A, C-B, and D-C, while
`mean_chg.emmc` generates the contrasts (B+C+D)/3 - A, (C+D)/2 -
(A+B)/2, and D - (A+B+C)/3. With `reverse = TRUE`, these differences go
in the opposite direction.

`eff.emmc` and `del.eff.emmc` generate contrasts that compare each level
with the average over all levels (in `eff.emmc`) or over all other
levels (in `del.eff.emmc`). These differ only in how they are scaled.
For a set of k EMMs, `del.eff.emmc` gives weight 1 to one EMM and weight
-1/(k-1) to the others, while `eff.emmc` gives weights (k-1)/k and -1/k
respectively, as in subtracting the overall EMM from each EMM. The
default multiplicity adjustment method is `"fdr"`. This is a
Bonferroni-based method and is slightly conservative; see
[`p.adjust`](https://rdrr.io/r/stats/p.adjust.html).

`nrmlz.emmc` is a wrapper that can be used with any other `.emmc`
function that will normalize the contrast coefficients so that the sum
of its squares equals 1. Just provide the root name of the function in
`family`, along with any other arguments to pass to it.

`wtcon.emmc` generates weighted contrasts based on the function
[`contrMat`](https://rdrr.io/pkg/multcomp/man/contrMat.html) function in
the multcomp package, using the provided `type` as documented there. If
the user provides `wts`, they have to conform to the length of `levs`;
however, if `wts` is not specified, `contrast` will fill-in what is
required, and usually this is safer (especially when `by != NULL` which
usually means that the weights are different in each `by` group).

`identity.emmc` simply returns the identity matrix (as a data frame),
minus any columns specified in `exclude`. It is potentially useful in
cases where a contrast function must be specified, but none is desired.

## Note

Caution is needed in cases where the user alters the ordering of results
(e.g., using the the `"[...]"` operator), because the contrasts
generated depend on the order of the levels provided. For example,
suppose `trt.vs.ctrl1` contrasts are applied to two `by` groups with
levels ordered (Ctrl, T1, T2) and (T1, T2, Ctrl) respectively, then the
contrasts generated will be for (T1 - Ctrl, T2 - Ctrl) in the first
group and (T2 - T1, Ctrl - T1) in the second group, because the first
level in each group is used as the reference level.

## Examples

``` r
warp.lm <- lm(breaks ~ wool*tension, data = warpbreaks)
warp.emm <- emmeans(warp.lm, ~ tension | wool)
contrast(warp.emm, "poly")
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
contrast(warp.emm, "trt.vs.ctrl", ref = "M")
#> wool = A:
#>  contrast estimate   SE df t.ratio p.value
#>  L - M      20.556 5.16 48   3.986  0.0005
#>  H - M       0.556 5.16 48   0.108  0.9863
#> 
#> wool = B:
#>  contrast estimate   SE df t.ratio p.value
#>  L - M      -0.556 5.16 48  -0.108  0.9863
#>  H - M     -10.000 5.16 48  -1.939  0.1077
#> 
#> P value adjustment: dunnettx method for 2 tests 
if (FALSE) { # \dontrun{
## Same when enhanced labeling is used:
contrast(warp.emm, "trt.vs.ctrl", 
         enhance.levels = "tension", ref = "tensionM")} # }

# Comparisons with grand mean
contrast(warp.emm, "eff")
#> wool = A:
#>  contrast estimate   SE df t.ratio p.value
#>  L effect    13.52 2.98 48   4.540  0.0001
#>  M effect    -7.04 2.98 48  -2.363  0.0333
#>  H effect    -6.48 2.98 48  -2.177  0.0344
#> 
#> wool = B:
#>  contrast estimate   SE df t.ratio p.value
#>  L effect     2.96 2.98 48   0.995  0.3247
#>  M effect     3.52 2.98 48   1.182  0.3247
#>  H effect    -6.48 2.98 48  -2.177  0.1033
#> 
#> P value adjustment: fdr method for 3 tests 
# Comparisons with a weighted grand mean
contrast(warp.emm, "eff", wts = c(2, 5, 3))
#> wool = A:
#>  contrast estimate   SE df t.ratio p.value
#>  L effect    16.28 3.61 48   4.509  0.0001
#>  M effect    -4.28 2.25 48  -1.903  0.0946
#>  H effect    -3.72 3.22 48  -1.156  0.2535
#> 
#> wool = B:
#>  contrast estimate   SE df t.ratio p.value
#>  L effect     2.56 3.61 48   0.708  0.4824
#>  M effect     3.11 2.25 48   1.384  0.2592
#>  H effect    -6.89 3.22 48  -2.139  0.1127
#> 
#> P value adjustment: fdr method for 3 tests 

# Compare only low and high tensions
# Note pairs(emm, ...) calls contrast(emm, "pairwise", ...)
pairs(warp.emm, exclude = 2)
#> wool = A:
#>  contrast estimate   SE df t.ratio p.value
#>  L - H       20.00 5.16 48   3.878  0.0003
#> 
#> wool = B:
#>  contrast estimate   SE df t.ratio p.value
#>  L - H        9.44 5.16 48   1.831  0.0733
#> 
# (same results using exclude = "M" or include = c("L","H") or include = c(1,3))

# Same contrasts as above but with normalized contrast coefficients
contrast(warp.emm, "nrmlz", family = "pairwise", include = c(1, 3))
#> wool = A:
#>  contrast estimate   SE df t.ratio p.value
#>  L - H       14.14 3.65 48   3.878  0.0003
#> 
#> wool = B:
#>  contrast estimate   SE df t.ratio p.value
#>  L - H        6.68 3.65 48   1.831  0.0733
#> 

### Setting up a custom contrast function
revhelmert.emmc <- function(levs, ...) {
    M <- as.data.frame(contr.helmert(levs)[rev(seq_along(levs)), ])
    names(M) <- paste(rev(levs)[-1],"vs later")
    attr(M, "desc") <- "reverse Helmert contrasts"
    M
}
contrast(warp.emm, "revhelmert")
#> wool = A:
#>  contrast   estimate   SE df t.ratio p.value
#>  M vs later   -0.556 5.16 48  -0.108  0.9147
#>  L vs later   40.556 8.93 48   4.540 <0.0001
#> 
#> wool = B:
#>  contrast   estimate   SE df t.ratio p.value
#>  M vs later   10.000 5.16 48   1.939  0.0584
#>  L vs later    8.889 8.93 48   0.995  0.3247
#> 

if (FALSE) { # \dontrun{
# See what is used for polynomial contrasts with 6 levels
emmeans:::poly.emmc(1:6)
} # }
```
