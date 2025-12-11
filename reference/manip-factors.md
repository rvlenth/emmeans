# Manipulate factors in a reference grid

These functions manipulate the levels of factors comprising a reference
grid by combining factor levels, splitting a factor's levels into
combinations of newly-defined factors, creating a grouping factor in
which factor(s) levels are nested, or permuting the order of levels of a
factor

## Usage

``` r
comb_facs(object, facs, newname = paste(facs, collapse = "."),
  drop = FALSE, ...)

split_fac(object, fac, newfacs, ...)

add_grouping(object, newname, refname, newlevs, ...)

add_submodels(object, ..., newname = "model")

permute_levels(object, fac, pos)
```

## Arguments

- object:

  An object of class `emmGrid`

- facs:

  Character vector. The names of the factors to combine

- newname:

  Character name of grouping factor to add (different from any existing
  factor in the grid)

- drop:

  Logical value. If `TRUE`, any levels of the new factor that are
  dropped if all occurrences in the newly reconstructed object have
  weight zero. If `FALSE`, all levels are retained. (This argument is
  ignored if there is no `.wgt.` column in `object@grid`.)

- ...:

  arguments passed to other methods

- fac:

  The name of a factor that is part of the grid in `object`

- newfacs:

  A named list with the names of new factors and their levels. The names
  must not already exist in the object, and the product of the lengths
  of the levels must equal the number of levels of `fac`.

- refname:

  Character name(s) of the reference factor(s)

- newlevs:

  Character vector or factor of the same length as that of the
  (combined) levels for `refname`. The grouping factor `newname` will
  have the unique values of `newlevs` as its levels. The order of levels
  in `newlevs` is the same as the order of the level combinations
  produced by [`expand.grid`](https://rdrr.io/r/base/expand.grid.html)
  applied to the levels of `refname` – that is, the first factor's
  levels change the fastest and the last one's vary the slowest.

- pos:

  Integer vector consisting of some permutation of the sequence `1:k`,
  where `k` is the number of levels of `fac`. This determines which
  position each level of `fac` will occupy after the levels are
  permuted; thus, if the levels of `fac` are `A,B,C,D`, and
  `pos = c(3,1,2,4)`, then the permuted levels will be `B,C,A,D`.

## Value

A modified object of class `emmGrid`

## The `comb_facs` function

`comb_facs` combines the levels of factors into a single factor in the
reference grid (similar to
[`interaction`](https://rdrr.io/r/base/interaction.html)). This new
factor replaces the factors that comprise it.

*Additional note:* The choice of whether to drop levels or not can make
a profound difference. If the goal is to combine factors for use in
`joint_tests`, we advise *against* `drop = TRUE` because that might
change the weights used in deriving marginal means. If combining factors
in a nested structure, dropping unused cases can considerably reduce the
storage required.

## The `split_fac` function

The levels in `newfacs` are expanded via
[`expand.grid`](https://rdrr.io/r/base/expand.grid.html) into
combinations of levels, and the factor `fac` is replaced by those factor
combinations. Unlike `add_grouping`, this creates a crossed, rather than
a nested structure. Note that the order of factor combinations is
systematic with the levels of first factor in `newfacs` varying the
fastest; and those factor combinations are assigned respectively to the
levels of `fac` as displayed in `str(object)`.

## The `add_grouping` function

This function adds a grouping factor to an existing reference grid or
other `emmGrid` object, such that the levels of one or more existing
factors (call them the reference factors) are mapped to a smaller number
of levels of the new grouping factor. The reference factors are then
nested in a new grouping factor named `newname`, and a new nesting
structure `refname %in% newname`. This facilitates obtaining marginal
means of the grouping factor, and contrasts thereof.

*Additional notes:* By default, the levels of `newname` will be ordered
alphabetically. To dictate a different ordering of levels, supply
`newlevs` as a `factor` having its levels in the desired order.

When `refname` specifies more than one factor, this can fundamentally
(and permanently) change what is meant by the levels of those individual
factors. For instance, in the `gwrg` example below, there are two levels
of `wool` nested in each `prod`; and that implies that we now regard
these as four different kinds of wool. Similarly, there are five
different tensions (L, M, H in prod 1, and L, M in prod 2).

## The `add_submodels` function

This function updates `object` with a named list of submodels specified
in `...`. These are `rbind`ed together and the corresponding rows for
each submodel are assigned a factor named `newname` with levels equal to
the names in `...`. This facilitates comparing estimates obtained from
different submodels. For this to work, the underlying model object must
be of a class supported by the `submodel` argument of
[`update.emmGrid`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md).

## The `permute_levels` function

This function permutes the levels of `fac`. The returned object has the
same factors, same `by` variables, but with the levels of `fac`
permuted. The order of the columns in `object@grid` may be altered.

NOTE: `fac` must not be nested in another factor. `permute_levels`
throws an error when `fac` is nested.

NOTE: Permuting the levels of a numeric predictor is tricky. For
example, if you want to display the new ordering of levels in
[`emmip()`](https://rvlenth.github.io/emmeans/reference/emmip.md), you
must add the arguments `style = "factor"` and `nesting.order = TRUE`.

## Examples

``` r
mtcars.lm <- lm(mpg ~ factor(vs)+factor(cyl)*factor(gear), data = mtcars)
(v.c.g <- ref_grid(mtcars.lm))
#>  vs cyl gear prediction    SE df
#>   0   4    3       21.7 4.420 23
#>   1   4    3       21.5 3.420 23
#>   0   6    3       19.9 3.690 23
#>   1   6    3       19.8 2.420 23
#>   0   8    3       15.1 0.987 23
#>   1   8    3       14.8 2.960 23
#>   0   4    4       27.1 3.040 23
#>   1   4    4       26.9 1.210 23
#>   0   6    4       19.9 2.210 23
#>   1   6    4       19.6 2.210 23
#>   0   8    4     nonEst    NA NA
#>   1   8    4     nonEst    NA NA
#>   0   4    5       28.3 2.790 23
#>   1   4    5       28.1 2.790 23
#>   0   6    5       19.7 3.420 23
#>   1   6    5       19.5 4.420 23
#>   0   8    5       15.4 2.420 23
#>   1   8    5       15.2 3.690 23
#> 
(v.cg <- comb_facs(v.c.g, c("cyl", "gear")))
#>  vs cyl.gear prediction    SE df
#>   0 4:3            21.7 4.420 23
#>   1 4:3            21.5 3.420 23
#>   0 6:3            19.9 3.690 23
#>   1 6:3            19.8 2.420 23
#>   0 8:3            15.1 0.987 23
#>   1 8:3            14.8 2.960 23
#>   0 4:4            27.1 3.040 23
#>   1 4:4            26.9 1.210 23
#>   0 6:4            19.9 2.210 23
#>   1 6:4            19.6 2.210 23
#>   0 8:4          nonEst    NA NA
#>   1 8:4          nonEst    NA NA
#>   0 4:5            28.3 2.790 23
#>   1 4:5            28.1 2.790 23
#>   0 6:5            19.7 3.420 23
#>   1 6:5            19.5 4.420 23
#>   0 8:5            15.4 2.420 23
#>   1 8:5            15.2 3.690 23
#> 
  
# One use is obtaining a single test for the joint contributions of two factors:
joint_tests(v.c.g)
#>  model term   df1 df2 F.ratio p.value note
#>  vs             1  23   0.005  0.9435     
#>  cyl            1  23   6.570  0.0174    e
#>  gear           1  23   0.755  0.3939    e
#>  cyl:gear       3  23   0.681  0.5725    e
#>  (confounded)   2  23  26.980 <0.0001     
#> 
#> e: df1 reduced due to non-estimability 

joint_tests(v.cg)
#>  model term df1 df2 F.ratio p.value note
#>  vs           1  23   0.005  0.9435     
#>  cyl.gear     7  23   4.402  0.0031    e
#> 
#> e: df1 reduced due to non-estimability 

# undo the 'comb_facs' operation:
split_fac(v.cg, "cyl.gear", list(cyl = c(4, 6, 8), gear = 3:5))
#>  vs cyl gear cyl.gear prediction    SE df
#>   0   4    3 4:3            21.7 4.420 23
#>   1   4    3 4:3            21.5 3.420 23
#>   0   6    3 6:3            19.9 3.690 23
#>   1   6    3 6:3            19.8 2.420 23
#>   0   8    3 8:3            15.1 0.987 23
#>   1   8    3 8:3            14.8 2.960 23
#>   0   4    4 4:4            27.1 3.040 23
#>   1   4    4 4:4            26.9 1.210 23
#>   0   6    4 6:4            19.9 2.210 23
#>   1   6    4 6:4            19.6 2.210 23
#>   0   8    4 8:4          nonEst    NA NA
#>   1   8    4 8:4          nonEst    NA NA
#>   0   4    5 4:5            28.3 2.790 23
#>   1   4    5 4:5            28.1 2.790 23
#>   0   6    5 6:5            19.7 3.420 23
#>   1   6    5 6:5            19.5 4.420 23
#>   0   8    5 8:5            15.4 2.420 23
#>   1   8    5 8:5            15.2 3.690 23
#> 

IS.glm <- glm(count ~ spray, data = InsectSprays, family = poisson)
IS.emm <- emmeans(IS.glm, "spray")
IS.new <- split_fac(IS.emm, "spray", list(A = 1:2, B = c("low", "med", "hi")))
str(IS.new)
#> 'emmGrid' object with variables:
#>     A = 1, 2
#>     B = low, med, hi
#> Transformation: “log” 

fiber.lm <- lm(strength ~ diameter + machine, data = fiber)
( frg <- ref_grid(fiber.lm) )
#>  diameter machine prediction    SE df
#>      24.1 A             40.4 0.724 11
#>      24.1 B             41.4 0.744 11
#>      24.1 C             38.8 0.788 11
#> 

# Suppose the machines are two different brands
brands <- factor(c("FiberPro", "FiberPro", "Acme"), levels = c("FiberPro", "Acme"))
( gfrg <- add_grouping(frg, "brand", "machine", brands) )
#>  diameter machine brand    prediction    SE df
#>      24.1 A       FiberPro       40.4 0.724 11
#>      24.1 B       FiberPro       41.4 0.744 11
#>      24.1 C       Acme           38.8 0.788 11
#> 

emmeans(gfrg, "machine")
#>  machine brand    emmean    SE df lower.CL upper.CL
#>  A       FiberPro   40.4 0.724 11     38.8     42.0
#>  B       FiberPro   41.4 0.744 11     39.8     43.1
#>  C       Acme       38.8 0.788 11     37.1     40.5
#> 
#> Confidence level used: 0.95 

emmeans(gfrg, "brand")
#>  brand    emmean    SE df lower.CL upper.CL
#>  FiberPro   40.9 0.531 11     39.7     42.1
#>  Acme       38.8 0.788 11     37.1     40.5
#> 
#> Results are averaged over the levels of: machine 
#> Confidence level used: 0.95 

### More than one reference factor
warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
gwrg <- add_grouping(ref_grid(warp.lm), 
    "prod",  c("tension", "wool"),  c(2, 1, 1,  1, 2, 1))
        # level combinations:         LA MA HA  LB MB HB

emmeans(gwrg, ~ wool * tension)   # some NAs due to impossible combinations
#>  wool tension prod emmean   SE df lower.CL upper.CL
#>  A    L       1    nonEst   NA NA       NA       NA
#>  B    L       1      28.2 3.65 48     20.9     35.6
#>  A    M       1      24.0 3.65 48     16.7     31.3
#>  B    M       1    nonEst   NA NA       NA       NA
#>  A    H       1      24.6 3.65 48     17.2     31.9
#>  B    H       1      18.8 3.65 48     11.4     26.1
#>  A    L       2      44.6 3.65 48     37.2     51.9
#>  B    L       2    nonEst   NA NA       NA       NA
#>  A    M       2    nonEst   NA NA       NA       NA
#>  B    M       2      28.8 3.65 48     21.4     36.1
#> 
#> Confidence level used: 0.95 

emmeans(gwrg, "prod")
#>  prod emmean   SE df lower.CL upper.CL
#>  1      24.6 1.92 48     20.8     28.5
#>  2      36.7 2.58 48     31.5     41.9
#> 
#> Results are averaged over the levels of: wool, tension 
#> Confidence level used: 0.95 

## Using 'add_submodels' to compare adjusted and unadjusted means
fibint.lm <- lm(strength ~ machine * diameter, data = fiber)
fibsub <- add_submodels(emmeans(fibint.lm, "machine"), 
    full = ~ ., additive = ~ . - machine:diameter, unadj = ~ machine)
#> NOTE: Results may be misleading due to involvement in interactions
emmeans(fibsub, pairwise ~ model | machine, adjust = "none")
#> $emmeans
#> machine = A:
#>  model    emmean    SE df lower.CL upper.CL
#>  full       40.2 0.777  9     38.5     42.0
#>  additive   40.4 0.760  9     38.7     42.1
#>  unadj      41.4 0.749  9     39.7     43.1
#> 
#> machine = B:
#>  model    emmean    SE df lower.CL upper.CL
#>  full       41.6 0.858  9     39.7     43.5
#>  additive   41.4 0.782  9     39.7     43.2
#>  unadj      43.2 0.749  9     41.5     44.9
#> 
#> machine = C:
#>  model    emmean    SE df lower.CL upper.CL
#>  full       38.5 0.966  9     36.3     40.7
#>  additive   38.8 0.827  9     36.9     40.7
#>  unadj      36.0 0.749  9     34.3     37.7
#> 
#> Confidence level used: 0.95 
#> 
#> $contrasts
#> machine = A:
#>  contrast         estimate    SE df t.ratio p.value
#>  full - additive    -0.160 0.162  9  -0.987  0.3492
#>  full - unadj       -1.178 0.207  9  -5.702  0.0003
#>  additive - unadj   -1.018 0.128  9  -7.966 <0.0001
#> 
#> machine = B:
#>  contrast         estimate    SE df t.ratio p.value
#>  full - additive     0.181 0.353  9   0.512  0.6209
#>  full - unadj       -1.600 0.418  9  -3.830  0.0040
#>  additive - unadj   -1.781 0.224  9  -7.966 <0.0001
#> 
#> machine = C:
#>  contrast         estimate    SE df t.ratio p.value
#>  full - additive    -0.263 0.499  9  -0.528  0.6105
#>  full - unadj        2.535 0.610  9   4.153  0.0025
#>  additive - unadj    2.798 0.351  9   7.966 <0.0001
#> 
#> 

# Permuting factor levels...
str(v.c.g)
#> 'emmGrid' object with variables:
#>     vs = 0, 1
#>     cyl = 4, 6, 8
#>     gear = 3, 4, 5
#> Some things are non-estimable (null space dim = 1)
str(permute_levels(v.c.g, "cyl", c(2,3,1)))
#> 'emmGrid' object with variables:
#>     vs = 0, 1
#>     cyl = 8, 4, 6
#>     gear = 3, 4, 5
#> Some things are non-estimable (null space dim = 1)
```
