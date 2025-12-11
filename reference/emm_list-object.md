# The `emm_list` class

An `emm_list` object is simply a list of
[`emmGrid`](https://rvlenth.github.io/emmeans/reference/emmGrid-class.md)
objects. Such a list is returned, for example, by
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md) with
a two-sided formula or a list as its `specs` argument. Several methods
for this class are provided, as detailed below. Typically, these methods
just quietly do the same thing as their `emmGrid` methods, using the
first element of the list. You can specify `which` to select a different
element, or just run the corresponding `emmGrid` method on
`object[[k]]`.

## Usage

``` r
# S3 method for class 'emm_list'
contrast(object, ..., which = NULL)

# S3 method for class 'emm_list'
pairs(x, ..., which = NULL)

# S3 method for class 'emm_list'
test(object, ..., which = seq_along(object))

# S3 method for class 'emm_list'
confint(object, ..., which = seq_along(object))

# S3 method for class 'emm_list'
plot(x, ..., which = 1)

# S3 method for class 'emm_list'
coef(object, ..., which = NULL)

# S3 method for class 'emm_list'
linfct(object, ..., which = seq_along(object))

# S3 method for class 'emm_list'
str(object, ...)

# S3 method for class 'emm_list'
summary(object, ..., which = seq_along(object))

# S3 method for class 'emm_list'
print(x, ...)

# S3 method for class 'emm_list'
as.data.frame(x, ...)

# S3 method for class 'summary_eml'
as.data.frame(x, row.names = NULL, optional = FALSE,
  which, ...)
```

## Arguments

- object, x:

  an object of class `emm_list`

- ...:

  additional arguments passed to corresponding `emmGrid` method. In
  addition, the user may include a logical argument `drop` that is akin
  to [`drop`](https://rdrr.io/r/base/drop.html) and the argument of the
  same name in subscripting matrices and data frames. When `drop` is
  `TRUE` (the default), then when the result is a `list` of length 1,
  the `list` structure is removed.

- which:

  integer vector specifying which elements to select; if `NULL`, we try
  to guess which elements make sense. Usually, this is all elements
  having names that start with ‘em’ or ‘ls’, or the first element if no
  matches are found. However, in `coef.emm_list`, these are the ones we
  *exclude*.

- row.names, optional:

  Required arguments of `as.data.frame`, ignored

## Value

a `list` of objects returned by the corresponding `emmGrid` method
(thus, often, another `emm_list` object). However, if `which` has length
1, the one result is not wrapped in a list.

`summary.emm_list` returns an object of class `summary_eml`, which is a
list of `summary_emm` objects.

The `as.data.frame` methods return a single data frame via
`as.data.frame(rbind(x))`. See also
[`rbind.emm_list`](https://rvlenth.github.io/emmeans/reference/rbind.emmGrid.md)
and
[`as.data.frame.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)

## Note

The `plot` method uses only the first element of `which`; the others are
ignored.

No `export` option is provided for printing an `emm_list` (see
[`print.emmGrid`](https://rvlenth.github.io/emmeans/reference/emmGrid-methods.md)).
If you wish to export these objects, you must do so separately for each
element in the list.

## Examples

``` r
mod <- lm(conc ~ source, data = pigs)
obj <- emmeans(mod, pairwise ~ source)

linfct(obj)
#> $emmeans
#>      (Intercept) sourcesoy sourceskim
#> [1,]           1         0          0
#> [2,]           1         1          0
#> [3,]           1         0          1
#> 
#> $contrasts
#>      (Intercept) sourcesoy sourceskim
#> [1,]           0        -1          0
#> [2,]           0         0         -1
#> [3,]           0         1         -1
#> 

coef(obj)     # done only for the contrasts
#>   source c.1 c.2 c.3
#> 1   fish   1   1   0
#> 2    soy  -1   0   1
#> 3   skim   0  -1  -1

contrast(obj, "consec")  # done only for the means
#>  contrast   estimate   SE df t.ratio p.value
#>  soy - fish     7.98 2.85 26   2.795  0.0181
#>  skim - soy     5.80 2.93 26   1.979  0.1047
#> 
#> P value adjustment: mvt method for 2 tests 

contrast(obj, "eff", drop = FALSE)   # kept as a list
#>  contrast    estimate   SE df t.ratio p.value
#>  fish effect   -7.255 1.66 26  -4.362  0.0005
#>  soy effect     0.725 1.66 26   0.436  0.6664
#>  skim effect    6.530 1.71 26   3.823  0.0011
#> 
#> P value adjustment: fdr method for 3 tests 
```
