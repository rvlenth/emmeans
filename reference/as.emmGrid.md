# Convert to and from `emmGrid` objects

These are useful utility functions for creating a compact version of an
`emmGrid` object that may be saved and later reconstructed, or for
converting old `ref.grid` or `lsmobj` objects into `emmGrid` objects.

## Usage

``` r
# S3 method for class 'emmGrid'
as.list(x, model.info.slot = FALSE, ...)

as.emm_list(object, ...)

as.emmGrid(object, ...)
```

## Arguments

- x:

  An `emmGrid` object

- model.info.slot:

  Logical value: Include the `model.info` slot? Set this to `TRUE` if
  you want to preserve the original call and information needed by the
  `submodel` option. If `FALSE`, only the nesting information (if any)
  is saved

- ...:

  In `as.emmGrid`, additional arguments passed to
  [`update.emmGrid`](https://rvlenth.github.io/emmeans/reference/update.emmGrid.md)
  before returning the object. This argument is ignored in
  `as.list.emmGrid`

- object:

  Object to be converted to class `emmGrid`. It may be a `list` returned
  by `as.list.emmGrid`, or a `ref.grid` or `lsmobj` object created by
  emmeans's predecessor, the lsmeans package. An error is thrown if
  `object` cannot be converted.

## Value

`as.list.emmGrid` returns an object of class `list`.

`as.emm_list` returns an object of class `emm_list`.

`as.emmGrid` returns an object of class `emmGrid`. However, in fact,
both `as.emmGrid` and `as.emm_list` check for an attribute in `object`
to decide whether to return an `emmGrid` or `emm_list)` object.

## Details

An `emmGrid` object is an S4 object, and as such cannot be saved in a
text format or saved without a lot of overhead. By using `as.list`, the
essential parts of the object are converted to a list format that can be
easily and compactly saved for use, say, in another session or by
another user. Providing this list as the arguments for
[`emmobj`](https://rvlenth.github.io/emmeans/reference/emmobj.md) allows
the user to restore a working `emmGrid` object.

## See also

[`emmobj`](https://rvlenth.github.io/emmeans/reference/emmobj.md)

## Examples

``` r
pigs.lm <- lm(log(conc) ~ source + factor(percent), data = pigs)
pigs.sav <- as.list(ref_grid(pigs.lm))

pigs.anew <- as.emmGrid(pigs.sav)
emmeans(pigs.anew, "source")
#>  source emmean     SE df lower.CL upper.CL
#>  fish     3.39 0.0367 23     3.32     3.47
#>  soy      3.67 0.0374 23     3.59     3.74
#>  skim     3.80 0.0394 23     3.72     3.88
#> 
#> Results are averaged over the levels of: percent 
#> Results are given on the log (not the response) scale. 
#> Confidence level used: 0.95 
```
