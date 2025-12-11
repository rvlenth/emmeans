# Compact letter displays

A method for
[`multcomp::cld()`](https://rdrr.io/pkg/multcomp/man/cld.html) is
provided for users desiring to produce compact-letter displays (CLDs).
This method uses the Piepho (2004) algorithm (as implemented in the
multcompView package) to generate a compact letter display of all
pairwise comparisons of estimated marginal means. The function obtains
(possibly adjusted) P values for all pairwise comparisons of means,
using the
[`contrast`](https://rvlenth.github.io/emmeans/reference/contrast.md)
function with `method = "pairwise"`. When a P value exceeds `alpha`,
then the two means have at least one letter in common.

## Usage

``` r
# S3 method for class 'emmGrid'
cld(object, details = FALSE, sort = TRUE, by,
  alpha = 0.05, Letters = c("1234567890", LETTERS, letters),
  reversed = decreasing, decreasing = FALSE, signif.sets = FALSE,
  delta = 0, ...)

# S3 method for class 'emm_list'
cld(object, ..., which = 1)
```

## Arguments

- object:

  An object of class `emmGrid`

- details:

  Logical value determining whether detailed information on tests of
  pairwise comparisons is displayed

- sort:

  Logical value determining whether the EMMs are sorted before the
  comparisons are produced. When `TRUE`, the results are displayed
  according to `reversed`.

- by:

  Character value giving the name or names of variables by which
  separate families of comparisons are tested. If NULL, all means are
  compared. If missing, the object's `by.vars` setting, if any, is used.

- alpha:

  Numeric value giving the significance level for the comparisons

- Letters:

  Character vector of letters to use in the display. Any strings of
  length greater than 1 are expanded into individual characters

- reversed, decreasing:

  Logical value (passed to
  [`multcompView::multcompLetters`](https://rdrr.io/pkg/multcompView/man/multcompLetters.html).)
  If `TRUE`, the order of use of the letters is reversed. Either
  `reversed` or `decreasing` may be specified, thus providing
  compatibility with both
  `multcompView::multcompLetters(..., reversed, ...)` and
  `multcomp::cld(..., decreasing, ...)`. In addition, if both `sort` and
  `reversed` are TRUE, the sort order of results is reversed.

- signif.sets:

  Logical value. If `FALSE` (and `delta = 0`), a ‘traditional’
  compact-letter display is constructed with groupings representing sets
  of estimates that are not statistically different. If `TRUE`, the
  criteria are reversed so that two estimates sharing the same symbol
  test as significantly different. See also `delta`.

- delta:

  Numeric value passed to
  [`test.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md).
  If this is positive, it is used as an equivalence threshold in the
  TOST procedure for two-sided equivalence testing. In the resulting
  compact letter display, two estimates share the same grouping letter
  only if they are found to be statistically equivalent – that is,
  groupings reflect actual *findings* of equivalence rather than failure
  to find a significant difference. When `delta` is nonzero,
  `signif.sets` is ignored.

- ...:

  Arguments passed to
  [`contrast`](https://rvlenth.github.io/emmeans/reference/contrast.md)
  (for example, an `adjust` method)

- which:

  Which element of the `emm_list` object to process (If length exceeds
  one, only the first one is used)

## Value

A
[`summary_emm`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)
object showing the estimated marginal means plus an additional column
labeled `.group` (when `signif.sets = FALSE`), `.signif.set` (when
`signif.sets = TRUE`), or `.equiv.set` (when `delta > 0`).

## Note

We warn that the default display encourages a poor practice in
interpreting significance tests. Such CLDs are misleading because they
visually group means with comparisons *P* \> `alpha` as though they are
equal, when in fact we have only failed to prove that they differ. A
better alternative if one wants to show groupings is to specify an
equivalence threshold `delta`; then groupings will be based on actual
findings of equivalence. Another way to display actual findings is to
set `signif.sets = TRUE`, so that estimates in the same group are those
found to be statistically *different*. Obviously, these different
options require different interpretations of the results; the
annotations and the label given the final column help guide how to
assess the results.

As further alternatives, consider
[`pwpp`](https://rvlenth.github.io/emmeans/reference/pwpp.md) (graphical
display of *P* values) or
[`pwpm`](https://rvlenth.github.io/emmeans/reference/pwpm.md) (matrix
display).

## References

Piepho, Hans-Peter (2004) An algorithm for a letter-based representation
of all pairwise comparisons, Journal of Computational and Graphical
Statistics, 13(2), 456-466.

## Examples

``` r
if(requireNamespace("multcomp"))
    emm_example("cld-multcomp")
#> 
#> --- Running code from 'system.file("extexamples", "cld-multcomp.R", package = "emmeans")'
#> 
#> > pigs.lm <- lm(log(conc) ~ source + factor(percent), 
#> +     data = pigs)
#> 
#> > pigs.emm <- emmeans(pigs.lm, "percent", type = "response")
#> 
#> > multcomp::cld(pigs.emm, alpha = 0.1, Letters = LETTERS)
#>  percent response   SE df lower.CL upper.CL .group
#>        9     31.4 1.28 23     28.8     34.1  A    
#>       12     37.5 1.44 23     34.7     40.6   B   
#>       15     39.0 1.70 23     35.6     42.7   B   
#>       18     42.3 2.24 23     37.9     47.2   B   
#> 
#> Results are averaged over the levels of: source 
#> Confidence level used: 0.95 
#> Intervals are back-transformed from the log scale 
#> P value adjustment: tukey method for comparing a family of 4 estimates 
#> Tests are performed on the log scale 
#> significance level used: alpha = 0.1 
#> NOTE: If two or more means share the same grouping symbol,
#>       then we cannot show them to be different.
#>       But we also did not show them to be the same. 
#> 
#> > multcomp::cld(pigs.emm, alpha = 0.1, signif.sets = TRUE)
#>  percent response   SE df lower.CL upper.CL .signif.set
#>        9     31.4 1.28 23     28.8     34.1  123       
#>       12     37.5 1.44 23     34.7     40.6  1         
#>       15     39.0 1.70 23     35.6     42.7   2        
#>       18     42.3 2.24 23     37.9     47.2    3       
#> 
#> Results are averaged over the levels of: source 
#> Confidence level used: 0.95 
#> Intervals are back-transformed from the log scale 
#> P value adjustment: tukey method for comparing a family of 4 estimates 
#> Tests are performed on the log scale 
#> significance level used: alpha = 0.1 
#> Estimates sharing the same symbol are significantly different 
#> 
#> > multcomp::cld(pigs.emm, delta = log(1.25), adjust = "sidak")
#>  percent response   SE df lower.CL upper.CL .equiv.set
#>        9     31.4 1.28 23     28.1     35.0  1        
#>       12     37.5 1.44 23     33.8     41.6   2       
#>       15     39.0 1.70 23     34.6     43.9   2       
#>       18     42.3 2.24 23     36.7     48.8    3      
#> 
#> Results are averaged over the levels of: source 
#> Confidence level used: 0.95 
#> Conf-level adjustment: sidak method for 4 estimates 
#> Intervals are back-transformed from the log scale 
#> P value adjustment: sidak method for 6 tests 
#> Statistics are tests of equivalence with a threshold of 0.22314 
#> P values are left-tailed 
#> Tests are performed on the log scale 
#> significance level used: alpha = 0.05 
#> Estimates sharing the same symbol test as equivalent 
#> 
    # Use emm_example("cld-multcomp", list = TRUE) # to just list the code
```
