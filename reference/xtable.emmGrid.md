# Using `xtable` for EMMs

These methods provide support for the xtable package, enabling polished
presentations of tabular output from
[`emmeans`](https://rvlenth.github.io/emmeans/reference/emmeans.md) and
other functions.

## Usage

``` r
# S3 method for class 'emmGrid'
xtable(x, caption = NULL, label = NULL, align = NULL,
  digits = 4, display = NULL, auto = FALSE, ...)

# S3 method for class 'summary_emm'
xtable(x, caption = NULL, label = NULL,
  align = NULL, digits = 4, display = NULL, auto = FALSE, ...)

# S3 method for class 'xtable_emm'
print(x, type = getOption("xtable.type", "latex"),
  include.rownames = FALSE, sanitize.message.function = footnotesize, ...)
```

## Arguments

- x:

  Object of class `emmGrid`

- caption:

  Passed to
  [`xtableList`](https://rdrr.io/pkg/xtable/man/xtableList.html)

- label:

  Passed to `xtableList`

- align:

  Passed to `xtableList`

- digits:

  Passed to `xtableList`

- display:

  Passed to `xtableList`

- auto:

  Passed to `xtableList`

- ...:

  Arguments passed to
  [`summary.emmGrid`](https://rvlenth.github.io/emmeans/reference/summary.emmGrid.md)

- type:

  Passed to
  [`print.xtable`](https://rdrr.io/pkg/xtable/man/print.xtable.html)

- include.rownames:

  Passed to `print.xtable`

- sanitize.message.function:

  Passed to `print.xtable`

## Value

The `xtable` methods return an `xtable_emm` object, for which its print
method is `print.xtable_emm` .

## Details

The methods actually use
[`xtableList`](https://rdrr.io/pkg/xtable/man/xtableList.html), because
of its ability to display messages such as those for P-value
adjustments. These methods return an object of class `"xtable_emm"` – an
extension of `"xtableList"`. Unlike other `xtable` methods, the number
of digits defaults to 4; and degrees of freedom and *t* ratios are
always formatted independently of `digits`. The `print` method uses
[`print.xtableList`](https://rdrr.io/pkg/xtable/man/xtableList.html),
and any `...` arguments are passed there.

## Examples

``` r
if(requireNamespace("xtable"))
    emm_example("xtable")
#> 
#> --- Running code from 'system.file("extexamples", "xtable.R", package = "emmeans")'
#> 
#> > pigsint.lm <- lm(log(conc) ~ source * factor(percent), 
#> +     data = pigs)
#> 
#> > pigsint.emm <- emmeans(pigsint.lm, ~percent | source)
#> 
#> > xtable::xtable(pigsint.emm, type = "response")
#> % latex table generated in R 4.5.2 by xtable 1.8-8 package
#> % Fri Feb 27 00:58:50 2026
#> \begin{table}[ht]
#> \centering
#> \begin{tabular}{rrrrrr}
#>   \hline
#> percent & response & SE & df & lower.CL & upper.CL \\ 
#>   \hline
#> \multicolumn{6}{l}{source = fish}\\
#> 9.0000 & 25.6683 & 2.1101 & 17 & 21.5810 & 30.5296 \\ 
#>   12.0000 & 30.8799 & 2.0727 & 17 & 26.8025 & 35.5777 \\ 
#>   15.0000 & 31.0193 & 2.5500 & 17 & 26.0801 & 36.8941 \\ 
#>   18.0000 & 32.3072 & 2.1685 & 17 & 28.0413 & 37.2222 \\ 
#>    \hline
#> \multicolumn{6}{l}{source = soy}\\
#> 9.0000 & 34.4135 & 2.3099 & 17 & 29.8695 & 39.6489 \\ 
#>   12.0000 & 39.6314 & 2.6601 & 17 & 34.3984 & 45.6606 \\ 
#>   15.0000 & 39.2286 & 2.6331 & 17 & 34.0487 & 45.1964 \\ 
#>   18.0000 & 42.9000 & 4.9874 & 17 & 33.5686 & 54.8254 \\ 
#>    \hline
#> \multicolumn{6}{l}{source = skim}\\
#> 9.0000 & 35.1821 & 2.3615 & 17 & 30.5365 & 40.5343 \\ 
#>   12.0000 & 43.1574 & 2.8968 & 17 & 37.4588 & 49.7230 \\ 
#>   15.0000 & 49.6316 & 4.0800 & 17 & 41.7287 & 59.0314 \\ 
#>   18.0000 & 59.8000 & 6.9522 & 17 & 46.7926 & 76.4232 \\ 
#>    \hline
#> \multicolumn{6}{l}{{\footnotesize Confidence level used: 0.95}}\\
#> 
#> \multicolumn{6}{l}{{\footnotesize Intervals are back-transformed from the log scale}}\\
#> \end{tabular}
#> \end{table}
#> 
    # Use emm_example("xtable", list = TRUE) # to just list the code
```
