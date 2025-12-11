# Effects of dietary protein on free plasma leucine concentration in pigs

A two-factor experiment with some observations lost

## Usage

``` r
pigs
```

## Format

A data frame with 29 observations and 3 variables:

- source:

  Source of protein in the diet (factor with 3 levels: fish meal,
  soybean meal, dried skim milk)

- percent:

  Protein percentage in the diet (numeric with 4 values: 9, 12, 15, and
  18)

- conc:

  Concentration of free plasma leucine, in mcg/ml

## Source

Windels HF (1964) PhD thesis, Univ. of Minnesota. (Reported as Problem
10.8 in Oehlert G (2000) *A First Course in Design and Analysis of
Experiments*, licensed under Creative Commons,
[http://users.stat.umn.edu/~gary/Book.html](http://users.stat.umn.edu/~gary/Book.md).)
Observations 7, 22, 23, 31, 33, and 35 have been omitted, creating a
more notable imbalance.

## Examples

``` r
  pigs.lm <- lm(inverse(conc) ~ source + factor(percent), data = pigs)
  emmeans(pigs.lm, "source")
#>  source emmean       SE df lower.CL upper.CL
#>  fish   0.0337 0.000926 23   0.0318   0.0356
#>  soy    0.0257 0.000945 23   0.0237   0.0276
#>  skim   0.0229 0.000994 23   0.0208   0.0249
#> 
#> Results are averaged over the levels of: percent 
#> Results are given on the inverse (not the response) scale. 
#> Confidence level used: 0.95 
```
