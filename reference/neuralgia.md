# Neuralgia data

These data arise from a study of analgesic effects of treatments of
elderly patients who have neuralgia. Two treatments and a placebo are
compared. The response variable is whether the patient reported pain or
not. Researchers recorded the age and gender of 60 patients along with
the duration of complaint before the treatment began.

## Usage

``` r
neuralgia
```

## Format

A data frame with 60 observations and 5 variables:

- `Treatment`:

  Factor with 3 levels `A`, `B`, and `P`. The latter is placebo

- `Sex`:

  Factor with two levels `F` and `M`

- `Age`:

  Numeric covariate – patient's age in years

- `Duration`:

  Numeric covariate – duration of the condition before beginning
  treatment

- `Pain`:

  Binary response factor with levels `No` and `Yes`

## Source

Cai, Weijie (2014) *Making Comparisons Fair: How LS-Means Unify the
Analysis of Linear Models*, SAS Institute, Inc. Technical paper
142-2014, page 12,
<http://support.sas.com/resources/papers/proceedings14/SAS060-2014.pdf>

## Examples

``` r
# Model and analysis shown in the SAS report:
neuralgia.glm <- glm(Pain ~ Treatment * Sex + Age, family = binomial(),
   data = neuralgia) 
pairs(emmeans(neuralgia.glm, ~ Treatment, at = list(Sex = "F")), 
    reverse = TRUE, type = "response", adjust = "bonferroni")
#> NOTE: Results may be misleading due to involvement in interactions
#>  contrast odds.ratio     SE  df null z.ratio p.value
#>  B / A         0.398  0.648 Inf    1  -0.566  1.0000
#>  P / A        16.892 22.300 Inf    1   2.141  0.0969
#>  P / B        42.492 63.400 Inf    1   2.511  0.0361
#> 
#> P value adjustment: bonferroni method for 3 tests 
#> Tests are performed on the log odds ratio scale 
```
