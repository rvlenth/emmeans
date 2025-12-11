# Dare to be un-"tidy"!

Users who use emmeans functions as part of a pipeline – or post-process
those results in some other way – are likely missing some important
information.

## Details

Your best bet is to display the actual results without any
post-processing. That's because `emmeans` and its relatives have their
own `summary` and `print` methods that display annotations that may be
helpful in explaining what you have. If you just pipe the results into
the next step, those annotations are stripped away and you never see
them. Statistical analysis is not just a workflow; it is a discipline
that involves care in interpreting intermediate results, and thinking
before moving on.

## Examples

``` r
neur.glm <- glm(Pain ~ Treatment + Sex + Age, family = binomial(),
            data = neuralgia)
            
### The actual results with annotations (e.g. ests are on logit scale):
emmeans(neur.glm, "Treatment")
#>  Treatment emmean    SE  df asymp.LCL asymp.UCL
#>  A          -1.40 0.664 Inf    -2.699   -0.0951
#>  B          -1.94 0.744 Inf    -3.403   -0.4867
#>  P           1.78 0.683 Inf     0.444    3.1198
#> 
#> Results are averaged over the levels of: Sex 
#> Results are given on the logit (not the response) scale. 
#> Confidence level used: 0.95 

### Post-processed results lose the annotations
if(requireNamespace("tibble")) {
    emmeans(neur.glm, "Treatment") |> tibble::as_tibble()
}
#> # A tibble: 3 × 6
#>   Treatment emmean    SE    df asymp.LCL asymp.UCL
#>   <fct>      <dbl> <dbl> <dbl>     <dbl>     <dbl>
#> 1 A          -1.40 0.664   Inf    -2.70    -0.0951
#> 2 B          -1.94 0.744   Inf    -3.40    -0.487 
#> 3 P           1.78 0.683   Inf     0.444    3.12  
```
