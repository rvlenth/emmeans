---
name: Question
about: I have a question about **emmeans**
title: ''
labels: question
assignees: ''

---

## Please ensure that you are in the right place
Some model classes have their **emmeans** support in another package -- 
e.g., the one that defines the model class itself. If so, you should
consider directing your question to that package, rather than to **emmeans**.

## Please show me the *actual* results from **emmeans** function(s)
Often there are useful annotation that you'll miss if you post-process your code.
Seeing the actual results may even answer your question!
``` r
# --- Actual results (YES):
emmeans(mod, ~ Sex)
##  Sex emmean    SE  df asymp.LCL asymp.UCL
##  F   -1.432 0.582 Inf    -2.572    -0.291
##  M    0.392 0.492 Inf    -0.573     1.356
## 
## Results are averaged over the levels of: Treatment 
## Results are given on the logit (not the response) scale. 
## Confidence level used: 0.95

# --- Filtered results (NO):
emmeans(mod, ~ Sex) %>% as_tibble()
## # A tibble: 2 × 6
##   Sex   emmean    SE    df asymp.LCL asymp.UCL
##   <fct>  <dbl> <dbl> <dbl>     <dbl>     <dbl>
## 1 F     -1.43  0.582   Inf    -2.57     -0.291
## 2 M      0.392 0.492   Inf    -0.573     1.36
```

## Explain your question
Tell me your concern, and illustrate it using enough code and output that I 
can reproduce it, and enough output (including annotations and messages) 
that I can tell what went wrong. Try to keep the example small and the
variable names simple.


## Ground rules
  * I really do expect you to look at the documentation and vignettes before
    sending questions.
  * Did you know that there is an index of vignette topics? That can be
    helpful for finding help on certain topics. See
    https://cran.r-project.org/web/packages/emmeans/vignettes/vignette-topics.html
  * More than one or two pipes is usually too many. I'd rather see the individual 
    steps and the results thereof.
  * Please examine the *direct* output from what you have tried. Often there are a few
    lines of annotation below the output. If you wrap or pipe your results with
    some kind of post-processing, you suppress those annotations. I consider
    that willful ignorance, and will not help you.
  * Please do not create different objects having the same name; 
    that causes confusion and makes it hard to refer to particular results.
