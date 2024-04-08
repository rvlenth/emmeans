---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---
Please replace the descriptions below with your content that relates to that topic.

## Please ensure that you are in the right place
Some model classes have their **emmeans** support in another package -- 
e.g., the one that defines the model class itself. If so, is the
bug really in that package, rather than in **emmeans**?

## Please show me the *actual* results from **emmeans** function(s)
Often there are useful annotation that you'll miss if you post-process your code.
Seeing the actual results may even explain what's happening!
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


## Describe the bug
A clear and concise description of what the bug is.

## To reproduce
Show me code and output that reproduces the bug. 
Please use a *small data set* (built-in if possible) and *simple variable names*.
And please, do not create different objects having the same name; 
to compare two or three different models or methods, give those objects
different names so we can talk about them.

**Show the actual output** from whatever function is in question -- do not
pipe it into some post-processing stuff. The actual output often shows some annotations that you don't see if you are trying to "tidy" it. Nothing would make me happier than to
avoid *all* tidyverse code.

## Expected behavior
A clear and concise description of what you expected to happen.

## Additional context
Add any other context about the problem here.

## Ground rules
  * I really do expect you to look at the documentation and vignettes before
    sending bug reports. Make sure that you have used things as they are documented.
  * I really do not want to see your whole workflow. Just show me the code for
    fitting the model in question (yes, *all* of the code for that including
    what libraries are needed), and complete output for where the bug occurs.
    Leave out things like plots and code/output from other packages unless
    it is relevant to reproducing the bug.
  * Show output as pre-formatted text, not a graphic screen shot.
  * Again, **do not ever re-use object names.** If you are comparing results for
    two or more models or datasets, assign them different names.
  * More than one or two pipes is usually too many. I'd rather see the individual 
    steps and the actual, unadorned results thereof. 
  * Did you know that there is an index of vignette topics? That can be
    helpful for clarifying certain issues. See
    https://cran.r-project.org/web/packages/emmeans/vignettes/vignette-topics.html
  * Please examine the *direct* output from what you have tried. Often there are a few
    lines of annotation below the output. If you wrap or pipe your results with
    some kind of post-processing, you suppress those annotations. I consider
    that willful ignorance, and will not help you.
