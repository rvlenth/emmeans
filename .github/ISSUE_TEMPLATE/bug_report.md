---
name: Bug report
about: Create a report to help us improve
title: ''
labels: bug
assignees: ''

---
## Please ensure that you are in the right place
Some model classes have their **emmeans** support in another package -- 
e.g., the one that defines the model class itself. If so, is the
bug really in that package, rather than in **emmeans**?

## Describe the bug
A clear and concise description of what the bug is.

## To reproduce
Show me code and output that reproduces the bug. 
Please use a *small data set* (built-in if possible) and *simple variable names*.
And please, do not create different objects having the same name; 
to compare two or three different models or methods, give those objects
different names so we can talk about them.

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
    steps and the results thereof. 
  * Did you know that there is an index of vignette topics? That can be
    helpful for clarifying certain issues. See
    https://cran.r-project.org/web/packages/emmeans/vignettes/vignette-topics.html
  * Please examine the output from what you have tried. Often there are a few
    lines of annotation below the output. If you wrap your results with 
    `as.data.frame()` once or twice to suppress those annotations, I
    consider that willful ignorance, and will not help you.
