library(testthat)
library(emmeans)

# I need to write a lot more tests. Those will come.
# We do have a ton of examples and vignettes, so a bunch of
# stuff DOES get tested more informally.

# I'm actually more concerned about testing support for all
# classes of model objects, whether in Suggests or not. Those tests
# will largely have to be in .Rbuildignore

test_check("emmeans")
