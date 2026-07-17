# Unbalanced dataset

This is a simulated unbalanced dataset with three factors and two
numeric variables. There are true relationships among these variables.
This dataset can be useful in testing or illustrating messy-data
situations. There are no missing data, and there is at least one
observation for every factor combination; however, the `"cells"`
attribute makes it simple to construct subsets that have empty cells.

## Usage

``` r
ubds
```

## Format

A data frame with 100 observations, 5 variables, and a special `"cells"`
attribute:

- A:

  Factor with levels 1, 2, and 3

- B:

  Factor with levels 1, 2, and 3

- C:

  Factor with levels 1, 2, and 3

- x:

  A numeric variable

- y:

  A numeric variable

In addition, `attr(ubds, "cells")` consists of a named list of length 27
with the row numbers for each combination of `A, B, C`. For example,
`attr(ubds, "cells")[["213"]]` has the row numbers corresponding to
levels `A == 2, B == 1, C == 3`. The entries are ordered by length, so
the first entry is the cell with the lowest frequency.

## Examples

``` r
 # Omit the three lowest-frequency cells
 low3 <- unlist(attr(ubds, "cells")[1:3]) 
 messy.lm <- lm(y ~ (x + A + B + C)^3, data = ubds, subset = -low3)
  
```
