##############################################################################
#    Copyright (c) 2012-2017 Russell V. Lenth                                #
#                                                                            #
#    This file is part of the emmeans package for R (*emmeans*)              #
#                                                                            #
#    *emmeans* is free software: you can redistribute it and/or modify       #
#    it under the terms of the GNU General Public License as published by    #
#    the Free Software Foundation, either version 2 of the License, or       #
#    (at your option) any later version.                                     #
#                                                                            #
#    *emmeans* is distributed in the hope that it will be useful,            #
#    but WITHOUT ANY WARRANTY; without even the implied warranty of          #
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           #
#    GNU General Public License for more details.                            #
#                                                                            #
#    You should have received a copy of the GNU General Public License       #
#    along with R and *emmeans*.  If not, see                                #
#    <https://www.r-project.org/Licenses/> and/or                            #
#    <http://www.gnu.org/licenses/>.                                         #
##############################################################################

# datasets (provided as .rda files in data/) -- this file is for documentation

### auto.noise ###
#' Auto Pollution Filter Noise
#' 
#' Three-factor experiment comparing pollution-filter noise for two filters,
#' three sizes of cars, and two sides of the car.
#' 
#' The data are from a statement by Texaco, Inc., to the Air and Water Pollution
#' Subcommittee of the Senate Public Works Committee on June 26, 1973.    
#' Mr. John McKinley, President of Texaco, cited an automobile filter developed
#' by Associated Octel Company as effective in reducing pollution. However, 
#' questions had been raised about the effects of filters on vehicle performance, 
#' fuel consumption, exhaust gas back pressure, and silencing. On the last 
#' question, he referred to the data included here as evidence that the silencing
#' properties of the Octel filter were at least equal to those of standard silencers.
#'  
#' @format A data frame with 36 observations on the following 4 variables.
#'   \describe{
#'   \item{\code{noise}}{Noise level in decibels (but see note) - a numeric vector.}
#'   \item{\code{size}}{The size of the vehicle - an ordered factor with
#'     levels \code{S}, \code{M}, \code{L}.}
#'   \item{\code{type}}{Type of anti-pollution filter - a factor with levels
#'     \code{Std} and \code{Octel}}
#'   \item{\code{side}}{The side of the car where measurement was taken -- a
#'     factor with levels \code{L} and \code{R}.}
#'   }
#' @source The dataset was obtained from the Data and Story Library (DASL)
#'   at Carnegie-Mellon University. Apparently it has since been removed. The
#'   original dataset was altered by assigning meaningful names to the factors
#'   and sorting the observations in random order as if this were the run order
#'   of the experiment.
#' @note While the data source claims that \code{noise} is measured in decibels,
#'   the values are implausible. I believe that these measurements are actually
#'   in tenths of dB (centibels?). Looking at the values in the dataset, note 
#'   that every measurement ends in 0 or 5, and it is reasonable to believe that
#'   measurements are accurate to the nearest half of a decibel.
#'   %%% Thanks to an email communication from a speech/hearing scientist
#' @examples
#' # (Based on belief that noise/10 is in decibel units)
#' noise.lm <- lm(noise/10 ~ size * type * side, data = auto.noise)
#' 
#' # Interaction plot of predictions
#' emmip(noise.lm, type ~ size | side)
#' 
#' # Confidence intervals
#' plot(emmeans(noise.lm, ~ size | side*type))
#' 
"auto.noise"

# This is where it used to be...
#  \url{http://lib.stat.cmu.edu/DASL/Datafiles/airpullutionfiltersdat.html}


### feedlot ###
#' Feedlot data
#' 
#' This is an unbalanced analysis-of-covariance example, where one covariate is
#' affected by a factor. Feeder calves from various herds enter a feedlot, where
#' they are fed one of three diets. The weight of the animal at entry is the
#' covariate, and the weight at slaughter is the response.
#' 
#' The data arise from a Western Regional Research Project conducted at New
#' Mexico State University. Calves born in 1975 in commercial herds entered a
#' feedlot as yearlings. Both diets and herds are of interest as factors. The
#' covariate, \code{ewt}, is thought to be dependent on \code{herd} due to
#' different genetic backgrounds, breeding history, etc. The levels of
#' \code{herd} ordered to similarity of genetic background.
#' 
#' Note: There are some empty cells in the cross-classification of 
#' \code{herd} and \code{diet}.
#' @format A data frame with 67 observations and 4 variables:
#' \describe{
#'   \item{\code{herd}}{a factor with levels \code{9} \code{16} \code{3}
#'     \code{32} \code{24} \code{31} \code{19} \code{36} \code{34} \code{35}
#'     \code{33}, designating the herd that a feeder calf came from.}
#'   \item{\code{diet}}{a factor with levels \code{Low} \code{Medium}
#'     \code{High}: the energy level of the diet given the animal.}
#'   \item{\code{swt}}{a numeric vector: the weight of the animal at slaughter.}
#'   \item{\code{ewt}}{a numeric vector: the weight of the animal at entry to the feedlot.}
#' }
#' @source Urquhart NS (1982) Adjustment in covariates when one factor affects
#'   the covariate. \emph{Biometrics} 38, 651-660.
#' @examples
#' feedlot.lm <- lm(swt ~ ewt + herd*diet, data = feedlot)
#' 
#' # Obtain EMMs with a separate reference value of ewt for each 
#' # herd. This reproduces the last part of Table 2 in the reference
#' emmeans(feedlot.lm,  ~ diet | herd,  cov.reduce = ewt ~ herd)
#' 
"feedlot"




### fiber ###
#' Fiber data
#' 
#' Fiber data from Montgomery Design (8th ed.), p.656 (Table 15.10). Useful as a
#' simple analysis-of-covariance example.
#' 
#' The goal of the experiment is to compare the mean breaking strength of fibers
#' produced by the three machines. When testing this, the technician also
#' measured the diameter of each fiber, and this measurement may be used as a
#' concomitant variable to improve precision of the estimates.
#' @format A data frame with 15 observations and 3 variables:
#' \describe{
#'   \item{\code{machine}}{a factor with levels \code{A} \code{B} \code{C}. 
#'     This is the primary factor of interest.}
#'   \item{\code{strength}}{a numeric vector. The response variable.}
#'   \item{\code{diameter}}{a numeric vector. A covariate.}
#' }
#' @source Montgomery, D. C. (2013) \emph{Design and Analysis of Experiments}
#'   (8th ed.). John Wiley and Sons, ISBN 978-1-118-14692-7.
#' @examples
#' fiber.lm <- lm(strength ~ diameter + machine, data=fiber)
#' ref_grid(fiber.lm)
#' 
#' # Covariate-adjusted means and comparisons
#' emmeans(fiber.lm, pairwise ~ machine)
#' 
"fiber"




### MOats ###
#' Oats data in multivariate form
#' 
#' This is the \code{Oats} dataset provided in the \pkg{nlme} package, but it is
#' rearranged as one multivariate observation per plot.
#' 
#' These data arise from a split-plot experiment reported by Yates (1935) and
#' used as an example in Pinheiro and Bates (2000) and other texts. Six blocks
#' were divided into three whole plots, randomly assigned to the three varieties
#' of oats. The whole plots were each divided into 4 split plots and randomized
#' to the four concentrations of nitrogen.
#' @format A data frame with 18 observations and 3 variables
#' \describe{
#'   \item{\code{Variety}}{a factor with levels \code{Golden Rain},
#'     \code{Marvellous}, \code{Victory}}
#'   \item{\code{Block}}{an ordered factor with levels \code{VI} < \code{V} <
#'     \code{III} < \code{IV} < \code{II} < \code{I}}
#'   \item{\code{yield}}{a matrix with 4 columns, giving the yields with
#'     nitrogen concentrations of 0, .2, .4, and .6.}
#' }
#' @source The dataset \code{\link[nlme]{Oats}} in the \pkg{nlme} package.
#' @references
#' Pinheiro, J. C. and Bates D. M. (2000) \emph{Mixed-Effects Models in S and
#' S-PLUS}, Springer, New York. (Appendix A.15)
#' 
#' Yates, F. (1935) Complex experiments, \emph{Journal of the Royal Statistical
#' Society} Suppl. 2, 181-247
#' @examples
#' MOats.lm <- lm (yield ~ Block + Variety, data = MOats)
#' MOats.rg <- ref_grid (MOats.lm, mult.name = "nitro")
#' emmeans(MOats.rg, ~ nitro | Variety)
"MOats"



### neuralgia ###
#' Neuralgia data
#' 
#' These data arise from a study of analgesic effects of treatments of elderly
#' patients who have neuralgia. Two treatments and a placebo are compared. The 
#' response variable is whether the patient reported pain or not. Researchers
#' recorded the age and gender of 60 patients along with the duration of
#' complaint before the treatment began.
#' 
#' @format A data frame with 60 observations and 5 variables:
#' \describe{
#'   \item{\code{Treatment}}{Factor with 3 levels \code{A}, \code{B}, and \code{P}.
#'     The latter is placebo}
#'   \item{\code{Sex}}{Factor with two levels \code{F} and \code{M}}
#'   \item{\code{Age}}{Numeric covariate -- patient's age in years}
#'   \item{\code{Duration}}{Numeric covariate -- duration of the condition before
#'     beginning treatment}
#'   \item{\code{Pain}}{Binary response factor with levels \code{No} and \code{Yes}}
#' }
#' @source Cai, Weijie (2014) \emph{Making Comparisons Fair: How LS-Means Unify 
#'   the Analysis of Linear Models}, SAS Institute, Inc. Technical paper 142-2014,
#'   page 12, 
#'   \url{http://support.sas.com/resources/papers/proceedings14/SAS060-2014.pdf}
#' @examples
#' # Model and analysis shown in the SAS report:
#' neuralgia.glm <- glm(Pain ~ Treatment * Sex + Age, family = binomial(),
#'    data = neuralgia) 
#' pairs(emmeans(neuralgia.glm, ~ Treatment, at = list(Sex = "F")), 
#'     reverse = TRUE, type = "response", adjust = "bonferroni")
#' 
"neuralgia"

### nutrition ###
#' Nutrition data
#' 
#' This observational dataset involves three factors, but where several factor 
#' combinations are missing. It is used as a case study in Milliken and Johnson,
#' Chapter 17, p.202. (You may also find it in the second edition, p.278.)
#' 
#' A survey was conducted by home economists ``to study how much
#' lower-socioeconomic-level mothers knew about nutrition and to judge the
#' effect of a training program designed to increase their knowledge of
#' nutrition.'' This is a messy dataset with several empty cells.
#' @format A data frame with 107 observations and 4 variables:
#' \describe{
#'   \item{\code{age}}{a factor with levels \code{1}, \code{2}, \code{3},
#'     \code{4}. Mother's age group.}
#'   \item{\code{group}}{a factor with levels \code{FoodStamps}, \code{NoAid}.
#'     Whether or not the family receives food stamp assistance.}
#'   \item{\code{race}}{a factor with levels \code{Black}, \code{Hispanic},
#'     \code{White}. Mother's race.}
#'   \item{\code{gain}}{a numeric vector (the response variable). Gain score
#'     (posttest minus pretest) on knowledge of nutrition.}
#' }
#' @source Milliken, G. A. and Johnson, D. E. (1984)
#' \emph{Analysis of Messy Data -- Volume I: Designed Experiments}. 
#' Van Nostrand, ISBN 0-534-02713-7.
#' @examples
#' nutr.aov <- aov(gain ~ (group + age + race)^2, data = nutrition)
#' 
#' # Summarize predictions for age group 3
#' nutr.emm <- emmeans(nutr.aov, ~ race * group, at = list(age="3"))
#'                    
#' emmip(nutr.emm, race ~ group)
#' 
#' # Hispanics seem exceptional; but this doesn't test out due to very sparse data
#' pairs(nutr.emm, by = "group")
#' pairs(nutr.emm, by = "race")
"nutrition"




### oranges ###
#' Sales of oranges
#' 
#' This example dataset on sales of oranges has two factors, two covariates, and
#' two responses. There is one observation per factor combination.
#' @format A data frame with 36 observations and 6 variables:
#' \describe{
#'   \item{\code{store}}{a factor with levels \code{1} \code{2} \code{3}
#'     \code{4} \code{5} \code{6}. The store that was observed.}
#'   \item{\code{day}}{a factor with levels \code{1} \code{2} \code{3}
#'     \code{4} \code{5} \code{6}. The day the observation was taken (same for
#'     each store).}
#'   \item{\code{price1}}{a numeric vector. Price of variety 1.}
#'   \item{\code{price2}}{a numeric vector. Price of variety 2.}
#'   \item{\code{sales1}}{a numeric vector. Sales (per customer) of variety 1.}
#'   \item{\code{sales2}}{a numeric vector. Sales (per customer) of variety 2.}
#' }
#' @source This is (or once was) available as a SAS sample dataset. 
#' @references
#' Littell, R., Stroup W., Freund, R. (2002) \emph{SAS For Linear Models} (4th
#' edition). SAS Institute. ISBN 1-59047-023-0.
#' @examples
#' # Example on p.244 of Littell et al.
#' oranges.lm <- lm(sales1 ~ price1*day, data = oranges)
#' emmeans(oranges.lm, "day")
#' 
#' # Example on p.246 of Littell et al.
#' emmeans(oranges.lm, "day", at = list(price1 = 0))
#' 
#' # A more sensible model to consider, IMHO (see vignette("interactions"))
#' org.mlm <- lm(cbind(sales1, sales2) ~ price1 * price2 + day + store, 
#'               data = oranges)
"oranges"




### pigs ###
#' Effects of dietary protein on free plasma leucine concentration in pigs
#' 
#' A two-factor experiment with some observations lost
#' 
#' @format A data frame with 29 observations and 3 variables:
#' \describe{
#'   \item{source}{Source of protein in the diet (factor with 3 levels: 
#'     fish meal, soybean meal, dried skim milk)}
#'   \item{percent}{Protein percentage in the diet (numeric with 4 values:
#'     9, 12, 15, and 18)}
#'   \item{conc}{Concentration of free plasma leucine, in mcg/ml}
#' }
#' @source Windels HF (1964) PhD thesis, Univ. of Minnesota. (Reported as
#'   Problem 10.8 in Oehlert G (2000) \emph{A First Course in Design and
#'   Analysis of Experiments}, licensed under Creative Commons,
#'   \url{http://users.stat.umn.edu/~gary/Book.html}.) Observations 7, 22, 23,
#'   31, 33, and 35 have been omitted, creating a more notable imbalance.
#' @examples
#'   pigs.lm <- lm(inverse(conc) ~ source + factor(percent), data = pigs)
#'   emmeans(pigs.lm, "source")
"pigs"

### ubds ###
#' Unbalanced dataset
#' 
#' This is a simulated unbalanced dataset with three factors
#' and two numeric variables. There are true relationships among these variables.
#' This dataset can be useful in testing or illustrating messy-data situations.
#' There are no missing data, and there is at least one observation for every 
#' factor combination; however, the \code{"cells"} attribute makes it simple
#' to construct subsets that have empty cells.
#' 
#' @format A data frame with 100 observations, 5 variables,
#'   and a special \code{"cells"} attribute:
#' \describe{
#'   \item{A}{Factor with levels 1, 2, and 3}
#'   \item{B}{Factor with levels 1, 2, and 3}
#'   \item{C}{Factor with levels 1, 2, and 3}
#'   \item{x}{A numeric variable}
#'   \item{y}{A numeric variable}
#' }
#' In addition, \code{attr(ubds, "cells")} consists of a named list of length 27 with the row numbers for
#' each combination of \code{A, B, C}. For example, 
#' \code{attr(ubds, "cells")[["213"]]} has the row numbers corresponding
#' to levels \code{A == 2, B == 1, C == 3}. The entries are ordered by
#' length, so the first entry is the cell with the lowest frequency.
#' @examples
#'  # Omit the three lowest-frequency cells
#'  low3 <- unlist(attr(ubds, "cells")[1:3]) 
#'  messy.lm <- lm(y ~ (x + A + B + C)^3, data = ubds, subset = -low3)
#'   

"ubds"