context("Reference grids")

pigs.lm = lm(log(conc) ~ source + factor(percent), data = pigs)
rg = ref_grid(pigs.lm)
rg1 = ref_grid(pigs.lm, at = list(source = "soy", percent = 12))

pigs.lm2 = update(pigs.lm, conc ~ source + percent, data = pigs)
rg2 = ref_grid(pigs.lm2)
rg2a = ref_grid(pigs.lm2, at = list(source = c("fish", "soy"), percent = 10))
rg2c = ref_grid(pigs.lm2, cov.reduce = FALSE)
rg2m = ref_grid(pigs.lm2, cov.reduce = min)

pigs.lm3 = update(pigs.lm2, . ~ source + source:factor(percent))
pigs = transform(pigs, sp = interaction(source, percent))
pigs.lm4 = update(pigs.lm2, . ~ source + sp)

test_that("Reference grid is constructed correctly", {
    expect_equal(nrow(rg@grid), 12)
    expect_equal(nrow(rg1@grid), 1)
    expect_equal(nrow(rg2@grid), 3)
    expect_equal(nrow(rg2a@grid), 2)
    expect_equal(nrow(rg2c@grid), 12)
    expect_equal(nrow(rg1@grid), 1)
    expect_equal(length(rg@levels), 2)
    expect_equal(rg2@levels$percent, mean(pigs$percent))
    expect_equal(rg2m@levels$percent, min(pigs$percent))
})

test_that("Reference grid extras are detected", {
    expect_equal(rg@misc$tran, "log")
    expect_true(is.null(rg2@misc$tran))
    expect_true(is.null(rg2@model.info$nesting))
    expect_is(ref_grid(pigs.lm3)@model.info$nesting, "list") # see note above
    expect_is(ref_grid(pigs.lm4)@model.info$nesting, "list") # see note above
})

colnames(ToothGrowth) <- c('len', 'choice of supplement', 'dose')
model <- stats::aov(`len` ~ `choice of supplement`, ToothGrowth)

test_that("Reference grid handles variables with spaces", {
    expect_output(str(ref_grid(model, ~`choice of supplement`)), "choice of supplement")
})

# models outside of data.frames
x = 1:10
y = rnorm(10)
mod1 = with(pigs, lm(log(conc) ~ source + factor(percent)))
test_that("ref_grid works with no data or subset", {
    expect_silent(ref_grid(lm(y ~ x)))
    expect_silent(ref_grid(mod1))
})

# Multivariate models
MOats.lm <- lm (yield ~ Block + Variety, data = MOats)
MOats.rg <- ref_grid (MOats.lm, mult.levs = list(
    trt = LETTERS[1:2], dose = as.character(1:2))
)
test_that("We can construct multivariate reference grid", {
    expect_equal(nrow(MOats.rg@grid), 72)
    expect_equal(length(MOats.rg@levels), 4)
})
MOats.nrg <- ref_grid (MOats.lm, mult.levs = list(nitro = c(0,.2,.4,.6)),
                       at = list(nitro = c(0.0001, .39998, .5)))
test_that("Fuzzy matching of numerical mult.levels works", {
    expect_equal(length(MOats.nrg@levels$nitro), 2)
    expect_equal(MOats.nrg@levels$nitro, c(0, .4), 0.001)
})

### Nuisance factors
MOats.rgn = ref_grid(MOats.lm, nuisance = "Block")
MOats.emm = emmeans(MOats.lm, ~ Variety * rep.meas)
test_that("We get same predictions with and without nuisance specs", {
    expect_equal(predict(MOats.rgn), predict(MOats.emm), 0.001)
})

# Missing and NA levels
miss.df = data.frame(x = factor(c("a", "a", "b", NA), levels = c("a", "b", "c", NA)), y = 1:4)
miss.lm = lm(y ~ x, data = miss.df)
miss.rg1 = ref_grid(miss.lm)
miss.rg2 = ref_grid(miss.lm, data = miss.df)
# Now try allowing NA levels
miss.dfa = transform(miss.df, x = factor(x, exclude = NULL))
miss.lma = lm(y ~ x, data = miss.dfa)
miss.rg1a = ref_grid(miss.lma)
miss.rg2a = ref_grid(miss.lma, data = miss.dfa)
test_that("Reference grid handles missing values", {
    expect_equal(length(miss.rg1@levels$x), 2)
    expect_equal(length(miss.rg2@levels$x), 2)
    expect_equal(length(miss.rg1a@levels$x), 3)
    expect_equal(length(miss.rg2a@levels$x), 3)
    expect_equal(inherits(with_emm_options(allow.na.levs = FALSE, ref_grid(miss.lma)), 
                 "try-error"), TRUE)
})

