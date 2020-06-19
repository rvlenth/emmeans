context("Estimated marginal means")

pigs.lm = lm(log(conc) ~ source + factor(percent), data = pigs)
rg = ref_grid(pigs.lm)
rgg = add_grouping(rg, "group", "source", c("1", "2", "2"))

test_that("Character interface works", {
    expect_equal(confint(emmeans(rg, "source"))$emmean, c(3.39, 3.67, 3.80), tol = 0.01)
    expect_equal(nrow(emmeans(rg, c("source", "percent"))@grid), 12)
    expect_equal(nrow(emmeans(rg, "source", by = "percent")@grid), 12)
    expect_equal(nrow(emmeans(rg, c("1"))@grid), 1)
})

test_that("Formula interface works", {
    expect_equal(confint(emmeans(rg, ~ source))$emmean, c(3.39, 3.67, 3.80), tol = 0.01)
    expect_equal(nrow(emmeans(rg, ~ source * percent)@grid), 12)
    expect_equal(nrow(emmeans(rg, ~ source | percent)@grid), 12)
    expect_equal(nrow(emmeans(rg, ~ 1)@grid), 1)
    expect_equal(nrow(emmeans(rg, ~ 1 | percent)@grid), 4)
    expect_equal(nrow(emmeans(rg, ~ percent | 1)@grid), 4)
})

# nesting
test_that("Nested EMMs work", {
    expect_equal(nrow(emmeans(rgg, ~ group)@grid), 2)
    expect_equal(nrow(emmeans(rgg, ~ source)@grid), 6) #   3 rows have 0 weight
    expect_equal(nrow(confint(emmeans(rgg, ~ source))), 3)
    expect_equal(colnames(emmeans(rgg, ~ source)@grid)[1:2], c("source","group"))
})

