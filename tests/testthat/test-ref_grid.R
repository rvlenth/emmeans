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
})

colnames(ToothGrowth) <- c('len', 'choice of supplement', 'dose')
model <- stats::aov(`len` ~ `choice of supplement`, ToothGrowth)

test_that("Reference grid handles variables with spaces", {
    ref_grid(model, ~`choice of supplement`)
})
