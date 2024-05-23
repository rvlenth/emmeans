context("Contrast function")

pigs.lm = lm(log(conc) ~ source + factor(percent), data = pigs)
rg = ref_grid(pigs.lm)
rgg = add_grouping(rg, "group", "source", c("1", "2", "2"))
emms = emmeans(rg, "source")
emmg = emmeans(rgg, "group")
emmgs = emmeans(rgg, "source")

pigs.lmi = lm(log(conc) ~ source * factor(percent), data = pigs)
rgi = ref_grid(pigs.lmi)

test_that("Non-nested contrasts work", {
    expect_equal(nrow(summary(contrast(emms))), 3)
    expect_equal(nrow(summary(contrast(emms, "consec"))), 2)
    expect_equal(nrow(summary(contrast(rg, by = "source"))), 12)
    expect_equal(nrow(summary(pairs(rg, by = "source"))), 18)
    expect_equal(nrow(summary(pairs(rg))), 66)
})

test_that("Nested contrasts work", {
    expect_equal(nrow(summary(contrast(emmg, "consec"))), 1)
    expect_warning(summary(contrast(emmgs, "consec", by = "group"))) 
        # warning from cov2cor due to mvt adjustment
    expect_equal(nrow(summary(pairs(emmgs, by = "group"))), 2)
})

test_that("Interaction contrasts work", {
    expect_equal(nrow(summary(contrast(rgi, interaction = TRUE))), 12)
    expect_equal(nrow(summary(contrast(rgi, interaction = "consec"))), 6)
    expect_equal(nrow(summary(contrast(rgi, interaction = c("consec", "pairwise")))), 12)
    expect_equal(nrow(summary(contrast(rgi, interaction = c("pairwise","consec")))), 9)
})


