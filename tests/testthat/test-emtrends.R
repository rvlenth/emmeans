context("Emtrends function")

pigs.lm = lm(log(conc) ~ source * poly(percent, 3), data = pigs)

test_that("emtrends works", {
    emt = emtrends(pigs.lm, ~ source, "percent")
    expect_equal(nrow(semt <- summary(emt)), 3)
    expect_equal(semt$percent.trend[1], .00429, tol = 0.0001)
    
    emtt = emtrends(pigs.lm, ~ source, "sqrt(percent)")
    expect_equal(summary(emtt)[["sqrt(percent).trend"]][1], .0309, tol = 0.001)
    
    emtp = emtrends(pigs.lm, ~ source, "percent", max.degree = 3)
    expect_equal(nrow(semtp <- summary(emtp)), 9)
    expect_equal(semtp$percent.trend[7], .001337, tol = 0.0001)
    
    emtpa = emtrends(pigs.lm, ~ source | percent, "percent", max.degree = 2,
                     at = list(percent = c(9,12,15,18)))
    expect_equal(nrow(summary(emtpa)), 24) # 3 sources * 2 degrees * 4 percents
})

