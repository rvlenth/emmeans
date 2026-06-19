context("response transformations")

# atanh transformation: mu.eta derivative must work so type = "response"
# can compute back-transformed estimates and delta-method SEs.
# Regression test for the typo `tanh^2(eta)` in the atanh mu.eta.
test_that("atanh response transformation supports type = 'response'", {
    rng <- range(iris$Sepal.Width)
    r <- 1.8 * (iris$Sepal.Width - rng[1]) / (rng[2] - rng[1]) - 0.9  # in (-0.9, 0.9)
    dat <- data.frame(r = r, Species = iris$Species)
    mod <- lm(atanh(r) ~ Species, data = dat)

    lk <- summary(emmeans(mod, "Species"))                  # link (atanh) scale
    rs <- summary(emmeans(mod, "Species"), type = "response")  # must not error

    # back-transform is tanh() of the link-scale means
    expect_equal(rs$response, tanh(lk$emmean), tolerance = 1e-8)
    # delta-method SEs are finite and positive
    expect_true(all(is.finite(rs$SE)) && all(rs$SE > 0))
})
