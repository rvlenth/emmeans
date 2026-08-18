context("recover.data methods")

library(mgcv)

set.seed(42)
dat <- gamSim(1, n = 400, dist = "binary", scale = 0.3, verbose = FALSE)
dat$group <- factor(sample(c("Control", "TrtA", "TrtB"), size = 400, replace = TRUE))

test_that("data recovery for mgcv gams", {
    fit1 <- gam(y ~ group + s(x0) + s(x1) + s(x2), data = dat, family = binomial)
    expect_no_error(recover_data(fit1))
    fit2 <- gam(y ~ group + s(x0) + s(x1) + offset(x2), data = dat, family = binomial)
    expect_no_error(recover_data(fit2))
    #dat2 <- recover_data(fit2)
    fit3 <- glm(y ~ group + offset(x2), data = dat, family = binomial)
    #dat3 <- recover_data(fit3) 
    expect_no_error(recover_data(fit3))
})

# these are not great tests at this time, not even sure what they are testing, a bit too meta
#names(fit1$model); delete.response(terms(fit1))
#names(fit2$model); delete.response(terms(fit2))
#names(fit3$model); delete.response(terms(fit3))

