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

data(penguins)
penguins$bill_len[[1]] = NA


test_that("data recovery for missing/subset data", {
  peng_lm = lm(bill_dep ~ log(bill_len) * species, data = subset(penguins, body_mass > 2500))
  expect_no_error(recover_data(peng_lm))
})

test_that("data recovery for variables with spaces", {
  colnames(ToothGrowth) <- c('len', 'choice of supplement', 'dose')
  tooth = stats::aov(`len` ~ `choice of supplement`, ToothGrowth)
  expect_no_error(recover_data(tooth))
})

test_that("data recovery for interactions", {
  org.quad <- lm(cbind(sales1, sales2) ~ poly(price1, price2, degree = 2) + day + store, data = oranges)
  org.int <- lm(cbind(sales1, sales2) ~ price1 * price2 + day + store, data = oranges)
  org.add <- lm(cbind(sales1, sales2) ~ price1 + price2 + day + store, data = oranges)
  expect_no_error(recover_data(org.quad))
  expect_no_error(recover_data(org.int))
  expect_no_error(recover_data(org.add))
})
