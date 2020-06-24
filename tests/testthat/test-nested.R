context("Nested structures")

set.seed(412.1948)
foo = data.frame(
    nest = factor(c(rep("A",8),rep("B",12),rep("C",4))),
    bird = factor(rep(1:6, each=4)),
    m = rep(c(2,3,0,-1,-2,5), each=4),
    b = rep(c(0,-1,1,2,2,3), each=4),
    x = 3 + 5*runif(24),
    e = 0.3*rnorm(24)
)
foo = transform(foo, resp = m*x + b + e)
foo1.lm = lm(resp ~ nest + bird*x, data = foo)
rg1 = ref_grid(foo1.lm)
foo2.lm = lm(resp ~ (nest + bird)*x, data = foo)
rg2 = ref_grid(foo2.lm)

test_that("nested EMM works", {
    emm1 = emmeans(rg1, "bird")
    expect_equal(nrow(summary(emm1)), 6)
    expect_equal(names(emm1@grid[1:2]), c("bird", "nest"))
    p1bn = predict(emmeans(emm1, "nest"))
    p1n = predict(emmeans(rg1, "nest"))
    p2n = predict(emmeans(rg2, "nest"))
    expect_equal(p1bn, p1n, tol = 1e-6)
    expect_equal(p1n, p2n, tol = 1e-6)
})

test_that("nested trends works", {
    emtb = emtrends(foo1.lm, "bird", "x")
    emtn = emtrends(foo1.lm, "nest", "x")
    emtbn = emmeans(emtb, "nest")
    expect_equal(nrow(summary(emtb)), 6)
    expect_equal(nrow(summary(emtn)), 3)
    expect_equal(predict(emtn), predict(emtbn), tol = 1e-6)
})

