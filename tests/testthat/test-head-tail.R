context("head/tail for emmGrid")

# head()/tail() for emmGrid should match base-R semantics, including a
# negative n. Previously the negative-n branches of head and tail were
# effectively swapped.
test_that("head()/tail() honor base-R negative-n semantics", {
    warp.rg <- ref_grid(lm(breaks ~ wool * tension, data = warpbreaks))
    ten <- as.character(warp.rg@grid$tension)   # 6 rows: L L M M H H
    expect_equal(length(ten), 6L)

    gten <- function(x) as.character(x@grid$tension)

    # positive n: first / last n rows
    expect_equal(gten(head(warp.rg, 3)), ten[1:3])
    expect_equal(gten(tail(warp.rg, 3)), ten[4:6])

    # negative n (base R): head(-2) drops the LAST 2; tail(-2) drops the FIRST 2
    expect_equal(gten(head(warp.rg, -2)), ten[1:4])
    expect_equal(gten(tail(warp.rg, -2)), ten[3:6])
})
