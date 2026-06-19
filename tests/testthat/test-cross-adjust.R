context("cross-group P-value adjustment")

# cross.adjust is documented to be ignored unless all 'by' groups are the
# same size. Previously the guard only checked the NUMBER of by groups, so
# unequal sizes caused matrix() recycling and silently wrong adjusted P values.
test_that("cross.adjust is ignored for unequal-size by groups", {
    warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
    prs <- pairs(emmeans(warp.lm, ~ tension | wool))
    uneq <- prs[-1]   # wool A: 2 comparisons, wool B: 3 -> unequal sizes

    base <- suppressMessages(
        test(uneq, by = "wool", adjust = "tukey", cross.adjust = "none"))

    # must not warn (no recycling) ...
    expect_warning(
        cross <- suppressMessages(
            test(uneq, by = "wool", adjust = "tukey", cross.adjust = "bonferroni")),
        regexp = NA)
    # ... and must equal the 'none' result (cross.adjust ignored as documented)
    expect_equal(cross$p.value, base$p.value)
})

# Positive control: with EQUAL-size by groups, cross.adjust must still be
# applied (this is the documented, working case - the fix must not disable it).
test_that("cross.adjust is still applied for equal-size by groups", {
    warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)
    eq <- pairs(emmeans(warp.lm, ~ tension | wool))   # wool A: 3, wool B: 3

    base  <- suppressMessages(
        test(eq, by = "wool", adjust = "tukey", cross.adjust = "none"))
    cross <- suppressMessages(
        test(eq, by = "wool", adjust = "tukey", cross.adjust = "bonferroni"))

    # cross.adjust must change (inflate) the P values across the two by groups
    expect_false(isTRUE(all.equal(cross$p.value, base$p.value)))
    expect_true(all(cross$p.value >= base$p.value - 1e-12))
})
