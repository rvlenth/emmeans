context("emmip")

# The lattice engine read its labels from a non-existent object (a typo:
# attr(emm, ...) instead of attr(emms, ...)), so axis/trace labels were
# silently dropped. They should match the labels used by the ggplot engine.
test_that("emmip lattice engine keeps axis labels", {
    skip_if_not_installed("lattice")
    warp.lm <- lm(breaks ~ wool * tension, data = warpbreaks)

    tmp <- tempfile(fileext = ".pdf"); pdf(tmp)
    on.exit({ dev.off(); unlink(tmp) })

    p_lat <- emmip(warp.lm, wool ~ tension, engine = "lattice")
    p_gg  <- emmip(warp.lm, wool ~ tension, engine = "ggplot")

    expect_s3_class(p_lat, "trellis")
    expect_false(is.null(p_lat$xlab))
    expect_false(is.null(p_lat$ylab))
    # lattice labels match the ggplot engine's labels
    expect_equal(p_lat$xlab, p_gg$labels$x)
    expect_equal(p_lat$ylab, p_gg$labels$y)
})
