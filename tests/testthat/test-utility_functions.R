context("utility functionality")

m1 <- lm(breaks ~ wool*tension, data = warpbreaks)
em1 <- emmeans(m1, ~ wool | tension)

emm_options(pval.digits = 8); invisible(capture.output(joint_tests(em1)))
pdig_1 = get_emm_option("pval.digits") 

emm_options(pval.digits = "eight"); invisible(capture.output(joint_tests(em1)))
pdig_2 = get_emm_option("pval.digits") 

emm_options(pval.digits = 1); invisible(capture.output(joint_tests(em1)))
pdig_3 = get_emm_option("pval.digits") 

emm_options(pval.digits = c(3, 6));  invisible(capture.output(joint_tests(em1)))
pdig_4 = get_emm_option("pval.digits")[1] 

test_that("digits for p-values process", {
  expect_equal(pdig_1, 6)
  expect_equal(pdig_2, 4)
  expect_equal(pdig_3, 2)
  expect_equal(pdig_4, 3)
})