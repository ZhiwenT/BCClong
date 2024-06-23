test_that("BayesT works", {
  data("example")
  fit.BCC <- example
  set.seed(20220929)
  res <- BayesT(fit.BCC)

  expect_equal(round(res$T.obs, 3),
               c(3573.428,3619.330,3520.798,3570.080,3603.569))
  expect_equal(round(res$T.rep, 3),
               c(3798.714,3536.665,3588.700,3587.186,3878.272))
})

test_that("BayesT fails", {
  data("example")
  fit.BCC <- example
  set.seed(20220929)

  expect_error(BayesT("a"))
})
