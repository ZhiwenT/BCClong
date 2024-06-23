test_that("model.selection.criteria works", {
  data("example1")
  fit.BCC <- example1
  res <- model.selection.criteria(fit.BCC, fast_version=T)

  expect_equal(res$DIC, 3026.2626)
  expect_equal(res$WAIC, 2884.56245)
  expect_equal(res$WBIC, 1436.34092)
})

test_that("model.selection.criteria fail", {
  data("example1")
  fit.BCC <- example1

  expect_error(model.selection.criteria(fit.BCC, fast_version=2),
               "fast_version should be either TRUE or FALSE")
  expect_error(model.selection.criteria("a", fast_version=T),
               "Object was created without names")
})
