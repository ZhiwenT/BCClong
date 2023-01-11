test_that("model.selection.criteria works", {
  filePath <- system.file("extdata", "example1.rds", package = "BCClong")
  fit.BCC <- readRDS(filePath)
  res <- model.selection.criteria(fit.BCC, fast_version=1)

  expect_equal(res$DIC, 3026.2626)
  expect_equal(res$WAIC, 2884.56245)
  expect_equal(res$WBIC, 1436.34092)
})

test_that("model.selection.criteria fail", {
  filePath <- system.file("extdata", "example1.rds", package = "BCClong")
  fit.BCC <- readRDS(filePath)

  expect_error(model.selection.criteria(fit.BCC, fast_version=2),
               "fast_version should be either 0 or 1")
  expect_error(model.selection.criteria("a", fast_version=1),
               "Object was created without names")
})
