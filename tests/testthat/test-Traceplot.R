test_that("traceplot works", {
  filePath <- system.file("extdata", "epil1.rds", package = "BCClong")
  fit.BCC <- readRDS(filePath)
  temp <- traceplot(fit=fit.BCC, parameter="PPI",
            ylab="pi",xlab="MCMC samples")
  expect_null(temp)
})

test_that("invailid parameter", {
  filePath <- system.file("extdata", "epil1.rds", package = "BCClong")
  fit.BCC <- readRDS(filePath)
  expect_error(traceplot(fit=fit.BCC, parameter="POI",
                    ylab="pi",xlab="MCMC samples"))

})
