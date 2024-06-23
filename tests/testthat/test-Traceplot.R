test_that("traceplot works", {
  data("epil1")
  fit.BCC <- epil1
  temp <- traceplot(fit=fit.BCC, parameter="PPI",
            ylab="pi",xlab="MCMC samples")
  expect_length(temp, 2)
})

test_that("invailid parameter", {
  data("epil1")
  fit.BCC <- epil1
  expect_error(traceplot(fit=fit.BCC, parameter="POI",
                    ylab="pi",xlab="MCMC samples"))

})
