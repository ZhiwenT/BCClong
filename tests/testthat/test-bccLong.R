test_that("BCC.multi works", {
  data("epil")
  dat <- epil
  set.seed(20220929)
  # example only, larger number of iteration required for accurate result
  fit.BCC <-  BCC.multi (
    mydat = list(dat$anxiety_scale,dat$depress_scale),
    dist = c("gaussian"),
    id = list(dat$id),
    time = list(dat$time),
    formula =list(y ~ time + (1|id)),
    num.cluster = 2,
    burn.in = 3,
    thin = 1,
    per =1,
    max.iter = 8)

  expect_type(fit.BCC, "list")
  expect_length(fit.BCC, 27)
  expect_length(fit.BCC$PPI, 10)
  expect_identical(round(fit.BCC$PPI[1,], 3), c(0.478, 0.522))
  expect_identical(round(fit.BCC$ALPHA[1,], 3), c(0.808, 0.959))
  expect_identical(dim(fit.BCC$GA[[1]]), as.integer(c(2,2,5)))
})


test_that("BCC.multi fail", {
  data("epil")
  dat <- epil
  set.seed(20220929)
  # example only, larger number of iteration required for accurate result
  expect_error(BCC.multi (
    mydat = list(dat$anxiety_scale,dat$depress_scale),
    dist = c("gaussian"),
    id = list(dat$id),
    time = list(dat$time),
    formula =list(y ~ time + (1|id)),
    num.cluster = 2,
    burn.in = 8,
    thin = 1,
    per =1,
    max.iter = 1), "invalid 'nrow' value")

  expect_error(BCC.multi (
    mydat = cbind(dat$anxiety_scale,dat$depress_scale),
    dist = c("gaussian"),
    id = list(dat$id),
    time = list(dat$time),
    formula =list(y ~ time + (1|id)),
    num.cluster = 2,
    burn.in = 2,
    thin = 1,
    per =1,
    max.iter = 5), "replacement has length zero")

  expect_error(BCC.multi (
    mydat = cbind(dat$anxiety_scale,dat$depress_scale),
    dist = c("gaussian"),
    id = list(dat$id),
    time = list(dat$time),
    num.cluster = 2,
    burn.in = 2,
    thin = 1,
    per =1,
    max.iter = 5), "argument \"formula\" is missing")

})
