test_that("trajplot works", {
  data("epil1")
  fit.BCC <- epil1
  # for local cluster
  p1 <- trajplot(fit=fit.BCC,feature.ind=1, which.cluster = "local.cluster",
           title= "Local Clustering",xlab="time (months)",
           ylab="anxiety",color=c("#00BA38", "#619CFF"))

  # for global cluster
  p2 <- trajplot(fit=fit.BCC,feature.ind=1,
          which.cluster = "global.cluster",
          title="Global Clustering",xlab="time (months)",
          ylab="anxiety",color=c("#00BA38", "#619CFF"))
  expect_equal(dim(p1$data), c(1789, 11))
  expect_equal(dim(p2$data), c(1789, 11))

  expect_type(p1$layers[[1]], "environment")
  expect_null(p1$scales$scales[[1]]$range$range)
  expect_identical(p1$labels$y, "anxiety")
  expect_type(p2$layers[[1]], "environment")
  expect_null(p2$scales$scales[[1]]$range$range)
  expect_identical(p2$labels$y, "anxiety")
})

test_that("errors",{
  data("epil1")
  fit.BCC <- epil1
  # for local cluster
  expect_error(trajplot(fit=fit.BCC,feature.ind=100, which.cluster = "local.cluster",
                 title= "Local Clustering",xlab="time (months)",
                 ylab="anxiety",color=c("#00BA38", "#619CFF")),
               "error in evaluating the argument")
  expect_error(trajplot(fit="fit.BCC",feature.ind=1, which.cluster = "local.cluster",
                        title= "Local Clustering",xlab="time (months)",
                        ylab="anxiety",color=c("#00BA38", "#619CFF")))
})

