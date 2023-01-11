#' Model selection
#'
#' A function that calculates DIC and WAIC for model selection
#'
#' @param fit an objective output from BCC.multi() function
#' @param fast_version if fast_verion=1 (default), then compute the DIC and WAIC using
#'    the first 100 MCMC samples (after burn-in and thinning) . If fast_version=0, then
#'    compute the DIC and WAIC using all MCMC samples (after burn-in and thinning)
#' @return Returns the calculated score
#' @examples
#' #import data
#' filePath <- system.file("extdata", "example1.rds", package = "BCClong")
#' fit.BCC <- readRDS(filePath)
#' res <- model.selection.criteria(fit.BCC, fast_version=1)
#' res
#'
#' @export
#' @import MASS
#' @import mclust
#' @import Rcpp
#' @importFrom LaplacesDemon WAIC
#' @useDynLib BCClong, .registration=TRUE

model.selection.criteria <- function(fit, fast_version=1){
  # calculate the log-likelihood
  log_lik <- LL(fit, fast_version = fast_version)
  Dev <- -2*colSums(log_lik)
  res <- WAIC(log_lik)
  WBIC <- -mean(colSums(log_lik))
  list(DIC = mean(Dev) + var(Dev)/2, WAIC = res$WAIC, WBIC=WBIC)
}

# [END]
