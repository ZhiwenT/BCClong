#' Model Selection
#'
#' A function that calculates DIC and WAIC
#'
#'
#' @param fit an objective output from BCC.multi() function
#' @param fast_version ???
#'
#'
#' @return Returns the calculated score
#'
#'
#' @export
#' @import MASS
#' @import mclust
#' @import Rcpp
#' @importFrom LaplacesDemon WAIC
#' @useDynLib BCClong, .registration=TRUE

#sourceCpp("C:/Likelihood_general_final.cpp")

model.selection.criteria <- function(fit, fast_version=1){
  # calculate the log-likelihood
  log_lik <- LL(fit, fast_version = fast_version)
  Dev <- -2*colSums(log_lik)
  res <- WAIC(log_lik)
  WBIC <- -mean(colSums(log_lik))
  list(DIC = mean(Dev) + var(Dev)/2, WAIC = res$WAIC, WBIC=WBIC)
}
