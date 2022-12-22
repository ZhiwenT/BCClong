#' Model Selection
#'
#' A function that calculates of DIC and WAIC
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
  LL <- LL(fit, fast_version = fast_version)
  Dev <- -2*colSums(LL)
  res <- WAIC(LL)
  WBIC <- -mean(colSums(LL))
  result <- list(DIC = mean(Dev) + var(Dev)/2, WAIC = res$WAIC, WBIC=WBIC)
  return(result)
}
