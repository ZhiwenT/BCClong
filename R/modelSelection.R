#' Model selection
#'
#' A function that calculates DIC and WAIC for model selection
#'
#' @param fit an objective output from BCC.multi() function
#' @param fast_version if fast_verion=TRUE (default), then compute the DIC and WAIC using
#'    the first 100 MCMC samples (after burn-in and thinning) . If fast_version=FALSE, then
#'    compute the DIC and WAIC using all MCMC samples (after burn-in and thinning)
#' @return Returns the calculated score
#' @examples
#' #import data
#' data(example1)
#' fit.BCC <- example1
#' res <- model.selection.criteria(fit.BCC, fast_version=TRUE)
#' res
#'
#' @export
#' @import MASS
#' @import mclust
#' @import Rcpp
#' @importFrom LaplacesDemon WAIC
#' @useDynLib BCClong, .registration=TRUE

model.selection.criteria <- function(fit, fast_version=TRUE){
  if(!is.logical(fast_version)){
    stop("fast_version should be either TRUE or FALSE")
  }
  # calculate the log-likelihood
  log_lik <- LL(fit, fast_version = as.numeric(fast_version))
  Dev <- -2*colSums(log_lik)
  res <- WAIC(log_lik)
  WBIC <- -mean(colSums(log_lik))
  list(DIC = mean(Dev) + var(Dev)/2, WAIC = res$WAIC, WBIC=WBIC)
}

# [END]
