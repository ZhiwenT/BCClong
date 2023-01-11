#' Trace plot function
#'
#' To visualize the MCMC chain for model parameters
#'
#' @param fit an objective output from BCC.multi() function.
#' @param cluster.indx a numeric value. For cluster-specific parameters,
#'                    specifying cluster.indx will generate the trace plot for
#'                    the corresponding cluster.
#' @param feature.indx a numeric value. For cluster-specific parameters,
#'                     specifying feature.indx will generate the trace
#'                     plot for the corresponding cluster.
#' @param parameter a character value. Specify which parameter for which the
#'                  trace plot will be generated. The value can be "PPI" for pi,
#'                  alpha for alpha, "GA" for gamma, "SIGMA.SQ.U" for Sigma
#'                  and "SIGMA.SQ.E" for sigma.
#' @param xlab Label for x axis
#' @param ylab Label for y axis
#' @param ylim The range for y axis
#' @param xlim The range for x axis
#' @param title Title for the trace plot
#' @return void function with no return value, only show plots
#' @examples
#' # get data from the package
#' filePath <- system.file("extdata", "epil1.rds", package = "BCClong")
#' fit.BCC <- readRDS(filePath)
#' traceplot(fit=fit.BCC, parameter="PPI",ylab="pi",xlab="MCMC samples")
#'
#' @export
#' @importFrom graphics plot
#' @useDynLib BCClong, .registration=TRUE

traceplot <- function(fit, cluster.indx=1, feature.indx=1, parameter="PPI",  xlab = NULL,
                      ylab = NULL,
                      ylim = NULL,
                      xlim = NULL,
                      title = NULL) {
  num.cluster <- fit$num.cluster
  num.sample <- (fit$max.iter - fit$burn.in)/fit$thin
  R <- fit$R
  x <- 1:num.sample

  if (!parameter %in% c("PPI", "ALPHA", "GA", "SIGMA.SQ.U", "SIGMA.SQ.E")){
    stop("invalid parameter")
  }

  if (parameter=="PPI"){
    opar <- par(mfrow=c(1,num.cluster))
    on.exit(par(opar))
    for (j in 1:  num.cluster){
      y <- fit$PPI[,j]
      plot(x,y,type="l",xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,
           main=paste0("Cluster ",j), lwd=1.5)}}

  if (parameter=="ALPHA") {
    opar <- par(mfrow=c(1,R))
    on.exit(par(opar))
    for (j in 1: R){
      y <- fit$ALPHA[,j]
      plot(x,y,type="l",xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,
           main=paste0("Feature ",j),  lwd=1.5)
    }
  }
  if (parameter=="GA") {
    dim.GA <- dim(fit$GA[[feature.indx]][cluster.indx,,])[1]
    opar <- par(mfrow=c(1,dim.GA))
    on.exit(par(opar))
    for (j in 1:dim.GA){
      y <- fit$GA[[feature.indx]][cluster.indx,j,]
      plot(x,y,type="l",xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,
           main=title, lwd=1.5)
    }
  }
  if (parameter=="SIGMA.SQ.U") {
    opar <- par(mfrow=c(1,num.cluster))
    on.exit(par(opar))
    for (j in 1:  num.cluster){
      y <- fit$SIGMA.SQ.U[[feature.indx]][cluster.indx,j,]
      plot(x,y,type="l",xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,
           main=title, lwd=1.5)
    }
  }
  if (parameter=="SIGMA.SQ.E") {
    if (fit$dist[feature.indx]=="gaussian"){
      opar <- par(mfrow=c(1,num.cluster))
      on.exit(par(opar))
      for (j in 1:  num.cluster){
        y <- fit$SIGMA.SQ.E[[feature.indx]][,cluster.indx]
        plot(x,y,type="l",xlab=xlab,ylab=ylab,xlim=xlim,ylim=ylim,
             main=title, lwd=1.5)
      }
    }
    else{
          message("SIGMA.SQ.E is not estimated for features with Binomial
                or Poisson distribution")
    }
  }
}

# [END]
