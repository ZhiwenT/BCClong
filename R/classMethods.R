#' Generic print method for BCC objects
#'
#' @param x An object of class BCC.
#' @param ... further arguments passed to or from other methods.
#' @return Void function prints model information, no object return
#' @examples
#' # get data from the package
#' data(epil2)
#' fit.BCC <- epil2
#' print(fit.BCC)
#' @export
#' @useDynLib BCClong, .registration=TRUE

print.BCC <- function(x, ...) {
  cat("Total number of individual:\n")
  print(x$N[1])
  cat("\n")
  cat("Number of features:\n")
  print(x$R[1])
  cat("\n")
  cat("Mean adherence parameters by features:\n")
  print(colMeans(x$ALPHA))
  cat("\n")
  cat("Cluster proportions for global clusters:\n")
  print(colMeans(x$PPI))
  cat("\n")
  cat("Globle clusters table:\n")
  print(table(x$cluster.global))
  cat("\n")
  cat("Local clusters tables:\n")
  for (i in 1:length(x$cluster.local))
    print(table(x$cluster.local[[i]]))
  cat("\n")
  cat("Available components:\n")
  print(names(x))
}

#' Generic summary method for BCC objects
#'
#' @param object An object of class BCC.
#' @param ... further arguments passed to or from other methods.
#' @return Void function summarize model information, no object return
#' @examples
#' # get data from the package
#' data(epil2)
#' fit.BCC <- epil2
#' summary(fit.BCC)
#' @export
#' @method summary BCC
#' @useDynLib BCClong, .registration=TRUE

summary.BCC <- function(object, ...){
  cat("Total number of individual:\n")
  print(object$N[1])
  cat("\n")

  cat("Number of features:\n")
  print(object$R[1])
  cat("\n")

  cat("Cluster proportions statistics for global clusters:\n")
  print(summary(object$PPI))
  cat("\n")

  cat("Globle clusters table:\n")
  print(table(object$cluster.global))
  cat("\n")

  cat("Adherence parameters statistics by feature:\n")
  print(summary(object$ALPHA))
  cat("\n")

  cat("Local clusters statistics by feature:\n")
  for (i in 1:length(object$summary.stat$GA)){
    cat("Cluster statistics for feature", i, ":\n")
    print(object$summary.stat$GA[[i]])
  }
  cat("\n")

  cat("Variance-covariance matrix statistics for random effects by feature:\n")
  for (i in 1:length(object$summary.stat$SIGMA.SQ.U)){
    cat("Variance-covariance matrix statistics for feature", i, ":\n")
    print(object$summary.stat$SIGMA.SQ.U[[i]])
  }
  cat("\n")

  cat("Residual variance of continuous features statistics by feature:\n")
  for (i in 1:length(object$summary.stat$SIGMA.SQ.E)){
    cat("Residual variance statistics for feature", i, ":\n")
    print(object$summary.stat$SIGMA.SQ.E[[i]])
  }
  cat("\n")

  cat("Local clusters tables by feature:\n")
  for (i in 1:length(object$cluster.local)){
    cat("Clusters table for feature", i, ":\n")
    print(table(object$cluster.local[[i]]))
  }
}

#' Generic plot method for BCC objects
#'
#' @param x An object of class BCC.
#' @param ... further arguments passed to or from other methods.
#' @return Void function plot model object, no object return
#' @examples
#' # get data from the package
#' data(epil1)
#' fit.BCC <- epil1
#' plot(fit.BCC)
#' @export
#' @method plot BCC
#' @useDynLib BCClong, .registration=TRUE

plot.BCC <- function(x, ...){
  ncluster <- x$num.cluster
  nfeature <- x$R
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par))
  par(ask = TRUE)
  for (i in 1:nfeature){
    temp <- trajplot(fit=x,feature.ind=i,
             which.cluster = "local.cluster",
             title= bquote(paste("Local Clustering (",hat(alpha)[1] ==
                                   .(round(x$alpha[1],2)),")")),
             xlab="time (months)",ylab=paste("Feature", i))
    plot(temp)
    temp1 <- trajplot(fit=x,feature.ind=i,
             which.cluster = "global.cluster",
             title="Global Clustering",xlab="time (months)",
             ylab=paste("Feature", i))
    plot(temp1)
  }
  traceplot(fit=x, parameter="PPI",ylab="pi",xlab="MCMC samples")
  traceplot(fit=x, parameter="ALPHA",ylab="alpha",xlab="MCMC samples")
  for (i in 1:ncluster){
    for (j in 1:nfeature){
      traceplot(fit=x,cluster.indx = i, feature.indx=j,parameter="GA",
                ylab="GA",xlab="MCMC samples", title = paste("Feature", j))
    }
  }
}

# [END]
