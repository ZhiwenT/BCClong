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
#' data(epil1)
#' fit.BCC <- epil1
#' traceplot(fit=fit.BCC, parameter="PPI",ylab="pi",xlab="MCMC samples")
#'
#' @export
#' @import ggplot2
#' @import gridExtra
#' @useDynLib BCClong, .registration=TRUE

traceplot <- function(fit, cluster.indx=1, feature.indx=1, parameter="PPI",
                      xlab = NULL, ylab = NULL, ylim = NULL, xlim = NULL, title = NULL) {
  num.cluster <- fit$num.cluster
  num.sample <- (fit$max.iter - fit$burn.in)/fit$thin
  R <- fit$R
  x <- 1:num.sample

  if (!parameter %in% c("PPI", "ALPHA", "GA", "SIGMA.SQ.U", "SIGMA.SQ.E")){
    stop("invalid parameter")
  }

  plot_list <- list()

  if (parameter == "PPI"){
    for (j in 1:num.cluster) {
      y <- fit$PPI[,j]
      df <- data.frame(x = x, y = y)
      p <- ggplot(df, aes(x = x, y = y)) +
        geom_line() +
        labs(title = paste0("Cluster ", j), x = xlab, y = ylab) +
        theme_minimal()
      if (!is.null(ylim)) p <- p + ylim(ylim)
      if (!is.null(xlim)) p <- p + xlim(xlim)
      plot_list[[j]] <- p
    }
  }

  if (parameter == "ALPHA") {
    for (j in 1:R) {
      y <- fit$ALPHA[,j]
      df <- data.frame(x = x, y = y)
      p <- ggplot(df, aes(x = x, y = y)) +
        geom_line() +
        labs(title = paste0("Feature ", j), x = xlab, y = ylab) +
        theme_minimal()
      if (!is.null(ylim)) p <- p + ylim(ylim)
      if (!is.null(xlim)) p <- p + xlim(xlim)
      plot_list[[j]] <- p
    }
  }

  if (parameter == "GA") {
    dim.GA <- dim(fit$GA[[feature.indx]][cluster.indx,,])[1]
    for (j in 1:dim.GA) {
      y <- fit$GA[[feature.indx]][cluster.indx, j, ]
      df <- data.frame(x = x, y = y)
      p <- ggplot(df, aes(x = x, y = y)) +
        geom_line() +
        labs(title = title, x = xlab, y = ylab) +
        theme_minimal()
      if (!is.null(ylim)) p <- p + ylim(ylim)
      if (!is.null(xlim)) p <- p + xlim(xlim)
      plot_list[[j]] <- p
    }
  }

  if (parameter == "SIGMA.SQ.U") {
    for (j in 1:num.cluster) {
      y <- fit$SIGMA.SQ.U[[feature.indx]][cluster.indx, j, ]
      df <- data.frame(x = x, y = y)
      p <- ggplot(df, aes(x = x, y = y)) +
        geom_line() +
        labs(title = title, x = xlab, y = ylab) +
        theme_minimal()
      if (!is.null(ylim)) p <- p + ylim(ylim)
      if (!is.null(xlim)) p <- p + xlim(xlim)
      plot_list[[j]] <- p
    }
  }

  if (parameter == "SIGMA.SQ.E") {
    if (fit$dist[feature.indx] == "gaussian") {
      for (j in 1:num.cluster) {
        y <- fit$SIGMA.SQ.E[[feature.indx]][, cluster.indx]
        df <- data.frame(x = x, y = y)
        p <- ggplot(df, aes(x = x, y = y)) +
          geom_line() +
          labs(title = title, x = xlab, y = ylab) +
          theme_minimal()
        if (!is.null(ylim)) p <- p + ylim(ylim)
        if (!is.null(xlim)) p <- p + xlim(xlim)
        plot_list[[j]] <- p
      }
    } else {
      message("SIGMA.SQ.E is not estimated for features with Binomial or Poisson distribution")
    }
  }

  do.call(grid.arrange, c(plot_list, nrow = 1))
}

# [END]
