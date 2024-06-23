#' Goodness of fit.
#'
#' This function assess the model goodness of fit by calculate the
#' discrepancy measure T(bm(y), bm(Theta)) with following steps
#' (a) Generate T.obs based on the MCMC samples
#' (b) Generate T.rep based on the posterior distribution of the parameters
#' (c) Compare  T.obs and T.rep, and calculate the P values.
#'
#' @param fit an objective output from BCC.multi() function
#' @return Returns a dataframe with length equals to 2 that contains
#'         observed and predict value
#' @examples
#' #import data
#' data(example)
#' fit.BCC <- example
#' BayesT(fit.BCC)
#'
#' @export
#' @importFrom stats rnorm rpois rbinom
#' @importFrom mvtnorm rmvnorm
#' @useDynLib BCClong, .registration=TRUE

BayesT <- function(fit){

  max.iter <- fit$max.iter
  burn.in <- fit$burn.in
  thin <- fit$thin

  N <- fit$N
  R <- fit$R
  K <- fit$K
  k <- fit$k
  dist <- fit$dist
  alpha <- fit$alpha
  num.cluster <- fit$num.cluster
  xt <- fit$dat[[1]]$time
  dat <- fit$dat
  summary.stat <- fit$summary.stat
  ZZ.LOCAL <- fit$ZZ.LOCAL
  THETA <- fit$THETA
  num.sample <-  (max.iter - burn.in)/thin

  T.obs <- NULL
  T.rep <- NULL
  for (count in 1:num.sample){ # at each iteration

    #--------------------------------------------------------------#
    # (a) Generate T.obs based on the MCMC samples
    #--------------------------------------------------------------#
    # count <- 1
    #--------------------------------------------------------#
    ga <- vector(mode = "list", length = R)
    sigma.sq.u 	<- vector(mode = "list", length = R)
    sigma.sq.e  <- vector(mode = "list", length = R)
    zz.local 	<- vector(mode = "list", length = R)
    theta 	<- vector(mode = "list", length = R)

    T.tmp <- 0
    for (s in 1:R){
      zz.local[[s]] <- fit$ZZ.LOCAL[[s]][count,]
      ga[[s]] <- fit$GA[[s]][,,count]
      sigma.sq.e[[s]] <- fit$SIGMA.SQ.E[[s]][count,]
      sigma.sq.u[[s]] <- fit$SIGMA.SQ.U[[s]][,,count]
      theta[[s]] <- fit$THETA[[s]][,,count]
      y <- dat[[s]]$y
      xt <- dat[[s]]$time
      ids <- dat[[s]]$id
      for (i in 1:N){
        m <- matrix(cbind(1,xt[which(ids==i)],I(xt[which(ids==i)]^2),
                          I(xt[which(ids==i)]^3))[,1:k[[s]]],ncol=k[[s]])
        mz <- matrix(cbind(1,xt[which(ids==i)],I(xt[which(ids==i)]^2),
                           I(xt[which(ids==i)]^3))[,1:K[[s]]],ncol=K[[s]])
        for (j in 1:num.cluster){
          g <- matrix(ga[[s]][j,],ncol=k[[s]]) %*% t(m)  +
            matrix(theta[[s]],nrow=N)[i,] %*% t(mz)

          if (dist[[s]] == "gaussian"){
            T.tmp <- T.tmp + (zz.local[[s]][i]==j)*
              sum(((y[which(ids==i)] - g)^2/sigma.sq.e[[s]][j]))
          }
          if (dist[[s]] == "poisson"){
            mu <- exp(g)
            var <- exp(g)
            T.tmp <- T.tmp + (zz.local[[s]][i]==j)*
              sum(((y[which(ids==i)] - mu)^2/var))
          }
          if (dist[[s]] == "binomial"){
            gt <- exp(g)
            mu <- gt/(1+gt)
            var <- mu*(1-mu)
            T.tmp <- T.tmp + (zz.local[[s]][i]==j)*
              sum(((y[which(ids==i)] - mu)^2/var))
          }
        }
      }
    }
    T.obs <- c(T.obs,T.tmp)

    #--------------------------------------------------------------#
    # (b) Generate T.rep based on the MCMC samples
    #--------------------------------------------------------------#
    T.tmpp <- NULL
    for (s in 1:R){
      zz.local[[s]] <- fit$ZZ.LOCAL[[s]][count,]
      ga[[s]] <- fit$GA[[s]][,,count]
      sigma.sq.e[[s]] <- fit$SIGMA.SQ.E[[s]][count,]
      sigma.sq.u[[s]] <- fit$SIGMA.SQ.U[[s]][,,count]
      theta[[s]] <- fit$THETA[[s]][,,count]
      y <- dat[[s]]$y
      xt <- dat[[s]]$time
      ids <- dat[[s]]$id
      nobss <- table(ids)
      for (i in 1:N){
        m <- matrix(cbind(1,xt[which(ids==i)],I(xt[which(ids==i)]^2),
                          I(xt[which(ids==i)]^3))[,1:k[[s]]],ncol=k[[s]])
        mz <- matrix(cbind(1,xt[which(ids==i)],I(xt[which(ids==i)]^2),
                           I(xt[which(ids==i)]^3))[,1:K[[s]]],ncol=K[[s]])
        for (j in 1:num.cluster){
          # generate the random effect from the multivariate normal distribution;
          if (K[[s]]==1){
            theta.new <-  rnorm(1,mean= 0,
                                sd= sqrt(summary.stat$SIGMA.SQ.U[[s]][1,,j]))
          }
          if (K[[s]] > 1){
            theta.new <-  mvtnorm::rmvnorm(1,
                          mean= rep(0,K[[s]]),
                          sigma= diag(summary.stat$SIGMA.SQ.U[[s]][1,,j]))
          }
          g <- matrix(summary.stat$GA[[s]][1,,][j,],ncol=k[[s]])  %*% t(m) +
            theta.new %*% t(mz)
          if (dist[[s]] == "gaussian"){
            if (K[[s]]==1){
              y.rep <- rnorm(1,mean=g,
                             sd= sqrt(summary.stat$SIGMA.SQ.E[[s]][1,j]))
            }
            if (K[[s]] > 1){
              y.rep <-  mvtnorm::rmvnorm(1,mean= g ,
                        sigma= summary.stat$SIGMA.SQ.E[[s]][1,j]*
                          diag(1,nobss[i]))
            }

            T.tmpp <- c(T.tmpp,(zz.local[[s]][i]==j)*sum(((y.rep - g)^2/
                                 summary.stat$SIGMA.SQ.E[[s]][1,j])))
          }
          if (dist[[s]] == "poisson"){
            mu <- exp(g)
            var <- exp(g)
            y.rep <- rpois(length(mu),lambda=mu)
            T.tmpp <- c(T.tmpp,
                        (zz.local[[s]][i]==j)*sum(((y.rep - mu)^2/var)))
          }
          if (dist[[s]] == "binomial"){
            gt <- exp(g)
            mu <- gt/(1+gt)
            var <- mu*(1-mu)
            y.rep <- rbinom(length(mu),size=1,prob=mu)
            T.tmpp <- c(T.tmpp,(zz.local[[s]][i]==j)*sum(((y.rep - mu)^2/var)))
          }
        }
      }
    }
    T.rep <- c(T.rep,sum(T.tmpp))
  }
  #result <- list(T.obs = T.obs, T.rep = T.rep)
  result <- data.frame(T.obs = T.obs, T.rep = T.rep)
  return(result)
}

# [END]
