#' Compute a Bayesian Consensus Clustering model for mixed-type longitudinal data
#'
#' A function that performs clustering on mixed-type (continuous, discrete and
#' categorical) longitudinal markers using Bayesian consensus clustering method
#' with MCMC sampling
#'
#'
#'
#' @param mydat List of R outcomes (R is the number of outcome, R>=2)
#' @param id id-variable: starting from 1 to N
#' @param time time variable
#' @param num.cluster number of cluster
#' @param formula fixed and random effect
#' @param dist "gaussian","poisson","binomial", distribution of the outcome
#' @param alpha.common 1 - common alpha, 0 - separate alphas for each outcome
#' @param initials List of initials for: zz, zz.local ga, sigma.sq.u, sigma.sq.e,
#'                 Default is NULL
#' @param sig.var 1 - unstructure random effect variance,
#'                0 - diagonal random effect variance structure,
#'                default is 0
#' @param sigma.sq.e.common 1 - estimate common residual variance across all groups,
#'                          0 - estimate distinct residual variance, default is 1
#' @param hyper.par Hyper parameter contains a list of variables includes
#'                  delta  = 1, a.star = 1, b.star = 1, aa0  = 1e-3, bb0 = 1e-3,
#'                  cc0 = 1e-3, ww0 = 0, vv0 = 1e3, dd0 = 1e-3, rr0 = 4, RR0 = 3
#' @param c.ga.tuning tuning parameter for MH algorithm (fixed effect parameters),
#'                    each parameter corresponds to an outcome/marker, default
#'                    value equals list(1,1,0.5)
#' @param c.theta.tuning tuning parameter for MH algorithm (random effect),
#'                       each parameter corresponds to an outcome/marker,
#'                       default value equals list(1,1,0.5)
#' @param adaptive.tuning adaptive tuning parameters, 1 - yes, 0 - no,
#'                        default is 1
#' @param align.clusters  assign clusters, default is 1
#' @param tuning.freq     tuning frequency, default is 20
#' @param initial.cluster.membership "mixAK" or "random" or "PAM" or "input" -
#'                                  input initial cluster membership for local
#'                                  clustering, default is "PAM"
#' @param input.initial.cluster.membership if use "input",
#'                                  option input.initial.cluster.membership
#'                                  must not be empty, default is NULL
#' @param initial.global.cluster.membership input initial cluster
#'                                  membership for global clustering
#'                                  default is NULL
#' @param seed.initial seed for initial clustering
#'                    (for initial.cluster.membership = "mixAK")
#'                    default is 2080
#' @param print.info print model information at each iteration, default is true
#' @param burn.in number of samples discarded
#' @param thin thinning
#' @param per output information every "per" iteration
#' @param max.iter maximum number of iteration
#'
#'
#' @return Returns a model contains clustering information
#'
#'
#' @export
#' @import label.switching
#' @import lme4
#' @import mclust
#' @import MCMCpack
#' @import mixAK
#' @importFrom mvtnorm rmvnorm
#' @import nnet
#' @import Rcpp
#' @import Rmpfr
#' @import truncdist
#' @import cluster
#' @importFrom coda geweke.diag
#' @importFrom stats binomial poisson sd var

BCC.multi <- function(
    mydat,                 # List of R outcomes (R is the number of outcome, R>=2)
    id,                    # id-variable: starting from 1 to N
    time,                  # time variable
    num.cluster,           # number of cluster
    formula,               # fixed and random effect
    dist,                  # "gaussian","poisson","binomial", distribution of the outcome
    alpha.common,          # 1 - common alpha, 0 - separate alphas for each outcome
    initials = NULL,       # List of initials for: zz, zz.local ga, sigma.sq.u, sigma.sq.e,
    sig.var  = 0,          # 1 - unstructure random effect variance, 0 - diagonal random effect variance structure
    sigma.sq.e.common = 1, # 1 - estimate common residual variance across all groups, 0 - estimate distinct residual variance
    hyper.par = list(
      delta  = 1,
      a.star = 1,
      b.star = 1,
      aa0    = 1e-3,
      bb0    = 1e-3,
      cc0    = 1e-3,
      ww0    = 0,
      vv0    = 1e3,
      dd0    = 1e-3,
      rr0    = 4,
      RR0    = 3
    ),
    c.ga.tuning     = list(1,1,0.5),      # tuning parameter for MH algorithm (fixed effect parameters), each parameter corresponds to an outcome/marker
    c.theta.tuning  = list(1,1,0.5),      # tuning parameter for MH algorithm (random effect), each parameter corresponds to an outcome/marker
    adaptive.tuning = 1,                  # adaptive tuning parameters, 1 - yes, 0 - no
    align.clusters  = 1,                  # assign clusters
    tuning.freq     = 20,                 # tuning frequency
    initial.cluster.membership = "PAM",        # "mixAK" or "random" or "PAM" or "input" - input initial cluster membership for local clustering
    input.initial.cluster.membership = NULL,   # if use "input", option input.initial.cluster.membership must not be empty
    initial.global.cluster.membership = NULL,  # input initial cluster membership for global clustering
    seed.initial    = 2080,               # seed for initial clustering (for initial.cluster.membership = "mixAK")
    print.info      = "TRUE",             # print model information at each iteration
    burn.in,                              # number of samples discarded
    thin,                                 # thinning
    per,                                  # output information every "per" interation
    max.iter                              # maximum number of iteration
  ) {

  # removing NA values;
  R   <- length(mydat)
  dat <- vector(mode = "list", length = R)
  for (s in 1:R){
       id[[s]] <-    id[[s]][is.na(mydat[[s]])==FALSE]
     time[[s]] <-  time[[s]][is.na(mydat[[s]])==FALSE]
    mydat[[s]] <- mydat[[s]][is.na(mydat[[s]])==FALSE] # note the order, this line is last
      dat[[s]] <- data.frame(cbind(
                    y     = mydat[[s]],
                    time  =  time[[s]],
                    time2 =  time[[s]]^2,
                    time3 =  time[[s]]^3,
                    id    =    id[[s]]
                  ))
  }

  n.obs <- lapply(id, function(x) as.vector(table(x)))
  N     <- length(unique(id[[1]]))  # sample size should be identifcal across markers
                                    # number of measurements and time points can be different

  #--------------------------------------------------------------#
  # starting values;
  #--------------------------------------------------------------#
  theta <- vector(mode="list", length=R)
  k     <- vector(mode="list", length=R)  # dimension of fixed effect
  K     <- vector(mode="list", length=R)  # dimension of random effect
  cf    <- 0.1
  for (s in 1:R) {
    if (dist[[s]] == "gaussian") {
      fit.glmm <- lmer(
        formula[[s]],
        data    = dat[[s]],
        control = lmerControl(
          optimizer = "bobyqa",
          optCtrl   = list(maxfun=2e5)
        )
      )
    } else if (dist[[s]] == "poisson") {
      fit.glmm <- glmer(
        formula[[s]],
        data    = dat[[s]],
        nAGQ    = 0,
        family  = poisson(link = "log"),
        control = glmerControl(
          optimizer = "bobyqa",
          optCtrl   = list(maxfun=2e5)
        )
      )
    } else if (dist[[s]] == "binomial") {
        fit.glmm <- glmer(
          formula[[s]],
          data    = dat[[s]],
          nAGQ    = 0,
          family  = binomial(link = "logit"),
          control = glmerControl(
            optimizer = "bobyqa",
            optCtrl   = list(maxfun=2e5)
          )
        )
    }
        k[[s]] <- length(fixef(fit.glmm))
    theta[[s]] <- cf * data.matrix(ranef(fit.glmm)$id)
        K[[s]] <- dim(theta[[s]])[2]
  }

  if (length(initials) == 0) {        # use default initial values
    # For local cluster membership (outcome-specific) - use PAM method to initialize
    my.cluster     <- vector(mode = "list", length = R)
    my.cluster.tmp <- NULL
    for (s in 1:R) {
      if (initial.cluster.membership == "mixAK") {
        if (dist[[s]] == "gaussian") {
          #my.cluster[[s]] <- cluster::pam(theta[[s]], num.cluster)$clustering
          set.seed(seed.initial)
          fit.mixAK <- mixAK::GLMM_MCMC(
            y                = dat[[s]][,"y"],
            dist             = c("gaussian"),
            id               = dat[[s]][,"id"],
            z                = list(y = dat[[s]][,"time"]),
            random.intercept = c(TRUE),
            prior.b          = list(Kmax = num.cluster),
            parallel         = TRUE
            # ,silent = TRUE
          )
          fit.mixAK <- NMixRelabel(
            fit.mixAK,
            type           = "stephens",
            keep.comp.prob = TRUE
            # ,silent = TRUE
          )
          my.cluster[[s]] <- apply(fit.mixAK[[1]]$poster.comp.prob, 1, which.max)
        } else if (dist[[s]] == "poisson") {
          set.seed(seed.initial)
          fit.mixAK <- mixAK::GLMM_MCMC(
            y                = dat[[s]][,"y"],
            dist             = c("poisson(log)"),
            id               = dat[[s]][,"id"],
            z                = list(y = dat[[s]][,"time"]),
            random.intercept = c(TRUE),
            prior.b          = list(Kmax = num.cluster),
            parallel         = TRUE
            # ,silent = TRUE
          )
          fit.mixAK <- NMixRelabel(
            fit.mixAK,
            type           = "stephens",
            keep.comp.prob = TRUE
            # ,silent = TRUE
          )
          my.cluster[[s]] <- apply(fit.mixAK[[1]]$poster.comp.prob, 1, which.max)
        } else if (dist[[s]] == "binomial") {
          set.seed(seed.initial)
          fit.mixAK <- mixAK::GLMM_MCMC(
            y                = dat[[s]][,"y"],
            dist             = c("binomial(logit)"),
            id               = dat[[s]][,"id"],
            z                = list(y = dat[[s]][,"time"]),
            random.intercept = c(TRUE),
            prior.b          = list(Kmax = num.cluster),
            parallel         = TRUE
            # ,silent = TRUE
          )
          fit.mixAK <- NMixRelabel(
            fit.mixAK,
            type           = "stephens",
            keep.comp.prob = TRUE
            # ,silent = TRUE
          )
          my.cluster[[s]] <- apply(fit.mixAK[[1]]$poster.comp.prob, 1, which.max)
        }
      }
      if (initial.cluster.membership == "PAM")    {my.cluster[[s]] <- cluster::pam(theta[[s]],num.cluster)$clustering}
      if (initial.cluster.membership == "random") {my.cluster[[s]] <- sample(1:num.cluster,N,replace=TRUE)}
      if (initial.cluster.membership == "input")  {my.cluster[[s]] <- input.initial.cluster.membership[[s]]}

      my.cluster.tmp <- cbind(my.cluster.tmp, my.cluster[[s]])
      mydf.clust     <- data.frame(
        id=1:length(my.cluster[[s]]),
        my.cluster=my.cluster[[s]]
      )
      dat[[s]] <- merge(
        dat[[s]],
        mydf.clust,
        by="id"
      )
    }

    # For regression coefficients
    fixed.effect <- vector(mode="list", length=R)
    for (s in 1:R) {
      fixed.effect[[s]] <- matrix(0, nrow=num.cluster, ncol=k[[s]])
      for (j in 1:num.cluster) {
        if (dist[[s]] == "gaussian") {
          fit.glmm <- lmer(
            formula[[s]],
            data   = dat[[s]][dat[[s]]$my.cluster==j,]
          )
          fixed.effect[[s]][j,] <- fixef(fit.glmm)
        } else if (dist[[s]] == "poisson") {
          fixed.effect[[s]][j,] <- suppressMessages(fixef(glmer(
            formula[[s]],
            data   = dat[[s]][dat[[s]]$my.cluster==j,],
            nAGQ   = 0,
            family = poisson(link = "log")
          )))
        } else if (dist[[s]] == "binomial") {
          fixed.effect[[s]][j,] <- suppressMessages(fixef(glmer(
            formula[[s]],
            data   = dat[[s]][dat[[s]]$my.cluster==j,],
            nAGQ   = 0,
            family = binomial(link = "logit")
          )))
        }
      }
    }

    alpha    <- rep(0.9,R)
    # initial cluster membership
    # initial cluster membership
    if (length(initial.global.cluster.membership)==0)  {zz <-  my.cluster[[1]]}  else{zz <- initial.global.cluster.membership }

    #zz       <- my.cluster[[1]]            # can use lcmm here
    #zz.local <- vector(mode = "list", length = R)
    #for (s in 1:R) {
    #  zz.local[[s]] <- my.cluster[[s]]
    #}
    zz.local <- my.cluster

    # regression coeffcients
    ga <- fixed.effect

    # for residual variance (for gussian distribution only)
                                      sigma.sq.e      <- vector(mode = "list", length = R)
    for (s in 1:R) {
      for (j in 1:num.cluster) {
        if (dist[[s]] == "gaussian")  sigma.sq.e[[s]] <- rbind(sigma.sq.e[[s]], 1)
        if (dist[[s]] == "poisson")   sigma.sq.e[[s]] <- rbind(sigma.sq.e[[s]], NA)
        if (dist[[s]] == "binomial")  sigma.sq.e[[s]] <- rbind(sigma.sq.e[[s]], NA)
      }
    }

    # dispersion parameters
                                      phi      <- vector(mode = "list", length = R)
    for (s in 1:R) {
      for (j in 1:num.cluster) {
        if (dist[[s]] == "gaussian")  phi[[s]] <- rbind(phi[[s]], sigma.sq.e[[s]][j])
        if (dist[[s]] == "poisson")   phi[[s]] <- rbind(phi[[s]], 1)
        if (dist[[s]] == "binomial")  phi[[s]] <- rbind(phi[[s]], 1)
      }
    }

    # starting values for random effect variance
        sigma.sq.u               <- vector(mode = "list", length = R)
        sigma.sq.u.inv           <- vector(mode = "list", length = R)

    for (s in 1:R) {
        sigma.sq.u    [[s]]      <- array(0, c(K[[s]], K[[s]], num.cluster))
        sigma.sq.u.inv[[s]]      <- array(0, c(K[[s]], K[[s]], num.cluster))
      for (j in 1:num.cluster) {
        sigma.sq.u    [[s]][,,j] <- var(theta[[s]][my.cluster.tmp[,s]==j,1:K[[s]]]);
        sigma.sq.u.inv[[s]][,,j] <- solve(sigma.sq.u[[s]][,,j]);
      }
    }
  } else {  # use specified initials
    alpha      <- initials$alpha
    zz         <- initials$zz
    zz.local   <- initials$zz.local
    ga         <- initials$ga
    sigma.sq.e <- initials$sigma.sq.e

    # dispersion parameters
                                      phi      <- vector(mode = "list", length = R)
    for (s in 1:R) {
      for (j in 1:num.cluster) {
        if (dist[[s]] == "gaussian")  phi[[s]] <- rbind(phi[[s]], sigma.sq.e[[s]][j])
        if (dist[[s]] == "poisson")   phi[[s]] <- rbind(phi[[s]], 1)
        if (dist[[s]] == "binomial")  phi[[s]] <- rbind(phi[[s]], 1)
      }
    }

        sigma.sq.u               <- initials$sigma.sq.u
        sigma.sq.u.inv           <- vector(mode = "list", length = R)
    for (s in 1:R) {
        sigma.sq.u.inv[[s]]      <- array(0,c(K[[s]], K[[s]], num.cluster))
      for (j in 1:num.cluster) {
        sigma.sq.u.inv[[s]][,,j] <- solve(sigma.sq.u[[s]][,,j])
      }
    }

    # Here theta over-write the previous assigned values, if initials are given
    theta <- vector(mode = "list", length = R)
    for (s in 1:R) {
      theta[[s]] <- matrix(0, ncol=K[[s]], nrow=N)
      for (i in 1:N) {
        for (j in 1:num.cluster) {
          if (K[[s]] == 1) theta[[s]][i,1:K[[s]]] <- rmvnorm(n = 1, mean = rep(0,K[[s]]), sigma = matrix(sigma.sq.u[[s]][,,j]))
          if (K[[s]] > 1)  theta[[s]][i,1:K[[s]]] <- rmvnorm(n = 1, mean = rep(0,K[[s]]), sigma =        sigma.sq.u[[s]][,,j])
        }
      }
    }

  }


  ppi <- rep(1/num.cluster,num.cluster)    # for overall clustering
  # complete initial values to pass to Cpp function
  initials.complete <- list(ppi=ppi, alpha=alpha, zz=zz, zz.local=zz.local,
                            ga=ga, sigma.sq.e=sigma.sq.e,sigma.sq.u=sigma.sq.u,
                            theta=theta);
  #mean(zz.local[[3]]==simdata$cluster.local[[3]])
  #--------------------------------------------------------------#
  # Hyper-parameters for the Prior Distributions
  #--------------------------------------------------------------#
  delta  <- hyper.par$delta
  a.star <- hyper.par$a.star;  b.star <- hyper.par$b.star
  #---- hyper-parameters for the residual variances - common to both dataset
  aa0    <- hyper.par$aa0; bb0 <- hyper.par$bb0
  a0     <- rep(list(rep(aa0,num.cluster)),R)
  b0     <- rep(list(rep(bb0,num.cluster)),R)
  #----- hyper-parameters for fixed effect variables - common to both dataset
  w0     <- vector(mode = "list", length = R)
  omega0 <- vector(mode = "list", length = R)
  #---- hyper-parameters for the random effect variances
  cc0    <- hyper.par$cc0; dd0 <- hyper.par$dd0
  c0     <- rep(list(rep(cc0,num.cluster)),R)
  d0     <- rep(list(rep(dd0,num.cluster)),R)

  #----  hyper-parameters for Wishart Distribution
  rr0    <- hyper.par$rr0;
  RR0    <- hyper.par$RR0
  ww0    <- hyper.par$ww0;
  vv0    <- hyper.par$vv0;
  r0     <- vector(mode = "list", length = R)
  R0     <- vector(mode = "list", length = R)
  for (s in 1:R) {
    for (j in 1:num.cluster) {
           q      <- length(ga[[s]][j,])
          w0[[s]] <- matrix(ww0, nrow=num.cluster, ncol=q)
      omega0[[s]] <- array(diag(vv0,q), dim=c(q,q,num.cluster))
          r0[[s]] <- rep(rr0, num.cluster)
          R0[[s]] <- array(diag(RR0,K[[s]]), c(K[[s]],K[[s]],num.cluster))
    }
  }

  #--------------------------------------------------------------#
  #--------------------------------------------------------------#
  #--------------------------------------------------------------#
  #--------------------------------------------------------------#
  # Storing the sample;
  LOG.LIK.ITER <- NULL
  PPI          <- NULL
  ZZ           <- NULL
  ALPHA        <- NULL
  ZZ.LOCAL     <- vector(mode = "list", length = R)
  GA           <- vector(mode = "list", length = R)
  GA.ACCEPT    <- vector(mode = "list", length = R)
  THETA        <- vector(mode = "list", length = R)
  THETA.ACCEPT <- vector(mode = "list", length = R)
  SIGMA.SQ.U   <- vector(mode = "list", length = R)
  SIGMA.SQ.E   <- vector(mode = "list", length = R)
  T.LOCAL      <- vector(mode = "list", length = R)
  T            <- NULL
  for (s in 1:R) {
        ZZ.LOCAL[[s]] <- matrix(0, nrow=(max.iter-burn.in)/thin, ncol=N)
      SIGMA.SQ.E[[s]] <- matrix(0, nrow=(max.iter-burn.in)/thin, ncol=num.cluster)
           THETA[[s]] <-  array(0, c(N,K[[s]],      (max.iter-burn.in)/thin))
         T.LOCAL[[s]] <-  array(0, c(N,num.cluster, (max.iter-burn.in)/thin))
       GA.ACCEPT[[s]] <- matrix(0, nrow=(max.iter-burn.in)/thin, ncol=num.cluster)
    THETA.ACCEPT[[s]] <- matrix(0, nrow=(max.iter-burn.in)/thin, ncol=N)
  }

  cat(paste(rep('-',60),sep='',collapse=''), '\n'); cat(paste(rep('-',60),sep='',collapse=''), '\n');
  cat('Running BCC Model (GLM)', '\n')
  cat(paste(rep('-',60),sep='',collapse=''), '\n'); cat(paste(rep('-',60),sep='',collapse=''), '\n');

  c.ga <- vector(mode = "list", length = R)
  for (s in 1:R) {
    c.ga[[s]] <- rep(c.ga.tuning[[s]],num.cluster)
  }
  c.theta <- c.theta.tuning

  begin <- proc.time()[1]

  sourceCpp("BCC.cpp")
  tryCatch({
    rst = BCC(
      dat, R,
      id, simplify2array(n.obs), N,
      num.cluster,
      dist,
      alpha.common,
      sigma.sq.e.common,
      align.clusters,
      unlist(k), unlist(K),
      sig.var,
      # initials
        ppi,
        alpha,
        zz,
        t(simplify2array(zz.local)),
        ga,
        lapply(sigma.sq.e, function(x) {if (is.null(x)) {numeric()} else {x}}),
        phi,
        sigma.sq.u,
        theta,
      # Hyper-parameters
        delta,
        a.star,
        b.star,
        aa0,
        bb0,
        t(simplify2array(a0)),
        t(simplify2array(b0)),
        w0,
        omega0,
        cc0,
        dd0,
        t(simplify2array(c0)),
        t(simplify2array(d0)),
        rr0,
        RR0,
        ww0,
        vv0,
        t(simplify2array(r0)),
        R0,
      # sample
        LOG.LIK.ITER,
        PPI,
        ZZ,
        ALPHA,
        ZZ.LOCAL,
        GA,
        GA.ACCEPT,
        THETA,
        THETA.ACCEPT,
        SIGMA.SQ.U,
        SIGMA.SQ.E,
        T.LOCAL,
        T,
      adaptive.tuning,
      tuning.freq,
      t(simplify2array(c.ga)),
      unlist(c.theta),
      burn.in,
      thin,
      per,
      max.iter,
      seed.initial
    )},
    error=function(cond) {
      message("Here's the original error message:")
      message(cond$message)
      # Choose a return value in case of error
      return(NULL)
    }
  )

  end = proc.time()[1]
  cat('It took', end - begin, 'seconds\n')
  run.time <-  end - begin

  rst$PPI   <- matrix(rst$PPI,   ncol=num.cluster, byrow=TRUE)
  rst$ZZ    <- matrix(rst$ZZ,    ncol=N,           byrow=TRUE)
  rst$ALPHA <- matrix(rst$ALPHA, ncol=R,           byrow=TRUE)
  for (s in 1:R) {
    rst$SIGMA.SQ.E  [[s]] <- matrix(rst$SIGMA.SQ.E  [[s]],ncol=num.cluster)
    rst$GA.ACCEPT   [[s]] <- matrix(rst$GA.ACCEPT   [[s]],ncol=num.cluster)
    rst$THETA.ACCEPT[[s]] <- matrix(rst$THETA.ACCEPT[[s]],ncol=N)
    rst$THETA       [[s]] <-  array(rst$THETA       [[s]],c(N, K[[s]], length(rst$THETA[[s]])/(N*K[[s]])))
    rst$ZZ.LOCAL    [[s]] <- matrix(rst$ZZ.LOCAL    [[s]],ncol=N)
    rst$T.LOCAL     [[s]] <-  array(rst$T.LOCAL     [[s]], c(N,num.cluster, length(rst$T.LOCAL[[s]])/(N*num.cluster)))
  }
  dimnames(rst$ZZ)    <- dimnames(ZZ)
  dimnames(rst$ALPHA) <- dimnames(ALPHA)
  PPI          <- rst$PPI
  ZZ           <- rst$ZZ
  T            <- rst$T
  ALPHA        <- rst$ALPHA
  SIGMA.SQ.E   <- rst$SIGMA.SQ.E
  GA.ACCEPT    <- rst$GA.ACCEPT
  THETA.ACCEPT <- rst$THETA.ACCEPT
  THETA        <- rst$THETA
  ZZ.LOCAL     <- rst$ZZ.LOCAL
  T.LOCAL      <- rst$T.LOCAL
  SIGMA.SQ.U   <- rst$SIGMA.SQ.U
  GA           <- rst$GA
  iter         <- rst$iter

  #--------------------------------------------------------------------------------------------------#
  cat(paste(rep('-',60),sep='',collapse=''), '\n'); cat(paste(rep('-',60),sep='',collapse=''), '\n');
  cat('Post-Processing Results', '\n')
  cat(paste(rep('-',60),sep='',collapse=''), '\n'); cat(paste(rep('-',60),sep='',collapse=''), '\n');
  #--------------------------------------------------------------------------------------------------#
  #----------------------------------------------------------------------------#
  #- Apply burn.in and thin
  num.sample <- length(seq(burn.in + 1,iter,thin))
  #----------------------------------------------------------------------------#
  # Address Label Switching Using Stephens' algorithm
  #----------------------------------------------------------------------------#
  if (num.cluster > 1) {
    # Apply Stephens' algorithm
    # for global cluster membership
    T.trans <- array(0,c(num.sample,N,num.cluster))
    for (j in 1:num.cluster){T.trans[,,j] <- t(T[,j,])}
    out.relabel <- label.switching(method="STEPHENS",z=ZZ,K=num.cluster,p=T.trans)$permutations$STEPHENS
    tp1 <- apply(T,c(1,2),mean); tp2 <- apply(tp1,1,sum)
    tp  <- cbind(tp1,tp2)
    tp[,1:num.cluster] <- tp[,1:num.cluster]/tp[,(num.cluster+1)]  # standardize
    postprob <- apply(tp[,1:num.cluster],1,max)
    #hist(postprob)

    # for local cluster membership
    T.LOCAL.trans     <- vector(mode = "list", length = R)
    out.relabel.local <- vector(mode = "list", length = R)
    for (s in 1:R){
      T.LOCAL.trans[[s]] <- array(0,c(num.sample,N,num.cluster))
      for (j in 1:num.cluster){T.LOCAL.trans[[s]][,,j] <- t(T.LOCAL[[s]][,j,])}
      out.relabel.local[[s]] <- label.switching(method="STEPHENS",z=ZZ.LOCAL[[s]],K=num.cluster,p=T.LOCAL.trans[[s]])$permutations$STEPHENS
    }

    # s
    # table(ZZ.LOCAL[[s]])
    # Post-processing the parameters according to swicthing label
    for (s in 1:R){
      for (j in 1:num.sample) {
        SIGMA.SQ.E[[s]][j,]  <- SIGMA.SQ.E[[s]][j, out.relabel.local[[s]][j,]     ]
        SIGMA.SQ.U[[s]][,,j] <- SIGMA.SQ.U[[s]][ , out.relabel.local[[s]][j,], j  ]
                GA[[s]][,,j] <-         GA[[s]][   out.relabel.local[[s]][j,],  ,j]
           T.LOCAL[[s]][,,j] <-    T.LOCAL[[s]][ , out.relabel.local[[s]][j,], j  ]
      }
    }

    # table(ZZ.LOCAL[[s]])

    # Compute the global and local cluster membership
    cluster.global <- apply(apply(T,c(1,2),mean),1,nnet::which.is.max)
    cluster.local  <- vector(mode = "list", length = R)
    for (s in 1:R) {
      cluster.local[[s]] <- apply(apply(T.LOCAL[[s]],c(1,2),mean),1,nnet::which.is.max)
      mycluster <- data.frame(id=unique(id[[s]]),cluster.global=cluster.global,cluster.local=cluster.local[[s]])
      dat[[s]] <- merge(dat[[s]],mycluster,by="id")
    }

    #--- adjusted adherance---
    my.alpha        <- apply(ALPHA,2,mean); my.alpha
    my.alpha.adjust <- (num.cluster*my.alpha - 1)/(num.cluster-1)
       alpha.adjust <- mean(my.alpha.adjust)
  }
  #-------------------------------------------------------------------------------------------#
  # Calculate Summary Statistics for Model Parameters (mean, sd, 95%CR and geweke statistics)
  #-------------------------------------------------------------------------------------------#
  res <- function(x) {c(mean=mean(x),sd=sd(x),quantile(x,c(0.025,0.975)), geweke.stat=as.vector(geweke.diag(x)$z[1]))}
         PPI.stat <- apply(PPI,2,res)
       ALPHA.stat <- apply(ALPHA,2,res)
  SIGMA.SQ.E.stat <- vector(mode = "list", length = R)
  SIGMA.SQ.U.stat <- vector(mode = "list", length = R)
  GA.stat <- vector(mode = "list", length = R)
  for (s in 1:R) {
    if (dist[[s]] == "gaussian") {
      SIGMA.SQ.E.stat[[s]] <- apply(SIGMA.SQ.E[[s]], 2,      res)
    }
      SIGMA.SQ.U.stat[[s]] <- apply(SIGMA.SQ.U[[s]], c(1,2), res)
              GA.stat[[s]] <- apply(GA[[s]],         c(1,2), res)
  }
  summary.stat <- list(
    PPI        = PPI.stat,
    ALPHA      = ALPHA.stat,
    GA         = GA.stat,
    SIGMA.SQ.U = SIGMA.SQ.U.stat,
    SIGMA.SQ.E = SIGMA.SQ.E.stat)
  summary.stat

  if (num.cluster == 1) {
    postprob <- cluster.global <- cluster.local <- PPI <- T <-
    ALPHA <- my.alpha <- alpha.adjust <- THETA.ACCEPT <- GA.ACCEPT <-  NULL;
  }
  # returning the parameters;
  list(
    dat            = dat,
    PPI            = PPI,
    ZZ             = ZZ,
    T              = T,
    ALPHA          = ALPHA,
    SIGMA.SQ.E     = SIGMA.SQ.E,
    SIGMA.SQ.U     = SIGMA.SQ.U,
    T.LOCAL        = T.LOCAL,
    ZZ.LOCAL       = ZZ.LOCAL,
    GA             = GA,
    THETA.ACCEPT   = THETA.ACCEPT,
    GA.ACCEPT      = GA.ACCEPT,
    alpha          = my.alpha,
    alpha.adjust   = alpha.adjust,
    postprob       = postprob,
    cluster.global = cluster.global,
    cluster.local  = cluster.local,
    k              = k,
    K              = K,
    sig.var        = sig.var,
    N              = N,
    R              = R,
    dist           = dist,
    num.cluster    = num.cluster,
    THETA          = THETA,
    #tot.num.par   = tot.num.par,
    max.iter       = max.iter,
    burn.in        = burn.in,
    thin           = thin,
    run.time       = run.time,
    summary.stat   = summary.stat)
}



# - AlignCllusters function - adapted from the codes of BCC original paper (Lock and Dunson 2013)
# currently unused function
AlignClusters <- function(Z1,Z2, type = 'vec') {
  if(type == 'vec') {
    for(k in 1:length(unique(Z1))) {
      Max = sum(Z1==k & Z2==k)/(.01+sum(Z2==k)+sum(Z1==k));
      for(tempk in  1:length(unique(Z2))) {
        if(sum(Z1==k & Z2==tempk)/(.01+sum(Z2==tempk)+sum(Z1==k)) > Max) {
          Max           = sum(Z1==k &Z2==tempk)/(.01+sum(Z2==tempk)+sum(Z1==k))
          dummy         = Z2==k
          Z2[Z2==tempk] = k
          Z2[dummy]     = tempk
        }
      }
    }
  } else if(type == 'mat') {
    for(k in 1:dim(Z1)[2]) {
      for(tempk in  1:dim(Z2)[2]) {
        Max             = sum(Z1==Z2)
        Z2dummy         = Z2
        Z2dummy[,k]     = Z2[,tempk]
        Z2dummy[,tempk] = Z2[,k]
        if(sum(Z1==Z2dummy)>Max)
          Z2 = Z2dummy
      }
    }
  }
  return(Z2)
}

#library(compiler)
#BCC.multic <- cmpfun(BCC.multi)

# [END]
