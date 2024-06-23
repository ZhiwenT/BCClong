#' Compute a Bayesian Consensus Clustering model for mixed-type longitudinal data
#'
#' This function performs clustering on mixed-type (continuous, discrete and
#' categorical) longitudinal markers using Bayesian consensus clustering method
#' with MCMC sampling
#'
#' @param mydat list of R longitudinal features (i.e., with a length of R),
#'              where R is the number of features. The data should be prepared
#'              in a long-format (each row is one time point per individual).
#' @param id a list (with a length of R) of vectors of the study id of
#'           individuals for each feature. Single value (i.e., a length of 1)
#'           is recycled if necessary
#' @param time a list (with a length of R) of vectors of time (or age) at which
#'             the feature measurements are recorded
#' @param center 1: center the time variable before clustering, 0: no centering
#' @param num.cluster number of clusters K
#' @param formula a list (with a length of R) of formula for each feature.
#'                Each formula is a twosided linear formula object describing
#'                both the fixed-effects and random effects part of the model,
#'                with the response (i.e., longitudinal feature) on the left
#'                of a ~ operator and the terms, separated by + operations,
#'                or the right. Random-effects terms are distinguished by
#'                vertical bars (|) separating expressions for design matrices
#'                from grouping factors.
#'                See formula argument from the lme4 package
#' @param dist a character vector (with a length of R) that determines the
#'             distribution for each feature. Possible values are "gaussian"
#'             for a continuous feature, "poisson" for a discrete feature
#'             (e.g., count data) using a log link and "binomial" for a
#'             dichotomous feature (0/1) using a logit link. Single value
#'             (i.e., a length of 1) is recycled if necessary
#' @param alpha.common 1 - common alpha, 0 - separate alphas for each outcome
#' @param initials List of initials for: zz, zz.local ga, sigma.sq.u, sigma.sq.e,
#'                 Default is NULL
#' @param sigma.sq.e.common 1 - estimate common residual variance across all groups,
#'                          0 - estimate distinct residual variance, default is 1
#' @param hyper.par hyper-parameters of the prior distributions for the model
#'                  parameters. The default hyper-parameters values will result
#'                  in weakly informative prior distributions.
#' @param c.ga.tunning tuning parameter for MH algorithm (fixed effect parameters),
#'                    each parameter corresponds to an outcome/marker, default
#'                    value equals NULL
#' @param c.theta.tunning tuning parameter for MH algorithm (random effect),
#'                       each parameter corresponds to an outcome/marker,
#'                       default value equals NULL
#' @param adaptive.tunning adaptive tuning parameters, 1 - yes, 0 - no,
#'                        default is 1
#' @param tunning.freq     tuning frequency, default is 20
#' @param initial.cluster.membership "mixAK" or "random" or "PAM" or "input" -
#'                                  input initial cluster membership for local
#'                                  clustering, default is "random"
#' @param input.initial.local.cluster.membership if use "input",
#'                                  option input.initial.cluster.membership
#'                                  must not be empty, default is NULL
#' @param input.initial.global.cluster.membership input initial cluster
#'                                  membership for global clustering
#'                                  default is NULL
#' @param seed.initial seed for initial clustering
#'                    (for initial.cluster.membership = "mixAK")
#'                    default is 2080
#' @param burn.in the number of samples disgarded.
#'                This value must be smaller than max.iter.
#' @param thin the number of thinning. For example, if thin = 10,
#'             then the MCMC chain will keep one sample every 10 iterations
#' @param per specify how often the MCMC chain will print the iteration number
#' @param max.iter the number of MCMC iterations.
#' @return Returns a BCC class model contains clustering information
#' @examples
#' # import dataframe
#' data(epil)
#' # example only, larger number of iteration required for accurate result
#' fit.BCC <-  BCC.multi (
#'        mydat = list(epil$anxiety_scale,epil$depress_scale),
#'        dist = c("gaussian"),
#'        id = list(epil$id),
#'        time = list(epil$time),
#'        formula =list(y ~ time + (1|id)),
#'        num.cluster = 2,
#'        burn.in = 3,
#'        thin = 1,
#'        per =1,
#'        max.iter = 8)
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
#' @import abind
#' @importFrom coda geweke.diag
#' @importFrom stats binomial poisson sd var
#' @importFrom utils capture.output
#' @useDynLib BCClong, .registration=TRUE

BCC.multi <- function(
    mydat,                 # List of R outcomes (R is the number of outcome, R>=2)
    id,                    # id-variable: starting from 1 to N
    time,                  # time variable
    center=1,              # 1 - center the time variable before clustering, 0 - no centering
    num.cluster,           # number of cluster
    formula,               # fixed and random effect
    dist,                  # "gaussian","poisson","binomial", distribution of the marker
    alpha.common = 0,          # 1 - common alpha, 0 - separate alphas for each marker
    initials = NULL,       # List of initials for: zz, zz.local ga, sigma.sq.u, sigma.sq.e,
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
    c.ga.tunning     = NULL,      # tunning parameter for MH algorithm (fixed effect parameters), each parameter corresponds to an outcome/marker
    c.theta.tunning  = NULL,      # tunning parameter for MH algorithm (random effect), each parameter corresponds to an outcome/marker
    adaptive.tunning = 0,                  # adaptive tunning parameters, 1 - yes, 0 - no
    tunning.freq     = 20,                 # tunning frequency
    initial.cluster.membership = "random",        # "mixAK" or "random" or "input" - input initial cluster membership for local clustering
    input.initial.local.cluster.membership = NULL,   # if use "input", option input.initial.local.cluster.membership must not be empty
    input.initial.global.cluster.membership = NULL,  # input initial cluster membership for global clustering
    seed.initial    = 2080,               # seed for initial clustering (for initial.cluster.membership = "mixAK")
    burn.in,                              # number of samples discarded
    thin,                                 # thinning
    per,                                  # output information every "per" interaction
    max.iter                              # maximum number of iteration
) {

  #-------------------------------------------------------------------------------------#
  # Set up
  create.new.id <- function(input_id){
    # Create new ID from 1 to N;
    subj <- unique(input_id)
    N <- length(subj)
    id.new <- NULL
    for (i in 1:N) {id.new   <- c(id.new,rep(i,length(input_id[input_id==subj[i]])))}
    return(id.new)}
  #-------------------------------------------------------------------------------------#

  R <- length(mydat)
  if (length(dist)==1) dist = rep(dist,R)
  if (length(id)==1) id = rep(id,R)
  if (length(time)==1) time = rep(time,R)
  if (length(formula)==1) formula = rep(formula,R)
  if (is.null(c.ga.tunning)==TRUE)  c.ga.tunning <- rep(list(1),R)
  if (is.null(c.theta.tunning)==TRUE)  c.theta.tunning <- rep(list(1),R)

  # removing NA values;
  dat <- vector(mode = "list", length = R)
  time.org <-  vector(mode = "list", length = R)
  for (s in 1:R){
    id[[s]] <-    id[[s]][is.na(mydat[[s]])==FALSE]
    time.org[[s]] <-  time[[s]][is.na(mydat[[s]])==FALSE]
    if (center==1){ time[[s]] <-  time[[s]][is.na(mydat[[s]])==FALSE]; time[[s]] <- time[[s]] - mean(time[[s]]) }
    mydat[[s]] <- mydat[[s]][is.na(mydat[[s]])==FALSE] # note the order, this line is last
  }

  # Find common id
  # (require each individual to have at least one observation for all markers)
  common.id <- NULL
  for (s in 1:R){
    if (s==1)  common.id <- unique(id[[s]]) else{
      common.id <- Reduce(intersect,list(common.id, unique(id[[s]])))}}
  common.id <- common.id[is.na(common.id)==FALSE]
  N <- length(common.id); N 			# number of individuals included in the analysis

  #---------------------------------------------#
  id.org <-  vector(mode = "list", length = R)
  id.new <-  vector(mode = "list", length = R)
  for (s in 1:R){
    id.org[[s]] <-    id[[s]][id[[s]] %in% common.id]; length(id.org[[1]])
    id.new[[s]] <- create.new.id(id.org[[s]]); length(id.new[[1]])
    time.org[[s]] <-  time.org[[s]][id[[s]] %in% common.id]; length(time.org[[1]])
    time[[s]] <-  time[[s]][id[[s]] %in% common.id]; length(time[[1]])
    mydat[[s]] <- mydat[[s]][id[[s]] %in% common.id]  # note the order, this line is last
    dat[[s]] <- data.frame(cbind(
      y     = mydat[[s]],
      time.org  =  time.org[[s]],
      time  =  time[[s]],
      time2 =  time[[s]]^2,
      time3 =  time[[s]]^3,
      id.org    = id.org[[s]],
      id    = id.new[[s]]))
  }

  n.obs <- lapply(id.org, function(x) as.vector(table(x)))  # number of measurements and time points can be different

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

  if (length(initials) == 0) {        			# use default initial values
    my.cluster     <- vector(mode = "list", length = R)
    my.cluster.tmp <- NULL
    for (s in 1:R) {
      if (initial.cluster.membership == "mixAK") {
        if (dist[[s]] == "gaussian") {
          set.seed(seed.initial)
          fit.mixAK <- mixAK::GLMM_MCMC(
            y                = dat[[s]][,"y"],
            dist             = c("gaussian"),
            id               = dat[[s]][,"id"],
            z                = list(y = dat[[s]][,"time"]),
            random.intercept = c(TRUE),
            prior.b          = list(Kmax = num.cluster),
            parallel         = TRUE
            ,silent = TRUE
          )
          fit.mixAK <- NMixRelabel(
            fit.mixAK,
            type           = "stephens",
            keep.comp.prob = TRUE
            ,silent = TRUE
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
            ,silent = TRUE
          )
          fit.mixAK <- NMixRelabel(
            fit.mixAK,
            type           = "stephens",
            keep.comp.prob = TRUE
            ,silent = TRUE
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
            ,silent = TRUE
          )
          fit.mixAK <- NMixRelabel(
            fit.mixAK,
            type           = "stephens",
            keep.comp.prob = TRUE
            ,silent = TRUE
          )
          my.cluster[[s]] <- apply(fit.mixAK[[1]]$poster.comp.prob, 1, which.max)
        }
      }
      if (initial.cluster.membership == "random") {my.cluster[[s]] <- sample(1:num.cluster,N,replace=TRUE)}
      if (initial.cluster.membership == "input")  {my.cluster[[s]] <- input.initial.local.cluster.membership[[s]]}

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
    if (length(input.initial.global.cluster.membership)==0)  {zz <-  my.cluster[[1]]}  else{zz <- input.initial.global.cluster.membership }
    zz.local <- my.cluster

    # regression coeffcients
    ga <- fixed.effect

    # for residual variance (for gaussian distribution only)
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
        solve.tmp <- try(solve(sigma.sq.u[[s]][,,j]), silent=TRUE)
        if (inherits(solve.tmp,"try-error")==FALSE)  {sigma.sq.u.inv[[s]][,,j] <- solve.tmp} else{sigma.sq.u.inv[[s]][,,j] <- 1e-5 }
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

    # Here theta over-write the previously assigned values, if initials are given
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

  #--------------------------------------------------------------#
  # Hyper-parameters for the Prior Distributions
  #--------------------------------------------------------------#
  if (length(hyper.par$delta) == 1) {delta  <- rep(hyper.par$delta,num.cluster)} else{delta = hyper.par$delta}
  a.star <- hyper.par$a.star;  b.star <- hyper.par$b.star
  #---- hyper-parameters for the residual variances - common to both datasets
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

  message(paste(rep('-',60),sep='',collapse='') ); message(paste(rep('-',60),sep='',collapse=''));
  message('Running BCC Model')
  message(paste(rep('-',60),sep='',collapse='')); message(paste(rep('-',60),sep='',collapse=''));

  c.ga <- vector(mode = "list", length = R)
  for (s in 1:R) {
    c.ga[[s]] <- rep(c.ga.tunning[[s]],num.cluster)
  }
  c.theta <- c.theta.tunning

  begin <- proc.time()[1]
  #sourceCpp("BCC.cpp")
  tryCatch({
    rst = BCC(
      dat, R,
      id, simplify2array(n.obs), N,
      num.cluster,
      dist,
      alpha.common,
      sigma.sq.e.common,
      unlist(k), unlist(K),
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
      adaptive.tunning,
      tunning.freq,
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
  message('It took ', end - begin, ' seconds')
  run.time <-  end - begin

  rst$PPI   <- matrix(rst$PPI,   ncol=num.cluster, byrow=TRUE)
  rst$ZZ    <- matrix(rst$ZZ,    ncol=N,           byrow=TRUE)
  if(num.cluster > 1) rst$ALPHA <- matrix(rst$ALPHA, ncol=R, byrow=TRUE) else{rst$ALPHA <- matrix(rst$ALPHA, ncol=1, byrow=TRUE);}
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
  message(paste(rep('-',60),sep='',collapse='')); message(paste(rep('-',60),sep='',collapse=''));
  message('Post-Processing Results')
  message(paste(rep('-',60),sep='',collapse='')); message(paste(rep('-',60),sep='',collapse=''));
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
    invisible(capture.output(out.relabel <- label.switching(method="STEPHENS",
                                                            z=ZZ,K=num.cluster,
                                                            p=T.trans)$permutations$STEPHENS))
    tp1 <- apply(T,c(1,2),mean); tp2 <- apply(tp1,1,sum)
    tp  <- cbind(tp1,tp2)
    tp[,1:num.cluster] <- tp[,1:num.cluster]/tp[,(num.cluster+1)]  # standardize
    postprob <- apply(tp[,1:num.cluster],1,max)


    # for local cluster membership
    T.LOCAL.trans     <- vector(mode = "list", length = R)
    out.relabel.local <- vector(mode = "list", length = R)
    for (s in 1:R){
      T.LOCAL.trans[[s]] <- array(0,c(num.sample,N,num.cluster))
      for (j in 1:num.cluster){T.LOCAL.trans[[s]][,,j] <- t(T.LOCAL[[s]][,j,])}
      invisible(capture.output(out.relabel.local[[s]] <-
                                 label.switching(method="STEPHENS",
                                                 z=ZZ.LOCAL[[s]],
                                                 K=num.cluster,
                                                 p=T.LOCAL.trans[[s]])$permutations$STEPHENS))
    }

    # Post-processing the parameters according to the switching label
    for (s in 1:R){
      for (j in 1:num.sample) {
        SIGMA.SQ.E[[s]][j,]  <- SIGMA.SQ.E[[s]][j, out.relabel.local[[s]][j,]     ]
        SIGMA.SQ.U[[s]][,,j] <- SIGMA.SQ.U[[s]][ , out.relabel.local[[s]][j,], j  ]
        GA[[s]][,,j] <-         GA[[s]][   out.relabel.local[[s]][j,],  ,j]
        T.LOCAL[[s]][,,j] <-    T.LOCAL[[s]][ , out.relabel.local[[s]][j,], j  ]
      }
    }

    # Compute the global and local cluster membership
    cluster.global <- apply(apply(T,c(1,2),mean),1,nnet::which.is.max)
    cluster.local  <- vector(mode = "list", length = R)
    for (s in 1:R) {
      cluster.local[[s]] <- apply(apply(T.LOCAL[[s]],c(1,2),mean),1,nnet::which.is.max)
      mycluster <- data.frame(id=1:N,cluster.global=cluster.global,cluster.local=cluster.local[[s]])
      dat[[s]] <- merge(dat[[s]],mycluster,by="id")
    }

    #--- adjusted adherence---
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
  #summary.stat

  if (num.cluster == 1) {
    postprob <- cluster.global <- cluster.local <- PPI <- T <-
      ALPHA <- my.alpha <- alpha.adjust <- THETA.ACCEPT <- GA.ACCEPT <-  NULL;
  }


  # setting class to the returning objects
  class(dat) <- "data"
  class(N) <- "data"
  class(R) <- "data"
  class(PPI) <- "MCMC_sample"
  class(ZZ) <- "MCMC_sample"
  class(ALPHA) <- "MCMC_sample"
  class(SIGMA.SQ.E) <- "MCMC_sample"
  class(SIGMA.SQ.U) <- "MCMC_sample"
  class(ZZ.LOCAL) <- "MCMC_sample"
  class(GA) <- "MCMC_sample"
  class(THETA.ACCEPT) <- "MCMC_sample"
  class(GA.ACCEPT) <- "MCMC_sample"
  class(my.alpha) <- "model_parameter"
  class(alpha.adjust) <- "model_parameter"
  class(postprob) <- "model_parameter"
  class(k) <- "model_parameter"
  class(K) <- "model_parameter"
  class(dist) <- "model_parameter"
  class(num.cluster) <- "model_parameter"
  class(THETA) <- "model_parameter"
  class(cluster.global) <- "cluster_membership"
  class(cluster.local) <- "cluster_membership"
  class(max.iter) <- "algorithm_parameter"
  class(burn.in) <- "algorithm_parameter"
  class(thin) <- "algorithm_parameter"
  class(run.time) <- "algorithm_parameter"
  class(summary.stat) <- "summary_statistics"

  # returning the parameters;
  res <- list(
    dat            = dat,
    N              = N,
    R              = R,
    PPI            = PPI,
    ZZ             = ZZ,
    ALPHA          = ALPHA,
    SIGMA.SQ.E     = SIGMA.SQ.E,
    SIGMA.SQ.U     = SIGMA.SQ.U,
    #T.LOCAL        = T.LOCAL,
    ZZ.LOCAL       = ZZ.LOCAL,
    GA             = GA,
    THETA.ACCEPT   = THETA.ACCEPT,
    GA.ACCEPT      = GA.ACCEPT,
    alpha          = my.alpha,
    alpha.adjust   = alpha.adjust,
    postprob       = postprob,
    k              = k,
    K              = K,
    dist           = dist,
    num.cluster    = num.cluster,
    THETA          = THETA,
    cluster.global = cluster.global,
    cluster.local  = cluster.local,
    max.iter       = max.iter,
    burn.in        = burn.in,
    thin           = thin,
    run.time       = run.time,
    summary.stat   = summary.stat)
  class(res) <- "BCC"
  return(res)
}

#library(compiler)
#BCC.multic <- cmpfun(BCC.multi)

# [END]
