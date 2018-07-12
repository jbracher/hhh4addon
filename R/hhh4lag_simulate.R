################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Simulate from a HHH4 model
###
### Copyright (C) 2012 Michaela Paul, 2013-2016,2018 Sebastian Meyer
### $Revision$
### $Date$
################################################################################


### Simulate-method for hhh4-objects

simulate.hhh4lag <- function (object, # result from a call to hhh4
                              nsim=1, # number of replicates to simulate
                              seed=NULL,
                              y.start=NULL, # initial counts for epidemic components
                              subset=1:nrow(object$stsObj),
                              coefs=coef(object), # coefficients used for simulation
                              components=c("ar","ne","end"), # which comp to include
                              simplify=nsim>1, # counts array only (no full sts)
                              ...)
{
  ## Determine seed (this part is copied from stats:::simulate.lm with
  ## Copyright (C) 1995-2012 The R Core Team)
  if(!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE))
    runif(1)                     # initialize the RNG if necessary
  if(is.null(seed))
    RNGstate <- get(".Random.seed", envir = .GlobalEnv)
  else {
    R.seed <- get(".Random.seed", envir = .GlobalEnv)
    set.seed(seed)
    RNGstate <- structure(seed, kind = as.list(RNGkind()))
    on.exit(assign(".Random.seed", R.seed, envir = .GlobalEnv))
  }
  ## END seed

  cl <- match.call()
  theta <- if (missing(coefs)) coefs else checkCoefs(object, coefs)
  control <- object$control #BJ

  ## lags
  lag.ar <- object$control$ar$lag
  lag.ne <- object$control$ne$lag
  # maxlag <- max(lag.ar, lag.ne) #BJ
  maxlag <- max(c(control$max_lag, control$ar$lag, control$ne$lag), na.rm = TRUE) #BJ
  minlag <- min(c(control$min_lag, control$ar$lag, control$ne$lag), na.rm = TRUE)

  ## initial counts
  nUnits <- object$nUnit
  if (is.null(y.start)) { # set starting value to mean observed (in subset!)
    y.means <- ceiling(colMeans(observed(object$stsObj)[subset,,drop=FALSE]))
    y.start <- matrix(y.means, maxlag, nUnits, byrow=TRUE)
  } else {
    if (is.vector(y.start)) y.start <- t(y.start)
    if (ncol(y.start) != nUnits)
      stop(sQuote("y.start"), " must have nUnits=", nUnits, " columns")
    if (nrow(y.start) < maxlag)
      stop("need 'y.start' values for lag=", maxlag, " initial time points")
  }

  ## get fitted exppreds nu_it, phi_it, lambda_it (incl. offsets, t in subset)
  exppreds <- get_exppreds_with_offsets_lag(object, subset = subset, theta = theta) #BJ

  ## extract overdispersion parameters (simHHH4 assumes psi->0 means Poisson)
  model <- terms.hhh4lag(object) #BJ
  psi <- surveillance:::splitParams(theta,model)$overdisp #BJ
  if (length(psi) > 1) # "NegBinM" or shared overdispersion parameters
    psi <- psi[model$indexPsi]

  ## weight matrix/array of the ne component
  neweights <- surveillance:::getNEweights(object, coefW(theta))

  ## set predictor to zero if not included ('components' argument)
  stopifnot(length(components) > 0, components %in% c("ar", "ne", "end"))
  getComp <- function (comp) {
    exppred <- exppreds[[comp]]
    if (comp %in% components) exppred else "[<-"(exppred, value = 0)
  }
  ar <- getComp("ar")
  ne <- getComp("ne")
  end <- getComp("end")

  ## simulate
  simcall <- quote(
    simHHH4lag(ar = ar, ne = ne, end = end, psi = psi, neW = neweights, start = y.start,
               lag.ar = lag.ar, lag.ne = lag.ne, funct_lag = control$funct_lag,
               par_lag = control$par_lag, min_lag = control$min_lag, max_lag = control$max_lag,
               use_distr_lag = control$use_distr_lag)
  )
  if (!simplify) {
    ## result template
    res0 <- object$stsObj[subset,]
    setObserved <- function (observed) {
      res0@observed[] <- observed
      res0
    }
    simcall <- call("setObserved", simcall)
  }
  res <- if (nsim==1) eval(simcall) else
    replicate(nsim, eval(simcall),
              simplify=if (simplify) "array" else FALSE)
  if (simplify) {
    dimnames(res)[1:2] <- list(subset, colnames(model$response))
    attr(res, "initial") <- y.start
    attr(res, "stsObserved") <- object$stsObj[subset,]
    class(res) <- "hhh4sims"
  }

  ## Done
  attr(res, "call") <- cl
  attr(res, "seed") <- RNGstate
  res
}


### Internal auxiliary function, which performs the actual simulation

simHHH4lag <- function(ar,     # lambda_it (nTime x nUnits matrix)
                       ne,     # phi_it (nTime x nUnits matrix)
                       end,    # nu_it (nTime x nUnits matrix, offset included)
                       psi,    # overdisp param(s) or numeric(0) (psi->0 = Poisson)
                       neW,    # weight matrix/array for neighbourhood component
                       start,  # starting counts (vector of length nUnits, or
                       # matrix with nUnits columns if lag > 1)
                       lag.ar, lag.ne, funct_lag, par_lag, min_lag, max_lag, use_distr_lag #BJ
)
{
  nTime <- nrow(end)
  nUnits <- ncol(end)

  ## simulate from Poisson or NegBin model
  rdistr <- if (length(psi)==0 ||
                isTRUE(all.equal(psi, 0, check.attributes=FALSE))) {
    rpois
  } else {
    psi.inv <- 1/psi   # since R uses different parametrization
    ## draw 'n' samples from NegBin with mean vector 'mean' (length=nUnits)
    ## and overdispersion psi such that Variance = mean + psi*mean^2
    ## where 'size'=1/psi and length(psi) == 1 or length(mean)
    function(n, mean) rnbinom(n, mu = mean, size = psi.inv)
  }

  ## if only endemic component -> simulate independently
  if (all(ar + ne == 0)) {
    return(matrix(rdistr(length(end), end), nTime, nUnits))
  }

  ## weighted sum of counts of other (neighbouring) regions
  ## params: y - vector with (lagged) counts of regions #BJ: distributed lags act here.
  ##         W - nUnits x nUnits adjacency/weight matrix (0=no neighbour)
  wSumNE <- if (is.null(neW) || all(neW == 0)) { # includes the case nUnits==1
    function (y, W) numeric(nUnits)
  } else function (y, W) .colSums(W * y, nUnits, nUnits)

  ## initialize matrices for means mu_i,t and simulated data y_i,t
  mu <- y <- matrix(0, nTime, nUnits)
  y <- rbind(start, y)
  nStart <- nrow(y) - nrow(mu)        # usually just 1 for lag=1

  ## simulate
  timeDependentWeights <- length(dim(neW)) == 3
  if (!timeDependentWeights) neWt <- neW
  for(t in seq_len(nTime)){
    if (timeDependentWeights) neWt <- neW[,,t]
    ## mean mu_i,t = lambda*y_i,t-1 + phi*sum_j wji*y_j,t-1 + nu_i,t
    Ylagged <- hhh4addon:::weightedSumAR(observed = y[nStart + t - (max_lag:0), , drop = FALSE], lag = lag.ar, #BJ
                                         funct_lag = funct_lag, par_lag = par_lag, min_lag = min_lag, max_lag = max_lag, #BJ
                                         use_distr_lag = use_distr_lag, sum_up = TRUE)[max_lag + 1, ] #BJ

    if(!is.null(neW)){
      Ylagged.ne <- hhh4addon:::weightedSumNE(y[nStart + t - (max_lag:0), , drop = FALSE], weights = neW, lag = lag.ne, #BJ
                                              funct_lag = funct_lag, #BJ
                                              par_lag = par_lag, #BJ
                                              min_lag = min_lag,
                                              max_lag = max_lag, #BJ
                                              use_distr_lag = use_distr_lag, #BJ
                                              sum_up = TRUE)[max_lag + 1, ]
    }else{
      Ylagged.ne <- 0
    }

    mu[t,] <-
      ar[t,] * Ylagged +
      ne[t,] * Ylagged.ne +
      end[t,]
    ## Sample from Poisson/NegBin with that mean
    y[nStart+t,] <- rdistr(nUnits, mu[t,])
  }

  ## return simulated data without initial counts
  y[-seq_len(nStart),,drop=FALSE]
}

# originally in hh4_plot #BJ
## extract exppreds multiplied with offsets
## note: theta = coef(object) would also work since psi is not involved here
get_exppreds_with_offsets_lag <- function (object,
                                           subset = seq_len(nrow(object$stsObj)),
                                           theta = object$coefficients)
{
  model <- terms.hhh4lag(object)
  means <- surveillance:::meanHHH(theta, model, subset = subset)
  res <- sapply(X = c("ar", "ne", "end"), FUN = function (comp) {
    exppred <- means[[paste0(comp, ".exppred")]]
    offset <- object$control[[comp]]$offset
    if (length(offset) > 1) offset <- offset[subset,,drop=FALSE]
    exppred * offset
  }, simplify = FALSE, USE.NAMES = TRUE)
  res
}
