### Updated version: Simulate-method for hhh4-objects
#' @import stats
#'@export
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

  #BJ: use lags or max_lags depending on setting
  # lags <- c(ar = ifelse(control$ar$use_distr_lag, control$ar$max_lag, control$ar$lag),
  #           ne = ifelse(control$ne$use_distr_lag, control$ne$max_lag, control$ne$lag))
  maxlag <- max(c(control$ar$max_lag, control$ar$lag,
                  control$ne$max_lag, control$ne$lag), na.rm = TRUE) #BJ
  #BJ maxlag <- max(lag.ar, lag.ne)

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

  ## get fitted components nu_it (with offset), phi_it, lambda_it, t in subset
  model <- hhh4addon:::terms.hhh4lag(object)
  means <- surveillance:::meanHHH(theta, model, subset=subset)
  psi <- surveillance:::splitParams(theta,model)$overdisp

  ## weight matrix/array of the ne component
  neweights <- surveillance:::getNEweights(object, surveillance:::coefW(theta))

  ## set predictor to zero if not included ('components' argument)
  stopifnot(length(components) > 0, components %in% c("ar", "ne", "end"))
  getComp <- function (comp) {
    sel <- if (comp == "end") "endemic" else paste(comp, "exppred", sep=".")
    res <- means[[sel]]
    if (!comp %in% components) res[] <- 0
    res
  }
  ar <- getComp("ar")
  ne <- getComp("ne")
  end <- getComp("end")

  ## simulate
  simcall <- quote(
    simHHH4(ar = ar, ne = ne, end = end, psi = psi, neW = neweights, start = y.start,
            lag.ar = lag.ar, funct_lag.ar = control$ar$funct_lag, par_lag.ar = control$ar$par_lag, max_lag.ar = control$ar$max_lag, use_distr_lag.ar = control$ar$use_distr_lag,
            lag.ne = lag.ne, funct_lag.ne = control$ne$funct_lag, par_lag.ne = control$ne$par_lag, max_lag.ne = control$ne$max_lag, use_distr_lag.ne = control$ne$use_distr_lag)
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



### updated version Internal auxiliary function, which performs the actual simulation

# THIS IS NOT READY YET

simHHH4 <- function(ar,     # lambda_it (nTime x nUnits matrix)
                    ne,     # phi_it (nTime x nUnits matrix)
                    end,    # nu_it (nTime x nUnits matrix, offset included)
                    psi,    # overdisp param(s) or numeric(0) (psi->0 = Poisson)
                    neW,    # weight matrix/array for neighbourhood component
                    start,  # starting counts (vector of length nUnits, or
                    # matrix with nUnits columns if lag > 1)
                    lag.ar, funct_lag.ar, par_lag.ar, max_lag.ar, use_distr_lag.ar,
                    lag.ne, funct_lag.ne, par_lag.ne, max_lag.ne, use_distr_lag.ne
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
    Ylagged <- hhh4addon:::weightedSumAR(observed = y[nStart + t - (max_lag.ar:0), , drop = FALSE], lag = lag.ar, #BJ
                                        funct_lag = funct_lag.ar, par_lag = par_lag.ar, max_lag = max_lag.ar, #BJ
                                        use_distr_lag = use_distr_lag.ar, sum_up = TRUE)[max_lag.ne + 1, ] #BJ

    Ylagged.ne <- hhh4addon:::weightedSumNE(y[nStart + t - (max_lag.ne:0), , drop = FALSE], weights = neW, lag = lag.ne, #BJ
                             funct_lag = funct_lag.ne, #BJ
                             par_lag = par_lag.ne, #BJ
                             max_lag = max_lag.ne, #BJ
                             use_distr_lag = use_distr_lag.ne, #BJ
                             sum_up = TRUE)[max_lag.ar + 1, ]
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
