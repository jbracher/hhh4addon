# Wrapper around hhh4_lag to fit par_lag
fit_par_lag <- function(stsObj, control, check.analyticals = FALSE, range_par){
  control$ar$use_distr_lag <- control$ne$use_distr_lag <- TRUE
  AICs <- rep(NA, length(range_par))
  best_mod <- hhh4_lag(stsObj, control, check.analyticals)
  for(i in 1:length(range_par)){
    control$ar$par_lag <- control$ne$par_lag <- range_par[i]
    mod_temp <- hhh4_lag(stsObj, control, check.analyticals)
    AICs[i] <- AIC(mod_temp) + 1 # +1 because of additional parameter
    if(AICs[i] < min(AICs[1:(i - 1)])){
      best_mod <- mod_temp
    }
  }
  return(list(best_mod = best_mod, AICs = AICs))
}


#' Fitting hhh4 models with distributed lags
#'
#' A modified version of \code{surveillance::hhh4} to allow for distributed lags.
#'
#' @param stsObj,control,check.analyticals As in \code{surveillance::hhh4}, but \code{control}
#' allows for some additional arguments (see details)
#'
#' In this modified version of \code{surveillance::hhh4}, distributed laga can be specified using the following
#' additional elements in the \code{ar} and \code{ne} parts of \code{control}:
#' \itemize{
#'   \item{\code{funct_lag}}{ Function to calculate weights for different lags; defaults to geometric lags (\code{hhh4addon:::geometric_lag}).
#'   Function needs to take the following parameters:
#'   \itemize{
#'   \item{lag1}{ Matrix containing the time-varying \eqn{\phi_{ijt}} which enter into \eqn{\phi_{ijqt} = \phi_{ijt}*w_q}}
#'   \item{par_lag}{ A scalar parameter or a list of parameters to steer \eqn{w_q}}
#'   }}
#'   \item{\code{max_lag}}{ Maximum number of lags after which}
#' }
#'
#' @export
hhh4_lag <- function (stsObj, control = list(
  ar = list(f = ~ -1,        # a formula "exp(x'lamba)*y_t-lag" (ToDo: matrix)
            offset = 1,      # multiplicative offset
            lag = 1,         # autoregression on y_i,t-lag; BJ: changed default to NA
            funct_lag = geometric_lag, par_lag = list(mu = 1), max_lag = 5, use_distr_lag = FALSE), #BJ: added arguments
  ne = list(f = ~ -1,        # a formula "exp(x'phi) * sum_j w_ji * y_j,t-lag"
            offset = 1,      # multiplicative offset
            lag = 1,         # regression on y_j,t-lag; BJ: changed default to NA
            funct_lag = geometric_lag, par_lag = list(mu = 1), max_lag = 5, use_distr_lag = FALSE, #BJ: added arguments
            weights = neighbourhood(stsObj) == 1,  # weights w_ji
            scale = NULL,    # such that w_ji = scale * weights
            normalize = FALSE), # w_ji -> w_ji / rowSums(w_ji), after scaling
  end = list(f = ~ 1,        # a formula "exp(x'nu) * n_it"
             offset = 1),    # optional multiplicative offset e_it
  family = c("Poisson", "NegBin1", "NegBinM"), # or a factor of length nUnit
  subset = 2:nrow(stsObj),   # epidemic components require Y_{t-lag}
  optimizer = list(stop = list(tol = 1e-5, niter = 100), # control arguments
                   regression = list(method = "nlminb"), # for optimization
                   variance = list(method = "nlminb")),  # <- or "Nelder-Mead"
  verbose = FALSE,           # level of reporting during optimization
  start = list(fixed = NULL, # list of start values, replacing initial
               random = NULL,  # values from fe() and ri() in 'f'ormulae
               sd.corr = NULL),
  data = list(t = stsObj@epoch - min(stsObj@epoch)), # named list of covariates
  keep.terms = FALSE  # whether to keep interpretControl(control, stsObj)
), check.analyticals = FALSE)
{
  ptm <- proc.time()

  ## Convert old disProg class to new sts class
  if (inherits(stsObj, "disProg")) {
    stsObj <- surveillance:::disProg2sts(stsObj)
  } else {
    stopifnot(inherits(stsObj, "sts"))
  }

  ## check control and set default values (for missing arguments)
  control <- hhh4addon:::setControl(control, stsObj) #BJ using the updated function

  ## get model terms
  model <- hhh4addon:::interpretControl(control, stsObj) #BJ using the updated function
  dimFixedEffects <- model$nFE + model$nd + model$nOverdisp
  dimRandomEffects <- model$nRE

  ## starting values
  #* -> better default values possible
  theta.start <- model$initialTheta
  Sigma.start <- model$initialSigma

  ## check if initial values are valid
  ## CAVE: there might be NA's in mu if there are missing values in Y
  mu <- surveillance:::meanHHH(theta.start, model, total.only=TRUE)
  if(any(mu==0, na.rm=TRUE) || any(is.infinite(mu)))
    stop("some mean is degenerate (0 or Inf) at initial values")

  ## check score vector and fisher information at starting values
  check.analyticals <- if (isTRUE(check.analyticals)) {
    if (length(theta.start) > 50) "maxLik" else "numDeriv"
  } else if (is.character(check.analyticals)) {
    match.arg(check.analyticals, c("numDeriv", "maxLik"), several.ok=TRUE)
  } else NULL
  if (length(check.analyticals) > 0L) {
    resCheck <- surveillance:::checkAnalyticals(model, theta.start, Sigma.start,
                                                methods=check.analyticals)
    return(resCheck)
  }

  ## maximize loglikelihood (penalized and marginal)
  myoptim <- surveillance:::fitHHH(theta=theta.start,sd.corr=Sigma.start, model=model,
                                   cntrl.stop       = control$optimizer$stop,
                                   cntrl.regression = control$optimizer$regression,
                                   cntrl.variance   = control$optimizer$variance,
                                   verbose=control$verbose)

  ## extract parameter estimates
  convergence <- myoptim$convergence == 0
  thetahat <- myoptim$theta
  if (dimRandomEffects>0) {
    Sigma.orig <- myoptim$sd.corr
    Sigma.trans <- surveillance:::getSigmai(head(Sigma.orig,model$nVar),
                                            tail(Sigma.orig,model$nCorr),
                                            model$nVar)
    dimnames(Sigma.trans) <-
      rep.int(list(sub("^sd\\.", "",
                       names(Sigma.orig)[seq_len(model$nVar)])), 2L)
  } else {
    Sigma.orig <- Sigma.trans <- NULL
  }

  ## compute covariance matrices of regression and variance parameters
  cov <- try(solve(myoptim$fisher), silent=TRUE)
  Sigma.cov <- if(dimRandomEffects>0) try(solve(myoptim$fisherVar), silent=TRUE)

  ## check for degenerate fisher info
  if(inherits(cov, "try-error")){ # fisher info is singular
    if (control$verbose)
      cat("WARNING: Final Fisher information matrix is singular!\n")
    convergence <- FALSE
  } else if(any(!is.finite(diag(cov))) || any(diag(cov)<0)){
    if (control$verbose)
      cat("WARNING: non-finite or negative covariance of regression parameters!\n")
    convergence <- FALSE
  }
  if (!convergence) {
    if (control$verbose) {
      cat("Penalized loglikelihood =", myoptim$loglik, "\n")
      thetastring <- paste(round(thetahat,2), collapse=", ")
      thetastring <- strwrap(thetastring, exdent=10, prefix="\n", initial="")
      cat("theta = (", thetastring, ")\n")
    }
    warning("Results are not reliable!", surveillance:::ADVICEONERROR)
  }

  #BJ: calculate distributed lags:
  distr_lag_ar <- if(control$ar$use_distr_lag){
    m1 <- matrix(1, nrow = control$ar$max_lag)
    control$ar$funct_lag(m1, par_lag = control$ar$par_lag, max_lag = control$ar$max_lag, sum_up = FALSE)[1,,]
  }else{ NA}
  distr_lag_ne <- if(control$ne$use_distr_lag){
    m1 <- matrix(1, nrow = control$ne$max_lag)
    control$ar$funct_lag(m1, par_lag = control$ar$par_lag, max_lag = control$ar$max_lag, sum_up = FALSE)[1,,]
  }else{NA}

  ## gather results in a list -> "hhh4" object
  result <- list(coefficients=thetahat,
                 se=if (convergence) sqrt(diag(cov)), cov=cov,
                 Sigma=Sigma.trans,     # estimated covariance matrix of ri's
                 Sigma.orig=Sigma.orig, # variance parameters on original scale
                 Sigma.cov=Sigma.cov,   # covariance matrix of Sigma.orig
                 call=match.call(),
                 dim=c(fixed=dimFixedEffects,random=dimRandomEffects),
                 loglikelihood=myoptim$loglik, margll=myoptim$margll,
                 convergence=convergence,
                 fitted.values=surveillance:::meanHHH(thetahat, model, total.only=TRUE),
                 control=control,
                 terms=if(control$keep.terms) model else NULL,
                 stsObj=stsObj,
                 lags=sapply(control[c("ar","ne")], function (comp)
                   if (comp$inModel) comp$lag else NA_integer_),
                 func_lag = list(ar = control$ar$funct_lag, ne = control$ne$funct_lag), #BJ
                 par_lag = list(ar = control$ar$par_lag, ne = control$ne$par_lag), #BJ
                 max_lags = sapply(control[c("ar", "ne")], #BJ
                                   function(comp) if (comp$inModel) comp$max_lag else NA_integer_), #BJ
                 use_distr_lag = c(ar = control$ar$use_distr_lag, ne = control$ne$use_distr_lag),
                 distr_lag = list(ar = distr_lag_ar, ne = distr_lag_ne),
                 nObs=sum(!model$isNA[control$subset,]),
                 nTime=length(model$subset), nUnit=ncol(stsObj),
                 ## CAVE: nTime is not nrow(stsObj) as usual!
                 runtime=proc.time()-ptm)
  if (!convergence) {
    ## add (singular) Fisher information for further investigation
    result[c("fisher","fisherVar")] <- myoptim[c("fisher","fisherVar")]
  }
  class(result) <- c("hhh4lag", "hhh4")
  return(result)
}
