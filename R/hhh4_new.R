################################################################################
### The following are modified versions of functions from the surveillance package
### and wrappers around them.
### See below the original copyright declaration.
################################################################################

################################################################################
### Copyright declaration from the surveillance package:
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### hhh4 is an extended version of algo.hhh for the sts-class
### The function allows the incorporation of random effects and covariates.
###
### Copyright (C) 2010-2012 Michaela Paul, 2012-2016 Sebastian Meyer
### $Revision: 1706 $
### $Date: 2016-05-03 16:09:49 +0200 (Die, 03. Mai 2016) $
################################################################################

#' Estimating the lag decay parameter of an \code{hhh4_lag} model using profile likelihood
#'
#' Wrapper around \code{hhh4_lag} to allow for profile likelihood estimation of the scalar parameter
#' governing the lag structure. \code{hhh4_lag} can fit models with fixed lag decay parameter; \code{fit_par_lag} loops
#' around it and tries a set of possible parameters provided in the argument \code{range_par}. NOTE: this will
#' soon be replaced by \code{profile_par_lag} which does the same, but using \code{optim..., method = "Brent", ...)}.
#'
#'
#' @param{stsObj,control,check.analyticals As in \code{surveillance::hhh4}, but \code{control}
#' allows for some additional elements in order to specify a distributed lag structure:
#' \itemize{
#'   \item{\code{funct_lag}}{ Function to compute the lag weights \eqn{u_q} (see details) depending on a scalar
#'   parameter \code{par_lag}. The function has to take the
#'   following arguments:
#'   \itemize{
#'   \item{\code{par_lag}}{ A scalar parameter to steer \eqn{u_q}. It should be specified in a way which allows it to
#'   take any value in the real numbers}
#'   \item{\code{min_lag,max_lag}}{ Minimum and maximum lags; e.g. \code{min_lag = 3, max_lag = 6} will assign all weights to lags 3 through 6.
#'   Usually \code{min_lag} is set to 1, higher values can be useful for direct forecasting at higher horizons.}
#'   }}
#'   \item{\code{min_lag, max_lag}}{ Specification of the arguments passed to funct_lag} to compute the distributed lags. Unlike in
#'   \code{hhh4_lag}, \code{par_lag} is not to be specified as it is estimated from the data.
#' }}
#'
#' In this modified version of \code{surveillance::hhh4}, distributed lags can be specified by
#' additional elements \code{control} argument:
#' \itemize{
#'   \item{\code{funct_lag}}{ Function to compute the lag weights \eqn{u_q} (see details) depending on a scalar
#'   parameter \code{par_lag}. The function has to take the
#'   following arguments:
#'   \itemize{
#'   \item{\code{par_lag}}{ A scalar parameter to steer \eqn{u_q}. It should be specified in a way which allows it to
#'   take any value in the real numbers}
#'   \item{\code{min_lag,max_lag}}{ Minimum and maximum lags; e.g. \code{min_lag = 3, max_lag = 6} will assign all weights to lags 3 through 6.
#'   Usually \code{min_lag} is set to 1, higher values can be useful for direct forecasting at higher horizons.}
#'   }}
#'   \item{\code{min_lag, max_lag}}{ Specification of the arguments passed to funct_lag} to compute the distributed lags. Unlike in
#'   \code{hhh4_lag}, \code{par_lag} is not to be specified as it is estimated from the data.
#'   }}
#'   \item{\code{max_lag}}{ Specification of the \code{max_lag} argument passed to funct_lag} to compute the lags.
#' }
#' Unlike in \code{hhh4_lag} the par_lag argument for \code{funct_lag} is not specified directly
#' by the user; instead the model is re-fit for each parameter value provided in \code{range_par}.
#'
#' @param range_par a vector of values to try for the \code{par_lag} argument of \code{funct_lag}
#' @param use_update should results from previous values in range_par be used as starting value for next iteration (via \code{update})?
#' @export
fit_par_lag <- function(stsObj, control, check.analyticals = FALSE, range_par, use_update = TRUE){
  AICs <- rep(NA, length(range_par))
  best_mod <- mod_temp <- NULL
  for(i in 1:length(range_par)){
    control$par_lag <- range_par[i]
    mod_temp <- if(is.null(mod_temp) || mod_temp$convergence == FALSE || use_update == FALSE){
      hhh4_lag(stsObj, control, check.analyticals)
    }else{
      update(mod_temp, par_lag = control$par_lag)
    }
    if(mod_temp$convergence == FALSE & use_update){ # catch convergence errors by trying to fit without update
      warning("Model with par_lag = ", range_par[i], " did not converge using update(). Refitting from scratch....")
      mod_temp <- hhh4_lag(stsObj, control, check.analyticals)
    }
    mod_temp$dim["fixed"] <- mod_temp$dim["fixed"] + 1 # + 1 for decay paramter
    AICs[i] <- AIC(mod_temp) # + 2 no longer necessary as +1 added above
    if(i == 1){ # keep first model in all cases
      best_mod <- mod_temp
    }else{ # otherwise only if it beats the previous ones.
      if(mod_temp$convergence){
        if(AICs[i] < min(AICs[1:(i - 1)], na.rm = TRUE)){
          best_mod <- mod_temp
        }
      }else{
        warning("Model with par_lag = ", range_par[i], " did not converge. AIC set to NA.")
      }
    }
  }
  if(which.min(AICs) == 1) warning("The minimum AIC is reached with the smallest value of par_lag. Consider lowering this value.")
  if(which.min(AICs) == length(AICs)) warning("The minimum AIC is reached with the largest value of par_lag. Consider increasing this value.")
  return(list(best_mod = best_mod, AICs = AICs))
}

#' Estimating the lag decay parameter of an \code{hhh4_lag} model using profile likelihood
#'
#' Wrapper around \code{hhh4_lag} to allow for profile likelihood estimation of the scalar parameter
#' governing the lag structure. \code{hhh4_lag} can fit models with fixed lag decay parameter; \code{fit_par_lag} loops
#' around it finds the optimal value of \code{par_lag}. See \code{?hhh4_lag} for details.
#'
#'
#' @param stsObj,control,check.analyticals As in \code{surveillance::hhh4}, but \code{control}
#' allows for some additional arguments in order to specify a distributed lag structure:
#' \itemize{
#'   \item{\code{funct_lag}}{ Function to compute the lag weights \eqn{u_q} (see details) depending on a scalar
#'   parameter \code{par_lag}. The function has to take the
#'   following arguments:
#'   \itemize{
#'   \item{\code{par_lag}}{ A scalar parameter to steer \eqn{u_q}. It should be specified in a way which allows it to
#'   take any value in the real numbers}
#'   \item{\code{min_lag,max_lag}}{ Minimum and maximum lags; e.g. \code{min_lag = 3, max_lag = 6} will assign all weights to lags 3 through 6.
#'   Usually \code{min_lag} is set to 1, higher values can be useful for direct forecasting at higher horizons.}
#'   }}
#'   \item{\code{min_lag, max_lag}}{ Specification of the arguments passed to funct_lag} to compute the distributed lags. Unlike in
#'   \code{hhh4_lag}, \code{par_lag} is not to be specified as it is estimated from the data.
#' }
#' @param start_par_lag A starting value for \code{par_lag}
#' @param lower_par_lag,upper_par_lag lower and upper limits for the value of par_lag; defaults to -10, 10
#' @param return_full_cov logical: should the full covariance matrix of the parameter estimates (including \code{par_lag})
#' be obtained numerically?
#' @param reltol_par_lag the relative tolerance passed to the \code{optim} call to identify \code{par_lag}
#'
#' @return If \code{return_full_cov == FALSE}: an \code{hhh4_lag} object. If \code{return_full_cov == TRUE} A list with two
#' elements: \code{best_mod} is the \code{hhh4_lag} fit for the best value of \code{par_lag}; \code{cov} is an extended covariance matrix for the regression parameters
#' which also includes par_lag.
#' @export
profile_par_lag <- function(stsObj, control, start_par_lag = 0.5, lower_par_lag = -10, upper_par_lag = 10,
                            return_full_cov = FALSE, reltol_par_lag = 1e-08, check.analyticals = FALSE){
  control$par_lag <- start_par_lag
  initial_fit <- hhh4_lag(stsObj = stsObj, control = control)
  profile_lik <- function(par_lag){
    # par_lag <- exp(logit_par_lag)/(1 + exp(logit_par_lag))
    control$par_lag <- par_lag
    # control$start <- initial_fit$coefficients
    mod_temp <- hhh4_lag(stsObj, control, check.analyticals)
    return(-mod_temp$loglikelihood)
  }
  opt_par_lag <- optim(par = 0.5, profile_lik, method = "Brent", lower = lower_par_lag, upper = upper_par_lag,
                       control = list(reltol = reltol_par_lag))$par
  control$par_lag <- opt_par_lag
  best_mod <- hhh4_lag(stsObj = stsObj, control = control)
  best_mod$dim["fixed"] <- best_mod$dim["fixed"] + 1 # + 1 for decay paramter
  if(return_full_cov){
    cov = numeric_fisher_hhh4lag(best_mod)
    best_mod$se_par_lag <- cov["par_lag", "par_lag"]
    return(list(best_mod = best_mod, cov = cov))
  }else{
    best_mod$se_par_lag <- NA
    return(best_mod)
  }
}

#' Numerical evaluation of the covariance matrix including the additional parameter
#' \code{par_lg}
#' @param best_mod an \code{hhh4lag} object; should be generated in \code{profile_par_lag} so that
#' the \code{par_lag} parameter is already optimized.
numeric_fisher_hhh4lag <- function(best_mod){
  stsObj <- best_mod$stsObj
  control <- best_mod$control
  coefficients <- c(best_mod$coefficients, par_lag = best_mod$par_lag)
  lik_vect <- function(coefficients){
    control$par_lag <- control$par_lag <- coefficients["par_lag"]
    mod <- hhh4addon:::interpretControl(control, stsObj = stsObj)
    surveillance:::penLogLik(coefficients[names(coefficients) != "par_lag"], NULL, model = mod)
  }
  hess <- numDeriv::hessian(lik_vect, coefficients)
  cov <- -solve(hess)
  colnames(cov) <- rownames(cov) <- names(coefficients)
  return(cov)
}


#' Fitting hhh4 models with distributed lags
#'
#' A modified version of \code{surveillance::hhh4} to allow for distributed
#' lags. Usually used from inside of the wrapper \code{fit_par_lag}.
#'
#' The standard \code{hhh4} function only allows for models with
#' first lags i.e. of the form
#' \deqn{mu_{it} = \lambda_{it}X_{i, t - 1} + \phi_{it}\sum_{j != i}w_{ji}X_{j, t - 1} + \nu_{it},}
#' see \code{?hhh4}. The extension \code{hhh4_lag} allows to specify
#' models of the form
#' \deqn{mu_{it} = \lambda_{it}\sum_{q= 1}^Q u_q X_{i, t - q} + \phi_{it}\sum_{j\neq i}sum_{q= 1}^Q w_{ji}u_q X_{j, t - q} + \nu_{it}.}
#' Here the first lags are now replaced by weighted sums of the Q
#' previous observations. The weights u_q, q = 1, ..., Q sum up to
#' 1 and need to be parametrizable by a single scalar parameter. The value of this parameter needs to be passed as \code{control$par_lag}.
#' Moreover, a function to obtain a vector of weights from \code{par_lag} needs to be provided in \code{control$funct_lag}.
#' Currently three such functions are implemented in the package:
#' \itemize{
#' \item{Geometric lags (function \code{geometric_lag}; the default).
#' These are specified as
#' \deqn{u0_q = \alpha * (1 - \alpha)^{q - 1}}
#' and \eqn{u_q = u0_q / sum_{q = 1}^Q u0_q}  for \eqn{q = 1, ..., Q}. The \code{par_lag} parameter corresponds to logit(\eqn{\alpha}),
#' i.e. the un-normalized weight of the first lag.}
#' \item{Poisson lags (function \code{poisson_lag}).
#' These are specified as
#' \deqn{u0_q =  \alpha^(q - 1)\exp(-\alpha)/(q - 1)!,}
#' and \eqn{u_q = u0_q / sum_{q = 1}^Q u0_q} for \eqn{q = 1, ..., Q}. Note that he Poisson distribution is shifted by one to
#' achieve a positive support. The \code{par_lag} parameter corresponds to log(\eqn{\alpha}).}
#' \item{A weighting only between first and second lags, i.e. \deqn{u_1 = \alpha, u_2 = 1 - \alpha}.
#' The \code{par_lag} parameter corresponds to logit(\eqn{\alpha}).}
#' }
#' Users can specify their own weighting functions as long as they take the arguments described above and return a vector of weights.
#'
#' @param stsObj,control,check.analyticals As in \code{surveillance::hhh4},
#' but with the following additional elements in the \code{control} argument in order to specify a distributed lag structure:
#' \itemize{
#'   \item{\code{funct_lag}}{ Function to compute the lag weights \eqn{u_q} (see details) depending on a scalar
#'   parameter \code{par_lag}. The function has to take the
#'   following arguments:
#'   \itemize{
#'   \item{\code{par_lag}}{ A scalar parameter to steer \eqn{u_q}. It should be specified in a way which allows it to
#'   take any value in the real numbers}
#'   \item{\code{min_lag,max_lag}}{ Minimum and maximum lags; e.g. \code{min_lag = 3, max_lag = 6} will assign all weights to lags 3 through 6.
#'   Usually \code{min_lag} is set to 1, higher values can be useful for direct forecasting at higher horizons.}
#'   }}
#'   \item{\code{par_lag, min_lag, max_lag}}{ Specification of the arguments
#'   passed to funct_lag} to compute the distributed  lags.
#' }
#' \code{hhh4_lag} requires \code{par_lag} to be pre-specified. Using the wrapper \code{fit_par_lag} it can also be estimated using a profile
#' likelihood approach.
#'
#' @examples
#' ## a simple univariate example:
#' data("salmonella.agona")
#' ## convert old "disProg" to new "sts" data class
#' salmonella <- disProg2sts(salmonella.agona)
#' # specify and fit model: fixed geometric lag structure
#' # with weight 0.8 for first lag
#' control_salmonella <- list(end = list(f = addSeason2formula(~ 1)),
#'                            ar = list(f = addSeason2formula(~ 1),
#'                                      par_lag = 0.8),
#'                            family = "NegBinM", subset = 6:312)
#' fit_salmonella <- hhh4_lag(salmonella, control_salmonella)
#' summary(fit_salmonella)
#'
#' @export
hhh4_lag <- function (stsObj, control = list(
  ar = list(f = ~ -1,        # a formula "exp(x'lamba)*y_t-lag" (ToDo: matrix)
            offset = 1,      # multiplicative offset
            lag = NA),         # autoregression on y_i,t-lag; BJ: changed default to NA
  ne = list(f = ~ -1,        # a formula "exp(x'phi) * sum_j w_ji * y_j,t-lag"
            offset = 1,      # multiplicative offset
            lag = NA,         # regression on y_j,t-lag; BJ: changed default to NA
            weights = neighbourhood(stsObj) == 1,  # weights w_ji
            scale = NULL,    # such that w_ji = scale * weights
            normalize = FALSE), # w_ji -> w_ji / rowSums(w_ji), after scaling
  end = list(f = ~ 1,        # a formula "exp(x'nu) * n_it"
             offset = 1),    # optional multiplicative offset e_it
  family = c("Poisson", "NegBin1", "NegBinM"), # or a factor of length nUnit
  funct_lag = geometric_lag, par_lag = 1, min_lag = 1, max_lag = 5,
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
  distr_lag <- control$funct_lag(par_lag = control$par_lag, min_lag = control$min_lag, max_lag = control$max_lag)

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
                 func_lag = control$funct_lag, #BJ
                 par_lag = control$par_lag, #BJ
                 max_lag = control$max_lag, #BJ
                 min_lag = control$min_lag, #BJ
                 distr_lag = distr_lag,
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
