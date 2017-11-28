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
#' around it and tries a set of possible parameters provided in the argument \code{range_par}.
#'
#'
#' @param stsObj,control,check.analyticals As in \code{surveillance::hhh4}, but \code{control}
#' allows for some additional arguments
#'
#' In this modified version of \code{surveillance::hhh4}, distributed lags can be specified by
#' additional elements in the \code{ar} and \code{ne} parts of \code{control}:
#' \itemize{
#'   \item{\code{funct_lag}}{ Function to calculate the a matrix of distributed lags from a matrix of first lags.
#'    Currently only geometric lags (\code{hhh4addon:::geometric_lag}) are available and set as default, see Details.
#'   The function has to take the following arguments:
#'   \itemize{
#'   \item{\code{lag1}}{ Matrix containing the first lags which would be used in a standard \code{hhh4} model.}
#'   \item{\code{par_lag}}{ A scalar parameter to steer \eqn{u_q}. For the geometric lags this is the un-normalized weight of the first lag.}
#'   \item{\code{max_lag}}{ Maximum number of lags.}
#'   \item{\code{sum_up}}{ Specifies how detailed the output of the function is - only for internal use.}
#'   }}
#'   \item{\code{max_lag}}{ Specification of the \code{max_lag} argument passed to funct_lag} to compute the lags.
#' }
#' Unlike in \code{hhh4_lag} the par_lag argument for \code{funct_lag} is not specified directly
#' by the user; instead the model is re-fit for each parameter value provided in \code{range_par}.
#'
#' @param range_par a vector of values to try for the \code{par_lag} argument of \code{funct_lag}
#' @export
fit_par_lag <- function(stsObj, control, check.analyticals = FALSE, range_par){
  control$ar$use_distr_lag <- control$ne$use_distr_lag <- TRUE
  AICs <- rep(NA, length(range_par))
  best_mod <- mod_temp <- NULL
  for(i in 1:length(range_par)){
    control$ar$par_lag <- control$ne$par_lag <- range_par[i]
    mod_temp <- if(is.null(mod_temp) || mod_temp$convergence == FALSE){
      hhh4_lag(stsObj, control, check.analyticals)
    }else{
      update(mod_temp, ar = control$ar, ne = control$ne)
    }
    mod_temp$dim["fixed"] <- mod_temp$dim["fixed"] + 1 # + 1 for decay paramter
    AICs[i] <- AIC(mod_temp) # + 2 no longer necessary as +1 added above
    if(i == 1){ # keep first model in all cases
      best_mod <- mod_temp
    }else{ # otherwise only if it beats the previous ones.
      if(AICs[i] < min(AICs[1:(i - 1)])){
        best_mod <- mod_temp
      }
    }
  }
  if(which.min(AICs) == 1) warning("The minimum AIC is reached with the smallest value of par_lag. Consider lowering this value.")
  if(which.min(AICs) == length(AICs)) warning("The minimum AIC is reached with the largest value of par_lag. Consider increasing this value.")
  return(list(best_mod = best_mod, AICs = AICs))
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
#' The simple first lags are replaced by weighted sums of the Q
#' previous observations. The weights u_q, q = 1, ..., Q sum up to
#' 1 and need to be parametrizable by a single scalar parameter \code{par_lag}.
#' This parameter is passed to a function \code{funct_lag} which takes
#' the first lags and transforms them into distributed lags. Currently
#' only geometric lags (function \code{geometric_lag}) are available.
#' These are specified as u0_q = p^q * (1 - p)^{q - 1} and u_q = u0_q /
#' sum_{q = 1}^Q u0_q. The \code{par_lag} parameter corresponds to u0_1,
#' i.e. the un-normalized weight of the first lag.
#'
#' @param stsObj,control,check.analyticals As in \code{surveillance::hhh4},
#' but the \code{control} argument allows for some additional specifications.
#'
#' Distributed lags can be specified by
#' additional elements in the \code{ar} and \code{ne} parts of \code{control}:
#' \itemize{
#'   \item{\code{use_distr_lag}}{ Logical: should distributed lags be used instead of ordinary lags as
#'   implemented in \code{surveillance}?}
#'   \item{\code{funct_lag}}{ Function to calculate the a matrix of
#'   distributed lags from a matrix of first lags. Currently only
#'   geometric lags (\code{hhh4addon:::geometric_lag}) are available
#'   and set as default, see Details. The function has to take the
#'   following arguments:
#'   \itemize{
#'   \item{\code{lag1}}{ Matrix containing the first lags which would
#'   be used in a standard \code{hhh4} model.}
#'   \item{\code{par_lag}}{ A scalar parameter to steer \eqn{u_q}. For
#'   the geometric lags this is the un-normalized weight of the first lag.}
#'   \item{\code{max_lag}}{ Maximum number of lags.}
#'   \item{\code{sum_up}}{ Specifies how detailed the output of the
#'   function is - only for internal use.}
#'   }}
#'   \item{\code{par_lag, max_lag}}{ Specification of the arguments
#'   passed to funct_lag} to compute the distributed  lags.
#' }
#' The current implementation requires the lag structure to be handled
#' the same way in the «code{ar} and the \code{ne} components.
#' The parameter \code{par_lag} can be estimated using a profile
#' likelihood approach. This is done using the wrapper \code{fit_par_lag}.
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
#'                                      par_lag = 0.8, use_distr_lag = TRUE),
#'                            family = "NegBinM", subset = 6:312)
#' fit_salmonella <- hhh4_lag(salmonella, control_salmonella)
#' summary(fit_salmonella)
#'
#' @export
hhh4_lag <- function (stsObj, control = list(
  ar = list(f = ~ -1,        # a formula "exp(x'lamba)*y_t-lag" (ToDo: matrix)
            offset = 1,      # multiplicative offset
            lag = 1,         # autoregression on y_i,t-lag; BJ: changed default to NA
            funct_lag = geometric_lag, par_lag = 1, max_lag = 5, use_distr_lag = FALSE), #BJ: added arguments
  ne = list(f = ~ -1,        # a formula "exp(x'phi) * sum_j w_ji * y_j,t-lag"
            offset = 1,      # multiplicative offset
            lag = 1,         # regression on y_j,t-lag; BJ: changed default to NA
            funct_lag = geometric_lag, par_lag = 1, max_lag = 5, use_distr_lag = FALSE, #BJ: added arguments
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
