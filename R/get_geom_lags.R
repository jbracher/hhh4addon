#' Display the function and parameters used for distributed lags
#'@export
#'@param hhh4Obj an object of class \code{hhh4}
#'@return A list containing the function and parameters characterizeing the lag weights
#' (for both the \code{ar} and \code{ne} components)
#'@examples
#' data("salmonella.agona")
#' ## convert old "disProg" to new "sts" data class
#' salmonella <- disProg2sts(salmonella.agona)
#' # specify and fit model: fixed geometric lag structure
#' # with weight 0.8 for first lag
#' control_salmonella <- list(end = list(f = addSeason2formula(~ 1)),
#'                            ar = list(f = addSeason2formula(~ 1),
#'                            par_lag = 0.8),
#'                            family = "NegBinM", subset = 6:312)
#' fit_salmonella <- hhh4_lag(salmonella, control_salmonella)
#' distr_lag(fit_salmonella)

distr_lag <- function(hhh4Obj){
  if(!("hhh4lag" %in% class(hhh4Obj))){
    stop("structure of distributed lags can only be extracted from objects of class hhh4lag.")
  }

  ret <- list(funct_lag = hhh4Obj$control$funct_lag,
              par_lag = hhh4Obj$control$par_lag,
              min_lag = hhh4Obj$control$min_lag,
              max_lag = hhh4Obj$control$max_lag)

  return(ret)
}

#' Transform matrix of first-order lagged observations to matrix of weighted sums of past observation
#'
#' This function modifies the design matrices from first-order lags to weighted lags with weighted structure. Used within \code{weightedSumNE}
#' and \code{weightedSumAR}.
#' @param lag1 a matrix of first lags as usually used in \code{hhh4}.
#' @param lag_weigths a vector of weights; the length of this vector determines the number of lags.
#' @param sum_up \code{sum_up = FALSE} returns a more detailed output; for debugging only.
get_weighted_lags <- function(lag1, lag_weights, sum_up = FALSE){

  if(any(is.na(lag_weights)) | any(is.nan(lag_weights))){
    stop(paste("Lag weights are NA or NaN:", paste(lag_weights, collapse = "; ")))
  }

  if(abs(sum(lag_weights) - 1) > 0.0001 | any(lag_weights < 0)){
    stop("Lag weights need to be positive and sum up to 1. Please make sure your lag weighting function only returns valid weights.")
  }
  max_lag <- length(lag_weights)
  n_units <- ncol(lag1)
  weighted_lags <- if (sum_up) {
    matrix(0, ncol = ncol(lag1), nrow = nrow(lag1))
  }
  else {
    array(NA, dim = c(dim(lag1), max_lag))
  }
  for (i in 1:max_lag) {
    lag_i <- rbind(matrix(NA, nrow = i - 1, ncol = n_units),
                   lag1[c(1:(nrow(weighted_lags) - i + 1)), , drop = FALSE])
    if (sum_up) {
      weighted_lags <- weighted_lags + lag_weights[i] * lag_i
    }
    else {
      weighted_lags[, , i] <- lag_weights[i] * lag_i
    }
  }
  weighted_lags
}


#' Function to obtain geometric weights

#' This function generates geometric weights which are subsequently used inside of \code{get_weighted_lags}. To be passed
#' to \code{hhh4_lag} or \code{profile_par_lag} as the \code{control$funct_lag} argument.
#' @param par_lag a parameter to steer the lag structure, here \eqn{logit(p)} where \eqn{p} is the parameter of
#' the geometric distribution characterizing the lag structure; see details of \code{hhh4lag} or \code{profile_par_lag}.
#' @param min_lag smallest lag to include; the support of the Poisson form starts only at \code{min_lag}. Defaults to 1.
#' @param max_lag highest lag to include; higher lags are cut off and he remaining weights standardized. Defaults to 5.
#' @export
geometric_lag <- function(par_lag, min_lag, max_lag){
  p_lag <- exp(par_lag)/(1 + exp(par_lag))
  weights0 <- c(rep(0, min_lag - 1), dgeom((min_lag:max_lag) -
                                             1, p_lag))
  weights <- weights0/sum(weights0)
  return(weights)
}

#' Function to obtain Poisson weights

#' This function generates Poisson weights which are subsequently used inside of \code{get_weighted_lags}. To be passed
#' to \code{hhh4_lag} or \code{profile_par_lag} as the \code{control$funct_lag} argument.
#' @param par_lag a parameter to steer the lag structure, here \eqn{log(\mu)} where \eqn{\mu} is the parameter of
#' the Poisson distribution characterizing the lag structure; see details of \code{hhh4lag} or \code{profile_par_lag}.
#' @param min_lag smallest lag to include; the support of the Poisson form starts only at \code{min_lag}. Defaults to 1.
#' @param max_lag highest lag to include; higher lags are cut off and he remaining weights standardized. Defaults to 5.
#' @export
poisson_lag <- function(par_lag, min_lag, max_lag){
  mu_lag <- exp(par_lag)
  weights0 <- c(rep(0, min_lag - 1), dpois((min_lag:max_lag) -
                                             1, mu_lag) + 0.0001) # avoids numerical problems
  weights <- weights0/sum(weights0)
  return(weights)
}

#' Function to obtain AR2 weights

#' This function generates AR2 weights which are subsequently used inside of \code{get_weighted_lags}. To be passed
#' to \code{hhh4_lag} or \code{profile_par_lag} as the \code{control$funct_lag} argument.
#' @param par_lag a parameter to steer the lag structure, here \eqn{logit(p)} where \eqn{p} is the weight of
#' the first lag; see details of \code{hhh4lag} or \code{profile_par_lag}.
#' @param min_lag smallest lag to include; the support of the Poisson form starts only at \code{min_lag}. Defaults to 1.
#' @param max_lag highest lag to include; higher lags are cut off and he remaining weights standardized. Defaults to 5.
#' @export
ar2_lag <- function(par_lag, min_lag, max_lag){
  p_lag <- exp(par_lag)/(1 + exp(par_lag))
  weights0 <- c(rep(0, min_lag - 1), p_lag, 1 - p_lag, rep(0, max_lag - min_lag - 1))
  weights <- weights0/sum(weights0)
  return(weights)
}


#' Function to obtain linearly decaying weights
#'
#' This function generates linearly decaying weights which are subsequently used inside of \code{get_weighted_lags}. To be passed
#' to \code{hhh4_lag} or \code{profile_par_lag} as the \code{control$funct_lag} argument. There are different ways of
#' specifying linearly decaying weights, the version implemented here is
#' \deqn{u0_q = max(1 - mq, 0)}
#' and \eqn{u_q = u0_q / sum_{q = 1}^Q u0_q} for \eqn{q = 1, ..., Q}.
#' @param par_lag a parameter to steer the lag structure, here \eqn{logit(m)} where \eqn{m} is the linear decay
#' factor; see details of \code{hhh4lag} or \code{profile_par_lag}.
#' @param min_lag smallest lag to include; the support of the Poisson form starts only at \code{min_lag}. Defaults to 1.
#' @param max_lag highest lag to include; higher lags are cut off and he remaining weights standardized. Defaults to 5.
#' @export
linear_lag <- function (par_lag, min_lag, max_lag){
  m_lag <- exp(par_lag)/(1 + exp(par_lag))
  weights0 <- pmax(c(rep(0, min_lag - 1),
                     1 - m_lag*(0:(max_lag - min_lag))), 0)
  weights <- weights0/sum(weights0)
  return(weights)
}

#' Function to obtain unrestricted lags
#'
#' This function generates \eqn{Q} lag weights without a parametric constraint. The weights are obtained via
#' a multinomial logit transformation where the first lag is the reference category.
#' @param par_lag theparameter vector of length \eqn{Q - 1}
#' @param min_lag smallest lag to include; the support of the Poisson form starts only at \code{min_lag}. Defaults to 1.
#' @param max_lag highest lag to include; higher lags are cut off and he remaining weights standardized. Defaults to 5.
#' @export
unrestricted_lag <- function(par_lag, min_lag, max_lag){
  if(length(par_lag) != max_lag - min_lag){
    stop("The starting value for par_lag needs to be of length max_lag - min_lag + 1.")
  }
  weights0 <- c(1, exp(par_lag))/sum(c(1, exp(par_lag)))
  weights <- c(rep(0, min_lag - 1), weights0)
  return(weights)
}
