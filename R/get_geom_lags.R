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

  list(funct_lag = hhh4Obj$control$funct_lag,
       par_lag = hhh4Obj$control$par_lag,
       min_lag = hhh4Obj$control$min_lag,
       max_lag = hhh4Obj$control$max_lag)

  return(ret)
}

#' Transform matrix of first-order lagged observations to matrix of weighted sums of past observation (with geometric weights)
#'
#' This function modifies the design matrices from first-order lags to weighted lags with geometric structure. To be passed to \code{hhh4lag} or
#' \code{profile_par_lag} in order to use geometric lags.
#' @param lag1 a matrix of first lags as usually used in \code{hhh4}.
#' @param par_lag a parameter to steer the lag structure, here \eqn{logit(p)} where \eqn{p} is the parameter of
#' the geometric distribution characterizing the lag structure; see details of \code{hhh4lag} or \code{profile_par_lag}
#' @param min_lag smallest lag to include; the support of the geometric form starts only at \code{min_lag}. Defaults to 1.
#' @param max_lag highest lag to include; higher lags are cut off and he remaining weights standardized. Defaults to 5.
#' @param sum_up \code{sum_up = FALSE} returns a more detailed output; for debugging only
#' @export
geometric_lag <- function(lag1, par_lag, min_lag, max_lag, sum_up = TRUE){ #BJ

  p_lag <- exp(par_lag)/(1 + exp(par_lag))

  weights0 <- c(rep(0, min_lag - 1), dgeom((min_lag:max_lag) - 1, p_lag)) #BJ
  weights <- weights0/sum(weights0) #BJ
  geom_lag <- if(sum_up){
    matrix(0, ncol = ncol(lag1), nrow = nrow(lag1))
  }else{
    array(NA, dim = c(dim(lag1), max_lag))
  }

  # first line contains NAs
  for(i in 1:max_lag){ #BJ
    lag_i <- lag1[c(rep(1, i - 1), 1:(nrow(geom_lag) - i + 1)), , drop = FALSE] #BJ
    if(sum_up){
      geom_lag <- geom_lag + weights[i]*lag_i #BJ
    }else{
      geom_lag[,,i] <- weights[i]*lag_i
    }
  } #BJ
  geom_lag #BJ
} #BJ


#' Transform matrix of first-order lagged observations to matrix of weighted sums of past observation (with Poisson weights)

#' This function modifies the design matrices from first-order lags to weighted lags with Poisson structure.
#' To be passed to \code{hhh4lag} or \code{profile_par_lag} in order to use Poisson lags.
#' @param lag1 a matrix of first lags as usually used in \code{hhh4}.
#' @param par_lag a parameter to steer the lag structure, here \eqn{log(\mu)} where \eqn{\mu} is the parameter of
#' the Poisson distribution characterizing the lag structure; see details of \code{hhh4lag} or \code{profile_par_lag}
#' @param min_lag smallest lag to include; the support of the Poisson form starts only at \code{min_lag}. Defaults to 1.
#' @param max_lag highest lag to include; higher lags are cut off and he remaining weights standardized. Defaults to 5.
#' @param sum_up \code{sum_up = FALSE} returns a more detailed output; for debugging only
#' @export
poisson_lag <- function(lag1, par_lag, min_lag, max_lag, sum_up = TRUE){ #BJ

  mu_lag <- exp(par_lag)

  weights0 <- c(rep(0, min_lag - 1), dpois((min_lag:max_lag) - 1, mu_lag)) #BJ
  weights <- weights0/sum(weights0) #BJ
  pois_lag <- if(sum_up){
    matrix(0, ncol = ncol(lag1), nrow = nrow(lag1))
  }else{
    array(NA, dim = c(dim(lag1), max_lag))
  }

  # first line contains NAs
  for(i in 1:max_lag){ #BJ
    lag_i <- lag1[c(rep(1, i - 1), 1:(nrow(pois_lag) - i + 1)), , drop = FALSE] #BJ
    if(sum_up){
      pois_lag <- pois_lag + weights[i]*lag_i #BJ
    }else{
      pois_lag[,,i] <- weights[i]*lag_i
    }
  } #BJ
  pois_lag #BJ
} #BJ

#' Transform matrix of first-order lagged observations to matrix of weighted sums of past observation
#' (with weights for first and second lags)

#' This function modifies the design matrices from first-order lags to weighted lsum of first and second lags.
#' To be passed to \code{hhh4lag} or \code{profile_par_lag} in order to use AR2-lags.
#' @param lag1 a matrix of first lags as usually used in \code{hhh4}.
#' @param par_lag a parameter to steer the lag structure, here \eqn{logit(p)} where \eqn{p} is the weight of
#' the first lag; see details of \code{hhh4lag} or \code{profile_par_lag}
#' @param min_lag smallest lag to include; the support of the geometric form starts only at \code{min_lag}. Defaults to 1.
#' @param max_lag highest lag to include; higher lags are cut off and he remaining weights standardized. Defaults to 5.
#' @param sum_up \code{sum_up = FALSE} returns a more detailed output; for debugging only
#' @export
ar2_lag <- function(lag1, par_lag, min_lag, max_lag, sum_up = TRUE){ #BJ

  p_lag <- exp(par_lag)/(1 + exp(par_lag))

  weights0 <- c(rep(0, min_lag - 1), p_lag, 1 - p_lag, rep(0, max_lag - min_lag - 1)) #BJ
  weights <- weights0/sum(weights0) #BJ
  ar2_lag <- if(sum_up){
    matrix(0, ncol = ncol(lag1), nrow = nrow(lag1))
  }else{
    array(NA, dim = c(dim(lag1), max_lag))
  }

  # first line contains NAs
  for(i in 1:max_lag){ #BJ
    lag_i <- lag1[c(rep(1, i - 1), 1:(nrow(ar2_lag) - i + 1)), , drop = FALSE] #BJ
    if(sum_up){
      ar2_lag <- ar2_lag + weights[i]*lag_i #BJ
    }else{
      ar2_lag[,,i] <- weights[i]*lag_i
    }
  } #BJ
  ar2_lag #BJ
} #BJ
