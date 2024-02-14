# #' This function generates (shifted) discrete Weibull weights which are subsequently used inside of \code{get_weighted_lags}. To be passed
# #' to \code{hhh4_lag} or \code{profile_par_lag} as the \code{control$funct_lag} argument.
# #' @param par_lag a parameter vector of length 2 to steer the lag structure, here \eqn{logit(q)} and \eqn{log(beta)},
# #' where \eqn{q} and \eqn{beta} are the parameters of the discrete Weibull distribution as implemented in the `extraDistr` package.
# #' @param min_lag smallest lag to include; the support of the Poisson form starts only at \code{min_lag}. Defaults to 1.
# #' @param max_lag highest lag to include; higher lags are cut off and he remaining weights standardized. Defaults to 5.
# #' @author Maria Dunbar, Johannes Bracher
# #' @export
# discrete_weibull_lag <- function(par_lag, min_lag, max_lag){
#   if(length(par_lag) != 2){
#     stop("The starting value for par_lag needs to be 2")
#   }
#   # ensure parameters are in the right domains:
#   # first parameter (q) is a probability. Note: this is not the standard parameterization
#   # described on Wikipedia.
#   p_lag <- numeric(2)
#   p_lag[1] <- exp(par_lag[1])/(1 + exp(par_lag[1]))
#   # second parameter (beta) is a positive real number
#   p_lag[2] <- exp(par_lag[2])
#
#   # compute weights:
#   weights0 <- c(rep(0, min_lag - 1),
#                 extraDistr::ddweibull(x = (min_lag : max_lag) - 1,
#                                            shape1 = p_lag[1],
#                                            shape2 = p_lag[2]))
#
#   # standardize weights:
#   weights <- weights0 / sum(weights0)
#   return(weights)
# }

#' #' This function generates (shifted) discrete gamma weights which are subsequently used inside of \code{get_weighted_lags}. To be passed
#' #' to \code{hhh4_lag} or \code{profile_par_lag} as the \code{control$funct_lag} argument.
#' #' @param par_lag a parameter vector of length 2 to steer the lag structure, here \eqn{log(shape)} and \eqn{log(rate)},
#' #' where \eqn{shape} and \eqn{rate} are the parameters of the discrete gamma distribution as implemented in the \code{extraDistr} package.
#' #' @param min_lag smallest lag to include; the support of the Poisson form starts only at \code{min_lag}. Defaults to 1.
#' #' @param max_lag highest lag to include; higher lags are cut off and he remaining weights standardized. Defaults to 5.
#' #' @author Maria Dunbar, Johannes Bracher
#' #' @export
# discrete_gamma_lag <- function(par_lag, min_lag, max_lag){
#   # ensure parameters are in the right domains. Both are strictly positive
#   p_lag <- exp(par_lag)
#   # compute weights:
#   # NB this implementation uses the inverse scale (rate) parameter
#   weights0 <- c(rep(0, min_lag - 1),
#                 extraDistr::ddgamma((min_lag : max_lag) - 1,
#                                     p_lag[1], p_lag[2],
#                                     log = FALSE))
#   # standardize weights:
#   weights <- weights0 / sum(weights0)
#   return(weights)
# }


#' This function generates discretized log-normal weights which are subsequently used inside of \code{get_weighted_lags}. To be passed
#' to \code{hhh4_lag} or \code{profile_par_lag} as the \code{control$funct_lag} argument.
#' @param par_lag a parameter vector of length 2 to steer the lag structure, here \eqn{meanlog} and \eqn{log(sdlog)},
#' where \eqn{meanlog} and \eqn{sdlog} are the parameters of the log-normal distribution.
#' @param min_lag smallest lag to include; the support of the Poisson form starts only at \code{min_lag}. Defaults to 1.
#' @param max_lag highest lag to include; higher lags are cut off and he remaining weights standardized. Defaults to 5.
#' @author Maria Dunbar, Johannes Bracher
#' @export
log_normal_lag <- function(par_lag, min_lag, max_lag){
  # ensure parameters are in the right domains. First is from R, second is strictly positive
  p_lag <- par_lag
  p_lag[2] <- exp(par_lag[2])
  # compute weights from increments of CDF:
  x_temp <- (min_lag - 1):max_lag
  weights0 <- c(rep(0, min_lag - 1),
                diff(plnorm(x_temp, p_lag[1], p_lag[2])))
  weights <- weights0 / sum(weights0)
  return(weights)
}
