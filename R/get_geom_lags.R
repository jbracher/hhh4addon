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
#'                            par_lag = 0.8, use_distr_lag = TRUE),
#'                            family = "NegBinM", subset = 6:312)
#' fit_salmonella <- hhh4_lag(salmonella, control_salmonella)
#' distr_lag(fit_salmonella)

distr_lag <- function(hhh4Obj){
  if(!("hhh4lag" %in% class(hhh4Obj))){
    stop("structure of distributed lags can only be extracted from objects of class hhh4lag.")
  }

  ret <- list()
  for(comp in c("ar", "ne")){
    ret[[comp]] <- if(hhh4Obj$control[[comp]]$use_distr_lag){
      list(funct_lag = hhh4Obj$control[[comp]]$funct_lag,
           par_lag = hhh4Obj$control[[comp]]$par_lag,
           max_lag = hhh4Obj$control[[comp]]$max_lag)
    }
  }
  return(ret)
}

# helper function: transform matrix of first lags to matrix of geometric lags
geometric_lag <- function(lag1, par_lag, max_lag, sum_up = TRUE){ #BJ

  if(is.list(par_lag)){
    p_lag <- par_lag$mu # extract parameter
  }else{
    p_lag <- par_lag
  }

  weights0 <- dgeom(0:(max_lag - 1), p_lag) #BJ
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
