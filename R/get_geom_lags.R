#' Display the function and parameters used for distributed lags
#'@export
#'@param hhh4Obj The hhh4 object for which the parameters used for distributed lags should be displayed.
distr_lag <- function(hhh4Obj){
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
