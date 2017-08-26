### updated version: calculate the weighted sum of counts of adjacent (or all other) regions
### i.e. the nTime x nUnit matrix with elements ne_ti = sum_j w_jit * y_jt
## W is either a nUnits x nUnits matrix of time-constant weights w_ji
## or a nUnits x nUnits x nTime array of time-varying weights

weightedSumNE <- function (observed, weights, lag, funct_lag, par_lag, max_lag , use_distr_lag, sum_up) #BJ
{
  dimY <- dim(observed)
  nTime <- dimY[1L]
  nUnits <- dimY[2L]
  tY <- t(observed)                     # -> nUnits x nTime

  res <- apply(weights, 2L, function (wi)
    ## if dim(weights)==2 (time-constant weights), length(wi)=nUnits,
    ## if dim(weights)==3, wi is a matrix of size nUnits x nTime
    .colSums(tY * wi, nUnits, nTime, na.rm=TRUE))

  if(use_distr_lag == FALSE){ #BJ case where 'lag' is used
    rbind(matrix(NA_real_, lag, nUnits),
          res[seq_len(nTime-lag),,drop=FALSE])
  }else{ #BJ case where distributed lags are used
    lag1 <- rbind(matrix(NA_real_, 1, nUnits), #BJ: force lag = 1 here
                  res[seq_len(nTime-1),,drop=FALSE]) #BJ force lag = 1 here
    funct_lag(lag1 = lag1, #BJ: transform to geometric lags
              par_lag = par_lag, #BJ
              max_lag = max_lag, #BJ
              sum_up = sum_up) #BJ
  }

}

# to keep uniform: define analoguously for AR
weightedSumAR <- function (observed, lag, funct_lag, par_lag, max_lag, use_distr_lag, sum_up) #BJ
{
  dimY <- dim(observed)
  nTime <- dimY[1L]
  nUnits <- dimY[2L]
  if(use_distr_lag){ #BJ: if distributed lags are used:
    Ym1 <- rbind(matrix(NA_integer_, 1, ncol(observed)), head(observed, nTime - 1)) #BJ: force to first lag
    funct_lag(lag1 = Ym1, par_lag = par_lag, max_lag = max_lag, sum_up = sum_up) #BJ transform to geometric lags
  }else{ #BJ: if regular lags are used
    rbind(matrix(NA_integer_, lag, nUnits), head(observed, nTime - lag))
  }
}
