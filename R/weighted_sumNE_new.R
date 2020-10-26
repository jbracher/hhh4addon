################################################################################
### The following are modified versions of functions from the surveillance package
### and wrappers around them.
### See below the original copyright declaration.
################################################################################

################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Helper functions for neighbourhood weight matrices in hhh4()
###
### Copyright (C) 2012-2016 Sebastian Meyer
### $Revision: 1687 $
### $Date: 2016-04-01 21:40:25 +0200 (Fre, 01. Apr 2016) $
################################################################################

### updated version: calculate the weighted sum of counts of adjacent (or all other) regions
### i.e. the nTime x nUnit matrix with elements ne_ti = sum_j w_jit * y_jt
## W is either a nUnits x nUnits matrix of time-constant weights w_ji
## or a nUnits x nUnits x nTime array of time-varying weights

weightedSumNE <- function(observed, weights, lag, funct_lag, par_lag, min_lag, max_lag, sum_up)
{
  dimY <- dim(observed)
  nTime <- dimY[1L]
  nUnits <- dimY[2L]

  if(!sum_up){
    warning("sum_up = FALSE has been deprecated")
    NULL
  }

  #BJ: weight with lag weights prior to applying matrix weights
  lag_weights <- funct_lag(par_lag = par_lag, min_lag = min_lag,
                           max_lag = max_lag)
  lag_weighted_observed <- hhh4addon:::get_weighted_lags(lag1 = observed, lag_weights = lag_weights,
                                                         sum_up = sum_up)

  #BJ: similar to surveillance implementation, but takes lag-weighted version of observations
  if (length(dim(weights)) == 2L) { # fast track for time-constant weights
    # if (any(isNA <- is.na(lag_weighted_observed)))
      # lag_weighted_observed[isNA] <- 0  # keep original na.rm = TRUE behaviour (for now)
      rbind(matrix(NA_real_, 1, nUnits),
            lag_weighted_observed[seq_len(nTime-1),,drop=FALSE] %*% weights)
  } else {
    tYlagged <- t(lag_weighted_observed[seq_len(nTime-1),,drop=FALSE])
    apply(weights[,,(1+1L):nTime,drop=FALSE], 2L, function (wi)
      ## wi and tYlagged are matrices of size nUnits x (nTime-lag)
      c(rep(NA_real_, 1),
        .colSums(tYlagged * wi, nUnits, nTime-1, na.rm=TRUE)))
  }
}

# to keep uniform: define analoguously for AR
weightedSumAR <- function (observed, lag, funct_lag, par_lag, min_lag, max_lag, sum_up) #BJ
{
  dimY <- dim(observed)
  nTime <- dimY[1L]
  nUnits <- dimY[2L]
  #BJ: distributed lags are used:
  Ym1 <- rbind(matrix(NA_integer_, 1, ncol(observed)), head(observed, nTime - 1)) #BJ: force to first lag
  lag_weights <- funct_lag(par_lag = par_lag, min_lag = min_lag, max_lag = max_lag)
  get_weighted_lags(lag1 = Ym1, #BJ: transform to geometric lags
                    lag_weights = lag_weights,
                    sum_up = sum_up) #BJ
}
