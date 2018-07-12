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

weightedSumNE <- function (observed, weights, lag, funct_lag, par_lag, min_lag, max_lag, sum_up) #BJ
{
  dimY <- dim(observed)
  nTime <- dimY[1L]
  nUnits <- dimY[2L]
  tY <- t(observed)                     # -> nUnits x nTime

  res <- apply(weights, 2L, function (wi)
    ## if dim(weights)==2 (time-constant weights), length(wi)=nUnits,
    ## if dim(weights)==3, wi is a matrix of size nUnits x nTime
    .colSums(tY * wi, nUnits, nTime, na.rm=TRUE))

  #BJ distributed lags are used
  lag1 <- rbind(matrix(NA_real_, 1, nUnits), #BJ: force lag = 1 here
                res[seq_len(nTime-1),,drop=FALSE]) #BJ force lag = 1 here
  funct_lag(lag1 = lag1, #BJ: transform to geometric lags
            par_lag = par_lag, #BJ
            min_lag = min_lag, #BJ
            max_lag = max_lag, #BJ
            sum_up = sum_up) #BJ
}

# to keep uniform: define analoguously for AR
weightedSumAR <- function (observed, lag, funct_lag, par_lag, min_lag, max_lag, sum_up) #BJ
{
  dimY <- dim(observed)
  nTime <- dimY[1L]
  nUnits <- dimY[2L]
  #BJ: distributed lags are used:
  Ym1 <- rbind(matrix(NA_integer_, 1, ncol(observed)), head(observed, nTime - 1)) #BJ: force to first lag
  funct_lag(lag1 = Ym1, par_lag = par_lag, min_lag = min_lag, max_lag = max_lag, sum_up = sum_up) #BJ transform to geometric lags

}
