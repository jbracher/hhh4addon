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

## Updated version: Create function (pars, type = "response") which
## returns the weighted sum of time-lagged counts of neighbours
## (or its derivates, if type = "gradient" or type = "hessian").
## For type="reponse", this is a nTime x nUnits matrix (like Y),
## otherwise a list of such matrices,
## which for the gradient has length length(pars) and
## length(pars)*(length(pars)+1)/2 for the hessian.
## If neweights=NULL (i.e. no NE component in model), the result is always 0.
## offset is a multiplicative offset for \phi_{it}, e.g., the population.
## scale is a nUnit-vector or a nUnit x nUnit matrix scaling neweights.
neOffsetFUN <- function (Y, neweights, scale, normalize,
                         nbmat, data, lag, funct_lag, par_lag, max_lag, use_distr_lag, #BJ: added arguments; may have to re-introduce 'lag'
                         sum_up = TRUE, offset = 1) #BJ: added argument sum_up
{
  if (is.null(neweights)) { # no neighbourhood component
    as.function(alist(...=, 0), envir=.GlobalEnv)
    ## dimY <- dim(Y)
    ## as.function(c(alist(...=),
    ##               substitute(matrix(0, r, c), list(r=dimY[1], c=dimY[2]))),
    ##             envir=.GlobalEnv)
  } else if (is.list(neweights)) { # parametric weights
    wFUN <- surveillance:::scaleNEweights.list(neweights, scale, normalize)
    function (pars, type = "response") {
      name <- switch(type, response="w", gradient="dw", hessian="d2w")
      weights <- wFUN[[name]](pars, nbmat, data)
      ## gradient and hessian are lists if length(pars$d) > 1L
      ## but can be single matrices/arrays if == 1 => _c_onditional lapply
      res <- surveillance:::clapply(weights, function (W)
        offset * hhh4addon:::weightedSumNE(observed = Y, weights = W, lag = lag,
                                          funct_lag = funct_lag, par_lag = par_lag, max_lag = max_lag,
                                          use_distr_lag = use_distr_lag, sum_up = sum_up)) # BJ: distr. lags now done inside of weightedSumNE
      ##<- clapply always returns a list (possibly of length 1)
      if (type=="response") res[[1L]] else res
    }
  } else { # fixed (known) weight structure (0-length pars)
    weights <- surveillance:::scaleNEweights.default(neweights, scale, normalize)
    env <- new.env(hash = FALSE, parent = emptyenv())  # small -> no hash
    if(!sum_up & is.matrix(offset)){ # transform ofset to array if provided as matrix.
      offset <- array(offset, dim = c(nrow(offset), ncol(offset), max_lag))
    }
    env$initoffset <- offset * hhh4addon:::weightedSumNE(Y, weights, lag = lag, #BJ
                                                        funct_lag = funct_lag, #BJ
                                                        par_lag = par_lag, #BJ
                                                        max_lag = max_lag, #BJ
                                                        use_distr_lag = use_distr_lag, #BJ
                                                        sum_up = sum_up) #BJ: distr. lags now done inside of weightedSumNE

    as.function(c(alist(...=), quote(initoffset)), envir=env)
  }
}
