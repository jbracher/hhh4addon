################################################################################
### The following are modified versions of functions from the surveillance package.
### See below the original copyright declaration.
################################################################################

################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Standard methods for hhh4-fits
###
### Copyright (C) 2010-2012 Michaela Paul, 2012-2016 Sebastian Meyer
### $Revision: 1697 $
### $Date: 2016-04-06 14:21:54 +0200 (Mit, 06. Apr 2016) $
################################################################################



## NOTE: we also apply print.hhh4 in print.summary.hhh4()
#' A modified version of \code{surveillance::print.hhh4}
#'
#' A modified version of \code{surveillance::print.hhh4} to deal with the added
#' features of the \code{hhh4lag} class.
#'@export
print.hhh4lag <- function (x, digits = max(3, getOption("digits")-3), ...)
{
  if (!x$convergence) {
    cat('Results are not reliable! Try different starting values.\n')
    return(invisible(x))
  }
  if (!is.null(x$call)) {
    cat("\nCall: \n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
  }
  if (x$dim["random"] > 0) {
    cat('Random effects:\n')
    .printREmat(if (is.null(x$REmat)) .getREmat(x) else x$REmat,
                digits = digits)
    cat("\nFixed effects:\n")
  } else if (x$dim["fixed"] > 0) {
    cat("Coefficients:\n")
  }
  if (x$dim["fixed"] > 0) {
    print.default(
      format(if (is.null(x$fixef)) surveillance:::fixef.hhh4(x, ...) else x$fixef,
             digits=digits),
      quote = FALSE, print.gap = 2)
  } else cat("No coefficients\n")
  cat("\n")
  # if(x$use_distr_lag){ #BJ
    wgts <- x$distr_lag
    cat(paste0("Distributed lags used (max_lag = ", length(wgts),
               "). Weights: "))
    cat(paste(round(wgts, 2), collapse = "; "))
    cat("\nUse distr_lag() to check the applied lag distribution and parameters.\n") #BJ
    cat("\n") #BJ
  #}
  invisible(x)
}

#' A modified version of \code{surveillance::summary.hhh4}
#'
#' A modified version of \code{surveillance::summary.hhh4} to deal with the added
#' features of the \code{hhh4lag} class.
#' @export
summary.hhh4lag <- function (object, maxEV = FALSE, ...)
{
  ## do not summarize results in case of non-convergence
  if (!object$convergence) {
    cat('Results are not reliable! Try different starting values.\n')
    return(invisible(object))
  }
  ret <- c(object[c("call", "convergence", "dim", "loglikelihood", "margll",
                    "lags", "nTime", "nUnit")],
           list(fixef = surveillance:::fixef.hhh4(object, se=TRUE, ...),
                ranef = surveillance:::ranef.hhh4(object, ...),
                REmat = surveillance:::.getREmat(object),
                AIC   = AIC(object),
                BIC   = BIC(object),
                # use_distr_lag = object$control$use_distr_lag,
                maxEV_range = if (maxEV) unique(range(getMaxEV(object))),
                distr_lag = object$distr_lag))
  class(ret) <- "summary.hhh4lag"
  return(ret)
}

#' A modified version of \code{print.summary.hhh4}
#'
#' A modified version of \code{print.summary.hhh4} to deal with the added
#' features of the \code{hhh4lag} class.
#' @export
print.summary.hhh4lag <- function (x, digits = max(3, getOption("digits")-3), ...)
{
  ## x$convergence is always TRUE if we have a summary
  print.hhh4lag(x) # also works for summary.hhh4-objects

  if (!is.null(x$maxEV_range))
    cat("Epidemic dominant eigenvalue: ",
        paste(sprintf("%.2f", x$maxEV_range), collapse = " -- "), "\n\n")
  if(x$dim["random"]==0){
    cat('Log-likelihood:  ',round(x$loglikelihood,digits=digits-2),'\n')
    cat('AIC:             ',round(x$AIC,digits=digits-2),'\n')
    cat('BIC:             ',round(x$BIC,digits=digits-2),'\n\n')
  } else {
    cat('Penalized log-likelihood: ',round(x$loglikelihood,digits=digits-2),'\n')
    cat('Marginal log-likelihood:  ',round(x$margll,digits=digits-2),'\n\n')
  }
  cat('Number of units:       ', x$nUnit, '\n')
  cat('Number of time points: ', x$nTime, '\n')
  if (!is.null(x$lags)) { # only available since surveillance 1.8-0
    if (!is.na(x$lags["ar"]) && x$lags["ar"] != 1)
      cat("Non-default autoregressive lag:  ", x$lags[["ar"]], "\n")
    if (!is.na(x$lags["ne"]) && x$lags["ne"] != 1)
      cat("Non-default neighbor-driven lag: ", x$lags[["ne"]], "\n")
  }
  cat("\n")
  invisible(x)
}

#' A modified version of \code{terms.hhh4}
#'
#' A modified version of \code{terms.hhh4} to deal with the added
#' features of the \code{hhh4lag} class.
#' @export
terms.hhh4lag <- function (x, ...)
{
  if (is.null(x$terms))
    hhh4addon:::interpretControl(x$control,x$stsObj) else x$terms
}

#' A modified version of \code{logLik.hhh4}
#'
#' A modified version of \code{terms.hhh4} to deal with the added
#' features of the \code{logLik.hhh4} class.
#' @export
logLik.hhh4lag <- function(object, ...)
{
  val <- if (object$convergence) object$loglikelihood else {
    warning("algorithm did not converge")
    NA_real_
  }
  attr(val, "df") <- if (object$dim["random"])
    NA_integer_ else object$dim[["fixed"]] # use "[[" to drop the name; # BJ: maybe add + 1 to account for additional lag parameter
  attr(val, "nobs") <- surveillance:::nobs.hhh4(object)
  class(val) <- "logLik"
  val
}


### refit hhh4-model
## ...: arguments modifying the original control list
## S: a named list to adjust the number of harmonics of the three components
## subset.upper: refit on a subset of the data up to that time point
## use.estimates: use fitted parameters as new start values

#' A modified version of \code{update.hhh4}
#'
#' A modified version of \code{update.hhh4} to deal with the added
#' features of the \code{hhh4lag} class.
#' @export
update.hhh4lag <- function (object, ..., S = NULL, subset.upper = NULL,
                         use.estimates = object$convergence, evaluate = TRUE)
{
  control <- object$control

  ## first modify the control list according to the components in ...
  extras <- list(...)
  control <- modifyList(control, extras)

  ## adjust start values
  control$start <- if (use.estimates) { # use parameter estimates
    surveillance:::hhh4coef2start(object)
  } else local({ # re-use previous 'start' specification
    ## for pre-1.8-2 "hhh4" objects,
    ## object$control$start is not necessarily a complete list:
    template <- eval(formals(hhh4)$control$start)
    template[] <- object$control$start[names(template)]
    template
  })
  ## and update according to an extra 'start' argument
  if (!is.null(extras[["start"]])) {
    if (!is.list(extras$start) || is.null(names(extras$start))) {
      stop("'start' must be a named list, see 'help(\"hhh4\")'")
    }
    control$start[] <- mapply(
      FUN = function (now, extra) {
        if (is.null(names(extra))) {
          extra
        } else { # can retain non-extra values
          now[names(extra)] <- extra
          now
        }
      },
      control$start, extras$start[names(control$start)],
      SIMPLIFY = FALSE, USE.NAMES = FALSE
    )
  }
  ## update initial values of parametric weight function
  if (use.estimates && length(coefW <- coefW(object)) &&
      ! "weights" %in% names(extras$ne)) { # only if function is unchanged
    control$ne$weights$initial <- coefW
  }

  ## adjust seasonality
  if (!is.null(S)) {
    stopifnot(is.list(S), !is.null(names(S)),
              names(S) %in% c("ar", "ne", "end"))
    control[names(S)] <- mapply(function (comp, S) {
      comp$f <- surveillance::addSeason2formula(removeSeasonFromFormula(comp$f),
                                  period = object$stsObj@freq, S = S)
      comp
    }, control[names(S)], S, SIMPLIFY=FALSE, USE.NAMES=FALSE)
  }

  ## restrict fit to those epochs of control$subset which are <=subset.upper
  if (surveillance:::isScalar(subset.upper)) {
    if (subset.upper < control$subset[1L])
      stop("'subset.upper' is smaller than the lower bound of 'subset'")
    control$subset <- control$subset[1L]:subset.upper
  }


  ## fit the updated model or just return the modified control list
  if (evaluate) {
    hhh4_lag(stsObj = object$stsObj, control = control)
  } else {
    control
  }
}

## adapted version: decompose the fitted mean of a "hhh4" model returning an array
## with dimensions (t, i, j), where the first j index is "endemic"

#' A wrapper around \code{decompose.hhh4lag} and \code{surveillance::decompose.hhh4}
#'
#' A wrapper around \code{decompose.hhh4lag} and \code{surveillance::decompose.hhh4}
#' to handle ordinary \code{hhh4} objects and objects of the new \code{hhh4lag} class.
#'
#' @export
decompose.hhh4 <- function(x, coefs = x$coefficients, ...){
  if(class(x)[1] == "hhh4lag"){
    decompose.hhh4lag(x, coefs = x$coefficients, ...)
  }else{
    surveillance::decompose.hhh4(x, coefs = x$coefficients, ...)
  }
}

#' A modified version of \code{decompose.hhh4}
#'
#' A modified version of \code{decompose.hhh4} to deal with the added
#' features of the \code{hhh4lag} class.
decompose.hhh4lag <- function (x, coefs = x$coefficients, ...)
{
  ## get three major components from meanHHH() function
  meancomps <- surveillance:::meanHHH(coefs, terms.hhh4lag(x))

  ## this contains c("endemic", "epi.own", "epi.neighbours")
  ## but we really want the mean by neighbour
  neArray <- c(meancomps$ne.exppred) * neOffsetArray.hhh4lag(x, coefW(coefs))
  ##<- ne.exppred is (t, i) and recycled for (t, i, j)
  stopifnot(all.equal(rowSums(neArray, dims = 2), meancomps$epi.neighbours,
                      check.attributes = FALSE))

  ## add autoregressive part to neArray
  diagidx <- cbind(c(row(meancomps$epi.own)),
                   c(col(meancomps$epi.own)),
                   c(col(meancomps$epi.own)))
  ## usually: neArray[diagidx] == 0
  neArray[diagidx] <- neArray[diagidx] + meancomps$epi.own

  ## add endemic component to the array
  res <- array(c(meancomps$endemic, neArray),
               dim = dim(neArray) + c(0, 0, 1),
               dimnames = with(dimnames(neArray), list(t=t, i=i, j=c("endemic",j))))
  stopifnot(all.equal(rowSums(res, dims = 2), meancomps$mean,
                      check.attributes = FALSE))
  res
}

## new version:
## get the w_{ji} Y_{j,t-1} values from a hhh4() fit
## (i.e., before summing the neighbourhood component over j)
## in an array with dimensions (t, i, j)

#' A modified version of \code{neOffsetArray}
#'
#' A modified version of \code{neOffsetArray} to deal with the added
#' features of the \code{hhh4lag} class.
neOffsetArray.hhh4lag <- function (object, pars = coefW(object),
                           subset = object$control$subset)
{
  ## initialize array ordered as (j, t, i) for apply() below
  res <- array(data = 0,
               dim = c(object$nUnit, length(subset), object$nUnit),
               dimnames = list(
                 "j" = colnames(object$stsObj),
                 "t" = rownames(object$stsObj)[subset],
                 "i" = colnames(object$stsObj)))
  # BJ: extract some elements of the control:
  control <- object$control

  ## calculate array values if the fit has an NE component
  if ("ne" %in% surveillance:::componentsHHH4(object)) {
    W <- surveillance:::getNEweights(object, pars = pars)
    Y <- observed(object$stsObj)
    # tm1 <- subset - object$control$ne$lag
    # is.na(tm1) <- tm1 <= 0
    # tYtm1 <- t(Y[tm1,,drop=FALSE])
    #BJ calculate lags using weightedSumAR instead of indexing as in original function
    tY_lagged <- t(hhh4addon:::weightedSumAR(observed = Y, lag = control$ar$lag, #BJ
                                          funct_lag = control$funct_lag, par_lag = control$par_lag,
                                          max_lag = control$max_lag, min_lag = control$min_lag,
                                          sum_up = TRUE)[subset, ]) #BJ
    # from now on everything continues as before
    res[] <- apply(W, 2L, function (wi) tY_lagged * wi)
    offset <- object$control$ne$offset
    res <- if (length(offset) > 1L) {
      offset <- offset[subset,,drop=FALSE]
      res * rep(offset, each = object$nUnit)
    } else {
      res * offset
    }
    ## stopifnot(all.equal(
    ##     colSums(res),  # sum over j
    ##     terms.hhh4(object)$offset$ne(pars)[subset,,drop=FALSE],
    ##     check.attributes = FALSE))
  }

  ## permute dimensions as (t, i, j)
  aperm(res, perm = c(2L, 3L, 1L), resize = TRUE)
}


#' A modified version of \code{residuals.hhh4}
#'
#' A modified version of \code{residuals.hhh4} to deal with the added
#' features of the \code{hhh4lag} class. Computes deviance residuals.
#' @export
residuals.hhh4lag <- function (object, type = c("deviance", "response"), ...)
{
  type <- match.arg(type)
  obs <- observed(object$stsObj)[object$control$subset,]
  fit <- fitted(object)
  if (type == "response")
    return(obs - fit)

  ## deviance residuals
  ## Cf. residuals.ah, it calculates:
  ## deviance = sign(y - mean) * sqrt(2 * (distr(y) - distr(mean)))
  ## pearson = (y - mean)/sqrt(variance)
  dev.resids <- if (identical(object$control$family, "Poisson")) {
    poisson()$dev.resids
  } else {
    size <- if (identical(object$control$family, "NegBin1")) {
      hhh4addon:::psi2size.hhh4(object, subset = NULL) # changed
    } else {
      hhh4addon:::psi2size.hhh4lag(object) # CAVE: a matrix -> non-standard "size" # changed
    }
    MASS:::negative.binomial(size)$dev.resids
  }

  di2 <- dev.resids(y=obs, mu=fit, wt=1)
  sign(obs-fit) * sqrt(pmax.int(di2, 0))
}

#' A modified version of \code{psi2size.hhh4}
#'
#' A modified version of \code{psi2size.hhh4} to deal with the added
#' features of the \code{hhh4lag} class. Extracts estimated overdispersion
#' in dnbinom() parametrization (and as matrix)
#' @export
##
psi2size.hhh4lag <- function (object, subset = object$control$subset, units = NULL)
{
  size <- sizeHHH(object$coefficients, hhh4addon:::terms.hhh4lag(object), subset = subset) # only change
  if (!is.null(size) && !is.null(units)) {
    if (is.null(subset)) {
      warning("ignoring 'units' (not compatible with 'subset = NULL')")
      size
    } else {
      size[, units, drop = FALSE]
    }
  } else {
    size
  }
}
