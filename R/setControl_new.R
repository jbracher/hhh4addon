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

## updated version: set default values for model specifications in control
setControl <- function (control, stsObj)
{
  stopifnot(is.list(control))
  nTime <- nrow(stsObj)
  nUnit <- ncol(stsObj)
  if(nTime <= 2) stop("too few observations")

  ## arguments in 'control' override any corresponding default arguments
  defaultControl <- eval(formals(hhh4_lag)$control) #BJ: adapt to new function name
  environment(defaultControl$ar$f) <- environment(defaultControl$ne$f) <-
    environment(defaultControl$end$f) <- .GlobalEnv
  control <- modifyList(defaultControl, control)

  ## check that component specifications are list objects
  for (comp in c("ar", "ne", "end")) {
    if(!is.list(control[[comp]])) stop("'control$", comp, "' must be a list")
  }

  { # case where 'min_lag' and 'max_lag' are used
    #BJ: check max_lag and mu_lag
    if (!surveillance:::isScalar(control$max_lag) || control$max_lag < (comp == "ar")) #BJ
      stop("'control$", comp, "$max_lag' must be a ", if (comp == "ar") "positive" else "non-negative", " integer")#BJ
    control$max_lag <- as.integer(control$max_lag)#BJ
    if (!surveillance:::isScalar(control$min_lag) || control$min_lag < (comp == "ar")) #BJ
      stop("'control$", comp, "$min_lag' must be a ", if (comp == "ar") "positive" else "non-negative", " integer")#BJ
    control$min_lag <- as.integer(control$min_lag)#BJ
  }

  ### check AutoRegressive component

  if (control$ar$isMatrix <- is.matrix(control$ar$f)) {
    ## this form is not implemented -> will stop() in interpretControl()
    if (any(dim(control$ar$f) != nUnit))
      stop("'control$ar$f' must be a square matrix of size ", nUnit)
    if (is.null(control$ar$weights)) { # use identity matrix
      control$ar$weights <- diag(nrow=nUnit)
    } else if (!is.matrix(control$ar$weights) ||
               any(dim(control$ar$weights) != nUnit)) {
      stop("'control$ar$weights' must be a square matrix of size ", nUnit)
    }
    control$ar$inModel <- TRUE
  } else if (inherits(control$ar$f, "formula")) {
    if (!is.null(control$ar$weights)) {
      warning("argument 'control$ar$weights' is not used")
      control$ar$weights <- NULL
    }
    # check if formula is valid
    control$ar$inModel <- surveillance:::isInModel(control$ar$f)
  } else {
    stop("'control$ar$f' must be either a formula or a matrix")
  }


  ### check NEighbourhood component

  if (!inherits(control$ne$f, "formula"))
    stop("'control$ne$f' must be a formula")
  control$ne$inModel <- surveillance:::isInModel(control$ne$f)

  if (control$ne$inModel) {
    if (nUnit == 1)
      stop("\"ne\" component requires a multivariate 'stsObj'")
    ## if ar$f is a matrix it includes neighbouring units => no "ne" component
    if (control$ar$isMatrix)
      stop("there must not be an extra \"ne\" component ",
           "if 'control$ar$f' is a matrix")
    ## check ne$weights specification
    surveillance:::checkWeights(control$ne$weights, nUnit, nTime,
                                neighbourhood(stsObj), control$data,
                                check0diag = control$ar$inModel)
    ## check optional scaling of weights
    if (!is.null(control$ne$scale)) {
      stopifnot(is.numeric(control$ne$scale))
      if (is.vector(control$ne$scale)) {
        stopifnot(length(control$ne$scale) == 1L ||
                    length(control$ne$scale) %% nUnit == 0,
                  !is.na(control$ne$scale))
      } else {
        surveillance:::checkWeightsArray(control$ne$scale, nUnit, nTime)
      }
    }
  } else {
    control$ne[c("weights", "scale", "normalize")] <- list(NULL, NULL, FALSE)
  }


  ### check ENDemic component

  if (!inherits(control$end$f, "formula"))
    stop("'control$end$f' must be a formula")
  control$end$inModel <- surveillance:::isInModel(control$end$f)


  ### check offsets

  for (comp in c("ar", "ne", "end")) {
    if (is.matrix(control[[comp]]$offset) && is.numeric(control[[comp]]$offset)){
      if (!identical(dim(control[[comp]]$offset), dim(stsObj)))
        stop("'control$",comp,"$offset' must be a numeric matrix of size ",
             nTime, "x", nUnit)
      if (any(is.na(control[[comp]]$offset)))
        stop("'control$",comp,"$offset' must not contain NA values")
    } else if (!identical(as.numeric(control[[comp]]$offset), 1)) {
      stop("'control$",comp,"$offset' must either be 1 or a numeric ",
           nTime, "x", nUnit, " matrix")
    }
  }


  ### stop if no component is included in the model

  if (length(comps <- surveillance:::componentsHHH4(list(control=control))) == 0L)
    stop("none of the components 'ar', 'ne', 'end' is included in the model")


  ### check remaining components of the control list

  if (is.factor(control$family)) {
    stopifnot(length(control$family) == nUnit)
    control$family <- droplevels(control$family)
    names(control$family) <- colnames(stsObj)
  } else {
    control$family <- match.arg(control$family, defaultControl$family)
  }

  if (!is.vector(control$subset, mode="numeric") ||
      !all(control$subset %in% seq_len(nTime)))
    stop("'control$subset' must be %in% 1:", nTime)
  #BJ: use lags or max_lags depending on setting
  lags <- c(ar = control$max_lag, #BJ
            ne = control$max_lag) # BJ
  maxlag <- suppressWarnings(max(lags[names(lags) %in% comps])) # could be -Inf
  if (control$subset[1L] <= maxlag) {
    warning("'control$subset' should be > ", maxlag, " due to epidemic lags")
  }

  if (!is.list(control$optimizer) ||
      any(! sapply(c("stop", "regression", "variance"),
                   function(x) is.list(control$optimizer[[x]]))))
    stop("'control$optimizer' must be a list of lists")

  control$verbose <- as.integer(control$verbose)
  if (length(control$verbose) != 1L || control$verbose < 0)
    stop("'control$verbose' must be a logical or non-negative numeric value")

  stopifnot(is.list(control$start))
  control$start <- local({
    defaultControl$start[] <- control$start[names(defaultControl$start)]
    defaultControl$start
  })
  if (!all(vapply(X = control$start,
                  FUN = function(x) is.null(x) || is.vector(x, mode="numeric"),
                  FUN.VALUE = TRUE, USE.NAMES = FALSE)))
    stop("'control$start' must be a list of numeric start values")

  stopifnot(length(control$keep.terms) == 1L, is.logical(control$keep.terms))

  ## Done
  return(control)
}
