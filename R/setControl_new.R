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

  ## check lags in "ar" and "ne" components
  for (comp in c("ar", "ne")) {

    #BJ: checking plausibility of lag arguments:
    if(!control[[comp]]$use_distr_lag){ #BJ case where 'lag' is used
      if (!surveillance:::isScalar(control[[comp]]$lag) || control[[comp]]$lag < (comp=="ar"))
        stop("'control$", comp, "$lag' must be a ",
             if (comp=="ar") "positive" else "non-negative", " integer")
      control[[comp]]$lag <- as.integer(control[[comp]]$lag)
      control[[comp]]$mu_lag <- NA # set mu_lag, max_lag to NA to avoid that they are used anywhere
      control[[comp]]$max_lag <- NA
    }else{ # case where 'mu_lag' and 'max_lag' are used
      #BJ: check max_lag and mu_lag
      if (!surveillance:::isScalar(control[[comp]]$max_lag) || control[[comp]]$max_lag < (comp == "ar")) #BJ
        stop("'control$", comp, "$max_lag' must be a ", if (comp == "ar") "positive" else "non-negative", " integer")#BJ
      control[[comp]]$max_lag <- as.integer(control[[comp]]$max_lag)#BJ
      # if(control[[comp]]$mu_lag < 1 | control[[comp]]$mu_lag > 0.8*control[[comp]]$max_lag){
      #   stop("mu_lag has to be larger than 1 and smaler than 0.8*max_lag in ar and ne.")
      # }
      control[[comp]]$lag <- NA # set lag to NA to avoid that it is used anywhere
    }
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
  lags <- c(ar = ifelse(control$ar$use_distr_lag, control$ar$max_lag, control$ar$lag), #BJ
            ne = ifelse(control$ne$use_distr_lag, control$ne$max_lag, control$ne$lag)) # BJ
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
