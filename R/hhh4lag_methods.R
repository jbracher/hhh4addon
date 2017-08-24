## NOTE: we also apply print.hhh4 in print.summary.hhh4()
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
  if(any(x$use_distr_lag)){
    cat("Distributed lags used; use distr_lag() to check the applied lag distribution and parameters.\n")
    cat("\n")
  }
  invisible(x)
}

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
                use_distr_lag = object$control$ar$use_distr_lag | object$control$ne$use_distr_lag,
                maxEV_range = if (maxEV) unique(range(getMaxEV(object)))))
  class(ret) <- "summary.hhh4lag"
  return(ret)
}

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

#' @export
terms.hhh4lag <- function (x, ...)
{
  if (is.null(x$terms))
    hhh4addon:::interpretControl(x$control,x$stsObj) else x$terms
}

# to do: predict.hhh4
# check whether necessary: update; decompose.hhh4 (this could be very useful!), neOffsetArray
