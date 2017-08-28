#' @export
plot.hhh4lag <- function (x,
                       type = c("fitted", "season", "maxEV", "maps", "ri", "neweights"),
                       ...)
{
  if(type[1] != "fitted") stop("This plot type is not implemented yet for objects of class hhh4lag")
  stopifnot(x$convergence)
  cl <- sys.call()  # not match.call() because plotHHH4_season() has no 'x'
  ## remove the type argument from the call
  if (is.null(names(cl)) && nargs() > 1L) { # unnamed call plot(x, type)
    cl[[3L]] <- NULL  # remove the second argument
  } else {
    cl$type <- NULL
  }
  cl[[1L]] <- as.name(paste("plotHHH4lag", match.arg(type), sep="_"))
  eval(cl, envir = parent.frame())
}
