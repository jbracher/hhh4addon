################################################################################
### The following are modified versions of functions from the surveillance package.
### See below the original copyright declaration.
################################################################################

################################################################################
### Part of the surveillance package, http://surveillance.r-forge.r-project.org
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Plot-method(s) for fitted hhh4() models
###
### Copyright (C) 2010-2012 Michaela Paul, 2012-2016 Sebastian Meyer
### $Revision: 1715 $
### $Date: 2016-05-17 10:01:28 +0200 (Die, 17. Mai 2016) $
################################################################################

#' @export
#' @method plot hhh4lag
plot.hhh4lag <- function (x,
                       type = c("fitted", "season", "maxEV", "maps", "ri", "neweights"),
                       ...)
{
  # if(!(type[1] %in% c("fitted", "seas"))) stop("This plot type is not implemented yet for objects of class hhh4lag")
  stopifnot(x$convergence)
  cl <- sys.call()  # not match.call() because plotHHH4_season() has no 'x'
  ## remove the type argument from the call
  if (is.null(names(cl)) && nargs() > 1L) { # unnamed call plot(x, type)
    cl[[3L]] <- NULL  # remove the second argument
  } else {
    cl$type <- NULL
  }
  # some of the plotting functions get confused if dim is increased by one due to fitting of par_lag.
  # therefore set dim[1] to to number of fixed effects excluding par_lag parameters:
  x$dim[1] <- length(fixef.hhh4lag(x))

  # cases where a specific function had to be implemented:
  if(type[1] %in% c("fitted")){
    cl[[1L]] <- as.name(paste("plotHHH4lag", match.arg(type), sep="_"))
    eval(cl, envir = parent.frame())
  }

  # cases where the hhh4 functions still work:
  if(type[1] %in% c("season", "neweights")){
    surveillance:::plot.hhh4(x = x, type = type)
  }

  # cases which are not implemented / cannot be implemented:
  if(type[1] == "maxEV"){
    stop("Plot type maxEV not implemented for hhh4lag objects (concept of maxEV not clearly defined.)")
  }

  if(type[1] == "maps"){
    plotHHH4lag_maps(x = x, type = type, ...)
  }

  if(type[1] %in% c("ri")){
    stop("Plot type ri currently not implemented for hhh4lag objects.)")
  }

  if(!type[1] %in% c("fitted", "season", "neweights", "maxEV", "maps", "ri")){
    stop("Plot type does not exist for hhh4lag objects.")
  }
}
