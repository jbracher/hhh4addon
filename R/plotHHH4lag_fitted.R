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
plotHHH4lag_fitted <- function (x, units = 1, names = NULL,
                             col = c("grey85", "blue", "orange"),
                             pch = 19, pt.cex = 0.6, pt.col = 1,
                             par.settings = list(),
                             legend = TRUE, legend.args = list(),
                             legend.observed = FALSE,
                             decompose = NULL, meanHHH = NULL, ...)
{
  if (is.null(units)) units <- seq_len(x$nUnit)
  if (!is.null(names)) stopifnot(length(units) == length(names))
  if (isTRUE(decompose)) decompose <- colnames(x$stsObj)

  ## get decomposed mean
  if (is.null(meanHHH)) {
    meanHHH <- if (is.null(decompose)) {
      surveillance:::meanHHH(x$coefficients, terms(x)) #BJ: removed hhh4. This was the only thing that prevented the function from working.
    } else {
      surveillance:::decompose.hhh4(x)
    }
  }

  ## check color vector
  col <- if (is.null(decompose) && length(col) == 4) {
    ## compatibility with surveillance < 1.10-0
    pt.col <- col[4L]
    rev(col[-4L])
  } else {
    surveillance:::plotHHH4_fitted_check_col_decompose(col, decompose, dimnames(meanHHH)[[3L]][-1L])
  }

  ## setup graphical parameters
  if (is.list(par.settings)) {
    par.defaults <- list(mfrow = sort(n2mfrow(length(units))),
                         mar = c(4,4,2,0.5)+.1, las = 1)
    par.settings <- modifyList(par.defaults, par.settings)
    opar <- do.call("par", par.settings)
    on.exit(par(opar))
  }

  ## legend options
  if (is.logical(legend)) legend <- which(legend)
  if (!is.list(legend.args)) {
    if (length(legend) > 0)
      warning("ignored 'legend' since 'legend.args' is not a list")
    legend <- integer(0L)
  }
  if (length(legend) > 0) {
    legendidx <- 1L + c(
      if (legend.observed && !is.na(pch)) 0L,
      if (is.null(decompose)) {
        which(c("ne","ar","end") %in% surveillance:::componentsHHH4(x))
      } else seq_along(col))
    default.args <- list(
      x="topright", col=c(pt.col,rev(col))[legendidx], lwd=6,
      lty=c(NA,rep.int(1,length(col)))[legendidx],
      pch=c(pch,rep.int(NA,length(col)))[legendidx],
      pt.cex=pt.cex, pt.lwd=1, bty="n", inset=0.02,
      legend=if (is.null(decompose)) {
        c("observed","spatiotemporal","autoregressive","endemic")[legendidx]
      } else c("observed", rev(decompose), "endemic")[legendidx]
    )
    legend.args <- modifyList(default.args, legend.args)
  }

  ## plot fitted values region by region
  meanHHHunits <- vector(mode="list", length=length(units))
  names(meanHHHunits) <- if (is.character(units)) units else colnames(x$stsObj)[units]
  for(i in seq_along(units)) {
    meanHHHunits[[i]] <- surveillance:::plotHHH4_fitted1(x, units[i], main=names[i],
                                          col=col, pch=pch, pt.cex=pt.cex, pt.col=pt.col,
                                          decompose=decompose, meanHHH=meanHHH, ...)
    if (i %in% legend) do.call("legend", args=legend.args)
  }
  invisible(meanHHHunits)
}
