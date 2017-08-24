plotHHH4lag_stationary <- function(hhh4Obj, unit, stat_mom = NULL){
  if(is.null(stat_mom)){
    stat_mom <- hhh4addon:::stationary_moments(hhh4Obj, return_mu_split = TRUE)
  }

  meanHHHunit <- stat_mom$mu_split[,unit,]
  inds <- rep_len(1:hhh4Obj$stsObj@freq, nrow(hhh4Obj$stsObj))
  meanHHHunit <- meanHHHunit[inds, ]

  tp <- hhh4Obj$stsObj@start[1] + (hhh4Obj$stsObj@start[2] + 1:length(inds))/hhh4Obj$stsObj@freq # all observation time points

  surveillance:::plotComponentPolygons(
    x = tp,
    y = meanHHHunit[,c("endemic", "epi.own", "epi.neighbours"),drop=FALSE])
}










plotHHH4_stationary1 <- function(x, unit=1, main=NULL,
                             col=c("grey85", "blue", "orange"),
                             pch=19, pt.cex=0.6, pt.col=1, border=col,
                             start=x$stsObj@start, end=NULL, xaxis=NULL,
                             xlim=NULL, ylim=NULL, xlab="", ylab="No. infected",
                             hide0s=FALSE, decompose=NULL, meanHHH=NULL)
{
  stsObj <- x$stsObj
  if (is.character(unit) &&
      is.na(unit <- match(.unit <- unit, colnames(stsObj))))
    stop("region '", .unit, "' does not exist")
  if (is.null(main)) main <- colnames(stsObj)[unit]
  if (isTRUE(decompose)) decompose <- colnames(stsObj)

  ## get observed counts
  obs <- observed(stsObj)[,unit]

  ## time range for plotting
  start0 <- surveillance:::yearepoch2point(stsObj@start, stsObj@freq, toleft=TRUE)
  start <- surveillance:::yearepoch2point(start, stsObj@freq)
  tp <- start0 + seq_along(obs)/stsObj@freq # all observation time points
  if (start < start0 || start > tp[length(tp)])
    stop("'start' is not within the time range of 'x$stsObj'")
  end <- if(is.null(end)) tp[length(tp)] else surveillance:::yearepoch2point(end,stsObj@freq)
  stopifnot(start < end)
  tpInRange <- which(tp >= start & tp <= end)            # plot only those
  tpInSubset <- intersect(x$control$subset, tpInRange)   # fitted time points

  ## use time indexes as x-values for use of addFormattedXAxis()
  if (is.list(xaxis)) {
    tp <- seq_along(obs)
    start <- tpInRange[1L]
    end <- tpInRange[length(tpInRange)]
  }

  stat_mom <- hhh4addon:::stationary_moments(x, return_mu_split = TRUE)
  meanHHHunit <- stat_mom$mu_split[,unit,]
  inds <- rep_len(1:stsObj@freq, length(tp))
  meanHHHunit <- meanHHHunit[inds, ]
  meanHHHunit <- meanHHHunit[tpInSubset, ]

  # ## get fitted component means
  # meanHHHunit <- if (is.null(decompose)) {
  #   if (is.null(meanHHH))
  #     meanHHH <- meanHHH(x$coefficients, terms.hhh4(x))
  #   sapply(meanHHH, "[", i=TRUE, j=unit)
  # } else {
  #   if (is.null(meanHHH))
  #     meanHHH <- decompose.hhh4(x)
  #   if (!setequal(decompose, dimnames(meanHHH)[[3L]][-1L]))
  #     stop("'decompose' must be (a permutation of) the fitted units")
  #   meanHHH[,unit,c("endemic",decompose)]
  # }
  # stopifnot(is.matrix(meanHHHunit), !is.null(colnames(meanHHHunit)),
  #           nrow(meanHHHunit) == length(x$control$subset))
  # meanHHHunit <- meanHHHunit[x$control$subset %in% tpInRange,,drop=FALSE]
  # if (any(is.na(meanHHHunit))) { # -> polygon() would be wrong
  #   ## could be due to wrong x$control$subset wrt the epidemic lags
  #   ## a workaround is then to set 'start' to a later time point
  #   stop("predicted mean contains missing values")
  # }

  ## check color vector
  col <- if (is.null(decompose) && length(col) == 4L) {
    ## compatibility with surveillance < 1.10-0
    pt.col <- col[4L]
    rev(col[-4L])
  } else {
    surveillance:::plotHHH4_fitted_check_col_decompose(col, decompose, colnames(meanHHHunit)[-1L])
  }

  ## establish basic plot window
  if (is.null(ylim)) ylim <- c(0, max(obs[tpInRange],na.rm=TRUE))
  plot(c(start,end), ylim, xlim=xlim, xlab=xlab, ylab=ylab, type="n",
       xaxt = if (is.list(xaxis)) "n" else "s")
  if (is.list(xaxis)) do.call("addFormattedXAxis", c(list(x = stsObj), xaxis))
  title(main=main, line=0.5)

  ## draw polygons
  if (is.null(decompose)) {
    non0 <- which(c("end", "ar", "ne") %in% surveillance:::componentsHHH4(x))
    surveillance:::plotComponentPolygons(
      x = tp[tpInSubset],
      y = meanHHHunit[,c("endemic", "epi.own", "epi.neighbours")[non0],drop=FALSE],
      col = col[non0], border = border[non0], add = TRUE)
  } else {
    non0 <- apply(X = meanHHHunit > 0, MARGIN = 2L, FUN = any)
    surveillance:::plotComponentPolygons(x = tp[tpInSubset], y = meanHHHunit[, non0, drop = FALSE],
                          col = col[non0], border = border[non0], add = TRUE)
  }

  ## add observed counts within [start;end]
  ptidx <- if (hide0s) intersect(tpInRange, which(obs > 0)) else tpInRange
  points(tp[ptidx], obs[ptidx], col=pt.col, pch=pch, cex=pt.cex)

  ## invisibly return the fitted component means for the selected region
  invisible(meanHHHunit)
}
