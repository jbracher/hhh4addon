###
### Maps of the fitted mean components averaged over time
###
#' @export
plotHHH4lag_maps <- function (x,
                           which = c("mean", "endemic", "epi.own", "epi.neighbours"),
                           prop = FALSE, main = which, zmax = NULL, col.regions = NULL,
                           labels = FALSE, sp.layout = NULL, ...,
                           map = x$stsObj@map, meanHHH = NULL)
{
  which <- match.arg(which, several.ok = TRUE)
  if (is.null(col.regions))
    col.regions <- surveillance:::.hcl.colors(10) # JB

  ## extract district-specific mean components
  if (is.null(meanHHH)) {
    meanHHH <- surveillance:::meanHHH(x$coefficients, hhh4addon:::terms.hhh4lag(x)) # JB
  }

  ## select relevant components and convert to an array
  meanHHH <- simplify2array(
    meanHHH[c("mean", "endemic", "epi.own", "epi.neighbours")],
    higher = TRUE)

  ## convert to proportions
  if (prop) {
    meanHHH[,,-1L] <- meanHHH[,,-1L,drop=FALSE] / c(meanHHH[,,1L])
  }

  ## select only 'which' components
  meanHHH <- meanHHH[,,which,drop=FALSE]

  ## check map
  map <- as(map, "SpatialPolygonsDataFrame")
  if (!all(dimnames(meanHHH)[[2L]] %in% row.names(map))) {
    stop("'row.names(map)' do not cover all fitted districts")
  }

  ## average over time
  comps <- as.data.frame(colMeans(meanHHH, dims = 1))

  ## attach to map data
  map@data <- cbind(map@data, comps[row.names(map),,drop=FALSE])

  ## color key range
  if (is.null(zmax)) {
    zmax <- if (prop) {
      ceiling(10*sapply(comps, max))/10
    } else ceiling(sapply(comps, max))
    ## sub-components should have the same color range
    .idxsub <- setdiff(seq_along(zmax), match("mean", names(zmax)))
    zmax[.idxsub] <- suppressWarnings(max(zmax[.idxsub]))
  }

  ## add sp.layout item for district labels
  if (!is.null(layout.labels <- layout.labels(map, labels))) {
    sp.layout <- c(sp.layout, list(layout.labels))
  }

  ## produce maps
  grobs <- mapply(
    FUN = function (zcol, main, zmax)
      spplot(map, zcol = zcol, main = main,
             at = seq(0, zmax, length.out = length(col.regions) + 1L),
             col.regions = col.regions, sp.layout = sp.layout, ...),
    zcol = names(comps), main = main, zmax = zmax,
    SIMPLIFY = FALSE, USE.NAMES = FALSE)
  if (length(grobs) == 1L) {
    gridExtra::grid.arrange(grobs[[1L]])
  } else {
    mfrow <- sort(n2mfrow(length(grobs)))
    gridExtra::grid.arrange(grobs = grobs, nrow = mfrow[1L], ncol = mfrow[2L])
  }
}
