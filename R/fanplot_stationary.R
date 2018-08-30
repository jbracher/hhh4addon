#' Display stationary distribution as a fanplot
#'
#' Plots a fanplot to display quantiles of (negative binomial approximations) of the week-wise stationary distributions
#'
#' @param stat_mom the stationary moments as returned by \code{stationary_moments_hhh4}
#' (and potentially aggregated using aggregate_prediction)
#' @param unit numeric denoting the unit to display
#' @param probs vector of probabilities: which quantiles shall be displayed in the fan plot?
#' @param interpolate_probs logical: smooth curves by simple interpolation of quantiles
#' @param add_pred_means logical: add line showing the the predictive means
#' @param timepoints vector giving the x-coordinates for the fanplot (generates \code{start} and \code{frequency} for \code{fanplot::fan})
#' @param fan.col,ln,ln.col,rlab,style graphical parameters passed on to \code{fanplot::fan}
#' @param  pt.col,pt.cex,l.col graphical parameters for display of observed values
#' @param means_col,mean_lty graphical parameters for display of predictive means
#' @param add logical: add to existing plot?
#' @param add_legend logical: shall a color key legend be added?
#' @param width_legend width of box for color key legend in user coordinates
#' @param probs_legend vecor of probabilities to display in the legend
#' @param hlines,vlines coordinates for horizontal and vertical grid lines
#' @param ylim limit for the y-axis, passed to \code{plot()}
#' @param xlab,ylab axis labels
#' @param return_matrix logical: return matrix passed to \code{fanplot::fan}; useful to make more sophisticated plots.
#' @param ... other arguments passed on to \code{plot()}
#'
#' @return Only if \code{return_matrix} set to \code{TRUE}: the matrix passed to fanplot::fan
#'
#' @examples
#' data("salmonella.agona")
#' ## convert old "disProg" to new "sts" data class
#' salmonella <- disProg2sts(salmonella.agona)
#' # specify and fit model
#' control_salmonella <- list(end = list(f = addSeason2formula(~ 1), lag = 1),
#'                            ar = list(f = addSeason2formula(~ 1), lag = 1),
#'                            family = "NegBinM")
#' fit_salmonella <- hhh4(salmonella, control_salmonella)
#' # obtain periodically stationary moments:
#' stat_mom <- stationary_moments(fit_salmonella)
#' # plot periodically stationary means:
#' fanplot_stationary(stat_mom, add_legend = TRUE)
#' # add paths of the six seasons in the data set:
#' for(i in 0:5){
#'  lines(1:52/52, salmonella@observed[(i*52 + 1):((i + 1)*52)], col = "blue")
#' }
#' legend("topleft", col = "blue", lty = 1, legend = "observed seasons")
#' @importFrom fanplot fan
#' @export
fanplot_stationary <- function(stat_mom, unit = 1, probs = 1:99/100,
                               interpolate_probs = TRUE, add_pred_means = TRUE,
                               fan.col = colorRampPalette(c("darkgreen", "gray90")), pt.col = "red", pt.cex = 0.3, l.col = "black",
                               mean_col = "black", mean_lty = "dashed",
                               ln = NULL, ln.col = "red", rlab = NULL, style = "fan", add = FALSE,
                               timepoints = 1:nrow(stat_mom$mu_matrix) / stat_mom$freq,
                               add_legend = FALSE, width_legend = 0.1*(max(timepoints) - min(timepoints)),
                               probs_legend = c(1, 25, 50, 75, 99)/100, hlines = NULL, vlines = NULL,
                               ylim = NULL,
                               xlab = "t", ylab  = "No. infected",
                               return_matrix = FALSE, ...){

  if(is.null(stat_mom$mu_matrix)){
    stop("stat_mom object does not contain element mu_matrix.")
  }

  mu <- stat_mom$mu_matrix[, unit]
  sigma2 <- stat_mom$var_matrix[, unit]
  # calculate size parameters for negative binomial approximation:
  psi_ci <- mu / (sigma2 / mu - 1)
  # store 5, 10, 15, ..., 90, 90% probs of this negbin approximation
  matr_cond <- matrix(NA, ncol = length(mu), nrow = length(probs))

  # calculate distributions:
  for(i in 1:ncol(matr_cond)){
    if(interpolate_probs){
      matr_cond[, i] <- interpolate_qnbinom(probs, mu = mu[i], size = psi_ci[i])
    }else{
      matr_cond[, i] <- qnbinom(probs, mu = mu[i], size = psi_ci[i])
    }
  }

  # new plot if add==FALSE
  if(!add){
    if(is.null(ylim)){
      ylim <- c(0, max(matr_cond))
    }
    plot(NULL, xlim = range(timepoints) + c(0, 2*add_legend*width_legend),
         ylim = ylim, xlab = xlab, ylab = ylab, ...)
    if(!is.null(hlines)){
        add_hlines(y = hlines, xlim = c(timepoints[1], max(timepoints) + width_legend), col = "lightgrey")
    }
    if(!is.null(vlines)){
      abline(v = vlines, col = "lightgrey")
    }
    if(add_legend){
      y_legend <- matrix(seq(from = 1.1*min(ylim), to = 0.95*max(ylim),
                             length.out = length(probs)),
                         ncol = 1)
      fan(y_legend, start = max(timepoints) + 1.5*width_legend, fan.col = fan.col,
          rlab = probs_legend)
      abline(v = max(timepoints) + width_legend)
    }
  }

  # add fan:
  par(bty = "n") # suppress box that fan adds by default
  # plot fan:
  fan(matr_cond, start = timepoints[1], frequency = 1/(timepoints[2] - timepoints[1]),
      fan.col = fan.col, ln = ln, ln.col = ln.col, rlab = rlab,
      data.type = "values", probs = probs,
      style = style)
  if(add_pred_means){
    lines(timepoints, mu, col = mean_col, lty = mean_lty)
  }
  # set par()$bty back to default
  par(bty = "o")
  if(return_matrix){
    return(matr_cond)
  }
}

add_hlines <- function(y, xlim, ...){
  for(i in y){
    lines(xlim, c(i, i), ...)
  }
}
