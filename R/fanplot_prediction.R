#' Interpolate between quantiles to avoid edgy display
#'
#' An auxiliary function used in \code{fanplot_prediction}, \code{fanplot_stationary}
#' @param p A vector of probabilities for which to obtain interpolated quantiles
#' @param ... other arguments passed on to \code{pnbinom}
interpolate_qnbinom <- function(p, ...){
  range_x <- qnbinom(p = range(p), ...)
  x <- seq(from = max(range_x[1] - 1, 0),
           to = range_x[2] + 1)
  probs <- c(0, pnbinom(x, ...))
  x <- c(0, x)
  linapprox <- approxfun(probs, x)
  linapprox(p) # + 0.5 # -0.5 makes that the linear approximation goes
  # through the middle of the points rather than below them
}

#' Display prediction as a fan plot
#'
#' Plots a fanplot to display quantiles of (negative binomial approximations) of the week-wise predictive distributions
#'
#' @param pred the prediction as returned by \code{longterm_prediction_hhh4}
#' (and potentially aggregated using aggregate_prediction)
#' @param unit numeric denoting the unit to display
#' @param probs vector of probabilities: which quantiles shall be displayed in the fan plot?
#' @param interpolate_probs logical: smooth curves by simple interpolation of quantiles
#' @param add_observed logical: shall observed values be added?
#' @param fan.col,ln,rlab graphical parameters passed on to \code{fanplot::fan}
#' @param  pt.col,pt.cex,l.col graphical parameters for display of observed values
#' @param add logical: add to existing plot?
#' @param add_legend logical: shall a color key legend be added?
#' @param width_legend width of box for color key legend in user coordinates
#' @param probs_legend vecor of probabilities to display in the legend
#' @param ylim limit for the y-axis, passed to \code{plot()}
#' @param xlab,ylab axis labels
#' @param return_matrix logical: return matrix passed to \code{fanplot::fan}; useful to make more sophisticated plots.
#' @param ... other arguments passed on to \code{plot()}
#'
#' @return Only if \code{return_matrix} set to \code{TRUE}: the matrix passed to fanplot::fan
#'
#' @examples
#' data("salmonella.agona")
#' # convert old "disProg" to new "sts" data class:
#' salmonella <- disProg2sts(salmonella.agona)
#' control_salmonella <- list(end = list(f = addSeason2formula(~ 1), lag = 1),
#'                            ar = list(f = addSeason2formula(~ 1), lag = 1),
#'                            family = "NegBinM", subset = 6:250)
#' fit_salmonella <- hhh4_lag(salmonella, control_salmonella) # fit model
#' # obtain prediction:
#' pred_mom <- predictive_moments(fit_salmonella, t_condition = 250, lgt = 52)
#' # plot the prediction only:
#' fanplot_prediction(pred_mom, add_legend = TRUE)
#' # or plot it along with the fit:
#' plot(fit_salmonella)
#' fanplot_prediction(pred_mom, add = TRUE) # add fan plot
#'
#' @importFrom fanplot fan
#' @export
fanplot_prediction <- function(pred, unit = 1, probs = 1:99/100,
                               interpolate_probs = TRUE, add_observed = TRUE, add_pred_means = TRUE,
                               fan.col = colorRampPalette(c("darkgreen", "gray90")),
                               pt.col = "red", pt.cex = 0.6, l.col = "black",
                               mean_col = "black", mean_lty = "dashed",
                               ln = NULL, rlab = NULL, add = FALSE,
                               add_legend = FALSE, width_legend = 0.1*(max(pred$timepoints) - min(pred$timepoints))/pred$freq,
                               probs_legend = c(1, 25, 50, 75, 99)/100,
                               ylim = NULL,
                               xlab = "t", ylab  ="No. infected",
                               return_matrix = FALSE, ...){

  if(is.null(pred$condition)){
    stop("Prediction object does not contain condition. This can happen when you
         aggregated the prediction using aggregate_prediction() but did not provide
         the aggregation_matrix_condition argument.")
  }

  if(is.null(pred$mu_matrix)){
    stop("Prediction object does not contain element mu_matrix. Either your aggregated
         prediction does not have a temporal structure or you did not specify
         by_week_in_matrix = TRUE in aggregate_prediction().")
  }

  mu <- pred$mu_matrix[, unit]
  sigma2 <- pred$var_matrix[, unit]
  # calculate size parameters for negative binomial approximation:
  psi_ci <- pmin(abs(mu / (sigma2 / mu - 1)), 10000)
  # store 5, 10, 15, ..., 90, 90% probs of this negbin approximation
  matr_cond <- matrix(NA, ncol = length(mu) + 1, nrow = length(probs))

  # need to handle the first time point separately:
  matr_cond[, 1] <- pred$condition[unit]
  # rest is calculated using probs():
  for(i in 2:ncol(matr_cond)){
    if(interpolate_probs){
      matr_cond[, i] <- interpolate_qnbinom(probs, mu = mu[i - 1], size = psi_ci[i - 1])
    }else{
      matr_cond[, i] <- qnbinom(probs, mu = mu[i - 1], size = psi_ci[i - 1])
    }
  }

  timepoints_calendar0 <- c(pred$timepoints_calendar[1] - 1/pred$freq, pred$timepoints_calendar)

  # new plot if add==FALSE
  if(!add){
    if(is.null(ylim)){
      ylim <- c(0, max(matr_cond))
    }

    plot(NULL, xlim = range(timepoints_calendar0) + c(0, 2*add_legend*width_legend),
         ylim = ylim, xlab = xlab, ylab = ylab, ...)
    if(add_legend){
      y_legend <- matrix(seq(from = 0.9*min(matr_cond), to = max(matr_cond),
                             length.out = length(probs)),
                         ncol = 1)
      fan(y_legend, start = max(timepoints_calendar0) + 1.5*width_legend, fan.col = fan.col,
          rlab = probs_legend)
      abline(v = max(timepoints_calendar0) + width_legend)
    }
  }

  # add fan:
  par(bty = "n") # suppress box that fan adds by default
  # plot fan:
  fan(matr_cond, start = timepoints_calendar0[1], frequency = pred$freq,
      fan.col = fan.col, ln = ln, rlab = rlab, data.type = "values", probs = probs)
  abline(v = timepoints_calendar0[1], lty = "dashed")
  # set par()$bty back to default
  par(bty = "o")
  # add means:
  if(add_pred_means){
    lines(timepoints_calendar0, c(pred$condition[nrow(pred$condition), unit], mu),
          col = mean_col, lty = mean_lty)
  }
  # add observed:
  if(add_observed){
    rlz <- c(pred$condition[nrow(pred$condition), unit], pred$realizations_matrix[, unit])
    # lines(timepoints_calendar0, rlz, lty = 2, col = l.col)
    points(timepoints_calendar0, rlz, pch = 19, cex = pt.cex, col = pt.col)
  }
  if(return_matrix){
    return(matr_cond)
  }
}
