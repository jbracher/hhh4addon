#' Plot predictive or stationary moments by unit
#'
#' Plot negative binomial approximation of predictive or stationary distributions.
#' Usually to be used with aggregated predictions (where columns correspond to regions
#' or age groups; no temporal structure kept).
#'
#' @param mom an object of class \code{predictive_moments_hhh4} or
#' \code{stationary_moments_hhh4}, usually an aggregated prediction
#' without temporal structure (aggregated using \code{aggregate_prediction})
#' @param probs probabilities displayed in fanplot (passed to \code{fanplot::fan})
#' @param add_observed should obseved values be added to the plot?
#' @param fan.col color palette for fanplot (passed to \code{fanplot::fan})
#' @param ln,rlab,style,space additional arguments passed to \code{fanplot::fan})
#' @param pt.col,pt.cex point color and size for observations
#' @param mean_col,mean_lty line color and type for predictive/stationary means
#' @param add_legend should a legend with the colour coding be added?
#' @param probs_legend probabilities to be displayed in the legend
#' @param ylim,main,xlab,las,axes usual plotting parameters
#'
#' @examples
#'
#' # load data:
#' data("noroBL")
#'
#' ########
#' # fit a bivariate model:
#' controlBL <- list(end = list(f = addSeason2formula(~ -1 + fe(1, unitSpecific = TRUE))),
#'                   ar = list(f = ~ -1 + fe(1, unitSpecific = TRUE)),
#'                   ne = list(f = ~ -1 + fe(1, unitSpecific = TRUE)),
#'                   family = "NegBinM", subset = 2:260) # not a very parsimonious parametrization, but feasible
#' fitBL <- hhh4(noroBL, control = controlBL)
#' pred_mom <- predictive_moments(fitBL, t_condition = 260, lgt = 52, return_Sigma = TRUE)
#' # Sigma is required in order to aggregate predictions.
#'
#' #########
#' # perform an aggregation over time: total burden in the two regions
#' aggr_matr_total_burden <- matrix(rep(c(1, 0, 0, 1), 52), nrow = 2,
#'                                  dimnames = list(c("Bremen", "Lower Saxony"),
#'                                                  NULL))
#' pred_mom_total_burden <- aggregate_moments(pred_mom, aggr_matr_total_burden)
#' plot_moments_by_unit(pred_mom_total_burden, main = "Total burden 2016", add_legend = TRUE)

#' @export

plot_moments_by_unit <- function (mom, probs = 1:99/100, add_observed = TRUE, add_pred_means = TRUE,
                                  fan.col = colorRampPalette(c("darkgreen", "gray90")),
                                  pt.col = "red", pt.cex = 0.3,
                                  mean_col = "black", mean_lty = "dashed",
                                  ln = NULL, rlab = NULL, style = "boxfan", space = 0.5, add_legend = FALSE,
                                  probs_legend = c(1, 25, 50, 75, 99)/100, ylim = NULL, main = NULL, xlab = NULL, las = NULL, axes = TRUE,
                                  ...)
{
  if (length(mom$mu) > 100 | mom$has_temporal_strucutre) {
    warning("You are plotting a lot of units. If you want to plot fanplots of over time use fanplot_prediction() or fanplot_stationary() instead.")
  }

  mu <- mom$mu
  sigma2 <- diag(mom$Sigma)
  psi_ci <- mu/(sigma2/mu - 1)
  matr_cond <- matrix(NA, ncol = length(mu), nrow = length(probs))
  for (i in 1:length(probs)) {
    matr_cond[i, ] <- qnbinom(probs[i], mu = mu, size = psi_ci)
  }
  if (is.null(ylim)) {
    ylim <- c(0.8, 1.1) * range(matr_cond)
  }
  plot(NULL, xlim = c(0, length(mu) + 1 + add_legend), ylim = ylim,
       axes = FALSE, xlab = xlab, ylab = "No. infected", main = main)
  if(axes){
    axis(1, at = 1:length(mu), labels = rownames(mom$mu), las = las)
    axis(2, las = las)
  }
  box
  myx <- seq(0, max(ylim), 500)
  for(i in 1:length(myx))
    lines(c(-1,length(mu) + 1), rep(myx[i], 2), col="grey", lty=2)
  fan(matr_cond, start = 1, frequency = 1, fan.col = fan.col,
      ln = ln, rlab = rlab, style = style, space = space, data.type = "values",
      probs = probs, ...)
  if(add_pred_means){
    for(i in 1:length(mu)){
      lines(i + c(-0.25, 0.25), rep(mu[i], 2), col = mean_col, lty = mean_lty)
    }
  }
  if (add_observed & "predictive_moments_hhh4" %in% class(mom)) {
    points(1:length(mu), mom$realizations, pch = 19, cex = pt.cex,
           col = pt.col)
  }
  if (add_legend) {
    y_legend <- matrix(seq(from = 1.06 * min(ylim), to = max(ylim),
                           length.out = length(probs)), ncol = 1)
    fan(y_legend, start = length(mu) + 1.5, fan.col = fan.col,
        rlab = probs_legend)
    abline(v = length(mu) + 1)
  }
}
