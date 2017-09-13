#' Calculate Dawid-Sebastiani score
#'
#' Calculate Dawid-Sebastiani score for a prediction returned by \code{predictive_moments}.
#'
#' The Dawid-Sebastiani score is defined as
#' \deqn{DSS = log(|\Sigma|) + t(Y_obs - \mu) \Sigma^{-1} (Y_obs - \mu)}
#' where \eqn{\mu} and \eqn{\Sigma} are the predictive mean and variance, repectively.
#' Y_obs represents the observation that has materialized.
#'
#' @param pred the prediction as returned by \code{longterm_prediction_hhh4}
#' (and potentially aggregated using aggregate_prediction)
#' @param detailed detailed or less detailed output?
#' @param scaled if \code{detailed == FALSE}: scale DSS with 2d?
#'
#' @return If \code{detailed == FALSE}: the (potentially scaled) Dawid-Sebastiani score. If
#' \code{detailed == TRUE}: a vector containing the following elements:
#' \itemize{
#' \item \code{dawid_sebastiani} the un-scaled Dawid-Sebastiani score
#' \item \code{term1} value of the log-determinant entering into the unn-scaled Dawid-Sebastiani score
#' \item \code{term2} value of the quadratic form entering into the un-scaled Dawid-Sebastiani score
#' \item \code{scaled_dawid_sebastiani} the scaled Dawid-Sebastiani score
#' \item \code{determinant_sharpness} the determinant sharpness (scaled version of \code{term1})
#' }
#'
#' @examples
#' ## a simple univariate example:
#' data("salmonella.agona")
#' ## convert old "disProg" to new "sts" data class
#' salmonella <- disProg2sts(salmonella.agona)
#' # specify and fit model: fixed geometric lag structure
#' # with weight 0.8 for first lag
#' control_salmonella <- list(end = list(f = addSeason2formula(~ 1)),
#'                            ar = list(f = addSeason2formula(~ 1),
#'                                      par_lag = 0.8, use_distr_lag = TRUE),
#'                            family = "NegBinM", subset = 6:312)
#' fit_salmonella <- hhh4_lag(salmonella, control_salmonella)
#' pred_salmonella <- predictive_moments(fit_salmonella, t_condition = 260,
#'                                       52, return_Sigma = TRUE)
#' ds_score_hhh4(pred_salmonella, detailed = TRUE)
#'
#' @export
ds_score_hhh4 <- function(pred, detailed = FALSE, scaled = TRUE){
  # Sigma needs to be available in momentsObj
  if(is.null(pred$Sigma)){
    stop("Aggregation requires that momentsObj$Sigma is available. Re-run calculation of momentsObj with return_Sigma = TRUE.")
  }
  dim <- length(pred$mu_vector)
  # calculate the different terms which are involved:
  term1 <- determinant(pred$Sigma, logarithm = TRUE)$modulus[1]
  term2 <- t(pred$realizations - pred$mu_vector) %*% solve(pred$Sigma) %*% (pred$realizations - pred$mu_vector)
  determinant_sharpness <- term1/(2*dim)
  dawid_sebastiani <- term1 + term2
  scaled_dawid_sebasiani <- dawid_sebastiani/(2*dim)
  # return detailed or simple result:
  if(detailed){
    return(c(dawid_sebastiani = dawid_sebastiani,
             term1 = term1,
             term2 = term2,
             scaled_dawid_sebastiani = scaled_dawid_sebasiani,
             determinant_sharpness = term1/(2*dim)))
  }else{
    return(ifelse(scaled, scaled_dawid_sebasiani, dawid_sebastiani))
  }
}
