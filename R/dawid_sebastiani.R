#' Calculate Dawid-Sebastiani score
#'
#' Calculate Dawid-Sebastiani score for a prediction returned by \code{longterm_prediction_hhh4}
#'
#' @param pred the prediction as returned by \code{longterm_prediction_hhh4}
#' (and potentially aggregated using aggregate_prediction)
#' @param detailed detailed (can be very large!) or less detailed output?
#' @param scaled if \code{detailed == FALSE}: scale DSS with 2d?
#'
#' @export
ds_score_hhh4 <- function(pred, detailed = FALSE, scaled = TRUE){
  dim <- length(pred$mu_vector)
  term1 <- determinant(pred$Sigma, logarithm = TRUE)$modulus[1]
  term2 <- t(pred$realizations - pred$mu_vector) %*% solve(pred$Sigma) %*% (pred$realizations - pred$mu_vector)
  determinant_sharpness <- term1/(2*dim)
  dawid_sebastiani <- term1 + term2
  scaled_dawid_sebasiani <- dawid_sebastiani/(2*dim)
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
