#' Aggregation of stationary or predictive moments
#'
#' Aggregation of stationary or predictive moments as calculated using
#' \code{stationary_moments} or \code{predictive_moments}.
#'
#' @param momentsObj an object of class \code{moments_hhh4} containing stationary
#'  or predictive moments, as returned by \code{stationary_moments} or
#'  \code{predictive_moments}
#' @param aggregation_matrix an aggregation matrix with either
#' \code{momentsObj$n_units} columns
#' (for aggregation across units while keeping the temporal structure; set option
#' \code{by_timepoint = TRUE} in this case)
#' or \code{length(momentsObj$mu_vector)} (for aggregation that does not preserve
#' the temporal structure; set option \code{by_timepoint = FALSE}).
#' @param by_timepoint logical: is aggregation only across units while preserving
#' the temporal structure? Note that the new  \code{moments_hhh4} object
#' cannot have the \code{condition} , \code{mu_matrix}, \code{var_matrix} and
#' \code{cov_array} elements if the temporal structure is given up.
#' @return An object of class \code{moments_hhh4} representing the new prediction.
#'
#' @examples
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
#' # plot predictions for two regions:
#' par(mfrow = 1:2)
#' fanplot_prediction(pred_mom, unit = 1, main = "Bremen")
#' fanplot_prediction(pred_mom, unit = 2, main = "Lower Saxony")
#'
#' #########
#' # aggregation 1: combine the two regions
#' aggr_matr_pool <- matrix(1, ncol = 2)
#' # specify by_timepoint = TRUE to keep the temporal structure and aggregate only
#' # counts from the same week:
#' pred_mom_pooled <- aggregate_moments(pred_mom, aggr_matr_pool, by_timepoint = TRUE)
#' fanplot_prediction(pred_mom_pooled, unit = 1, ylim = c(0, 500), main = "Aggregation over regions")
#'
#' #########
#' # aggregation 2: total burden in the two regions
#' aggr_matr_total_burden <- matrix(rep(c(1, 0, 0, 1), 52), nrow = 2,
#'                                  dimnames = list(c("Bremen", "Lower Saxony"),
#'                                                  NULL))
#' pred_mom_total_burden <- aggregate_moments(pred_mom, aggr_matr_total_burden)
#' plot_moments_by_unit(pred_mom_total_burden, main = "Total burdens")
#'
#' #########
#' # works also with stationary moments:
#' stat_mom <- stationary_moments(fitBL, return_Sigma = TRUE)
#' stat_mom_pooled <- aggregate_moments(stat_mom, aggr_matr_pool, by_timepoint = TRUE)
#' stat_mom_total_burden <- aggregate_moments(stat_mom, aggr_matr_total_burden, by_timepoint = FALSE)
#' fanplot_stationary(stat_mom_pooled)
#' plot_moments_by_unit(stat_mom_total_burden, main = "Total burdens")
#'
#' @export
aggregate_moments <- function(momentsObj, aggregation_matrix, by_timepoint = FALSE){

  # Sigma needs to be available in momentsObj
  if(is.null(momentsObj$Sigma)){
    stop("Aggregation requires that momentsObj$Sigma is available. Re-run calculation of momentsObj with return_Sigma = TRUE.")
  }

  n_units <- ncol(momentsObj$mu_matrix)
  n_timepoints <- nrow(momentsObj$mu_matrix)
  if(is.null(rownames(aggregation_matrix))){
    rownames(aggregation_matrix) <- paste0("V", 1:nrow(aggregation_matrix))
  }

  # check validity of arguments
  if(by_timepoint & (ncol(aggregation_matrix) != n_units)){
    stop("If by_timepoint == TRUE: ncol(aggregation_matrix) needs to be equal to the number of units in momentsObj.")
  }
  if(!by_timepoint & (ncol(aggregation_matrix) != ncol(momentsObj$Sigma))){
    stop("If by_timepoint == FALSE: ncol(aggregation_matrix) needs to be equal to ncol(momentsObj$Sigma).")
  }

  # if by_timepoint == TRUE: construct large transformation matrix for all timepoints combined
  if(by_timepoint){
    n_units_new <- nrow(aggregation_matrix)
    names_new_units <- rownames(aggregation_matrix)

    aggregation_matrix0 <- aggregation_matrix

    aggregation_matrix <- matrix(0, ncol = n_units*n_timepoints, nrow = n_units_new*n_timepoints)
    rownames(aggregation_matrix) <- paste0(rep(names_new_units, n_timepoints), ";t=",
                                           rep(momentsObj$timepoints, each = n_units_new))
    colnames(aggregation_matrix) <- names(momentsObj$mu)

    for(i in 1:n_timepoints){
      inds_rows <- seq(to = n_units_new*i, length.out = n_units_new)
      inds_cols <- seq(to = n_units*i, length.out = n_units)
      aggregation_matrix[inds_rows, inds_cols] <- aggregation_matrix0
    }
  }

  # initialize return object:
  ret <- list()

  # actual aggregation:
  ret$mu_vector <- aggregation_matrix %*% momentsObj$mu_vector
  ret$Sigma <- aggregation_matrix%*%momentsObj$Sigma%*%t(aggregation_matrix)
  # names(ret$mu_vector) <- colnames(ret$Sigma) <- rownames(ret$Sigma) <- rownames(aggregation_matrix)

  # M is not always availbable
  if(!is.null(momentsObj$M)){
    ret$M <- NULL
    warning("Aggregation of M currently not implemented.")
  }

  # realizations are not always available
  if(!is.null(momentsObj$realizations)){
    ret$realizations <- aggregation_matrix %*% as.vector(t(momentsObj$realizations))
    names(ret$realizations) <- names(ret$mu_vector)
  }

  # some quantities can only be aggregated for by_timepoint
  if(by_timepoint){
    ret$mu_matrix <- matrix(ret$mu_vector, nrow = n_timepoints, byrow = TRUE)
    ret$var_matrix <- matrix(diag(ret$Sigma), nrow = n_timepoints, byrow = TRUE)
    colnames(ret$mu_matrix) <- colnames(ret$var_matrix) <- names_new_units
    rownames(ret$mu_matrix) <- rownames(ret$var_matrix) <- rownames(momentsObj$mu_matrix)
    # realizations are only available in prediction case
    if(!is.null(momentsObj$realizations)){
      ret$realizations_matrix <- matrix(ret$realizations, nrow = n_timepoints, byrow = TRUE)
      dimnames(ret$realizations_matrix) <- dimnames(ret$mu_matrix)
    }
    # cov_array is not always available
    if(!is.null(momentsObj$cov_array)){
      ret$cov_array <- array(dim = c(n_units_new, n_units_new, n_timepoints))
      dimnames(ret$cov_array) <- (dimnames(ret$mu_matrix)[c(2, 2, 1)])
      for(i in 1:n_timepoints){
        inds_Sigma <- seq(to = i*n_units_new, length.out = n_units_new)
        ret$cov_array[,,i] <- ret$Sigma[inds_Sigma, inds_Sigma]
      }
    }

    ret$start <- momentsObj$start
    ret$freq <- momentsObj$freq
    ret$timepoints <- momentsObj$timepoints


    # certain timing informationd and conditions can only be handled for predictive moments:
    if("predictive_moments_hhh4" %in% class(momentsObj)){
      ret$n_seasons <- momentsObj$n_seasons
      ret$timepoints_calendar <- momentsObj$timepoints_calendar
      ret$t_condition <- momentsObj$t_condition
      ret$condition <- aggregation_matrix0 %*% t(momentsObj$condition)
    }
  }

  ret$type <- momentsObj$type
  ret$has_temporal_structure <- by_timepoint

  class(ret) <- class(momentsObj)

  return(ret)
}
