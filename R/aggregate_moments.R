#' Aggregation of stationary or predictive moments
#' 
#' Aggregation of stationary or predictive moments as calculated from \code{hhh4} objects using \code{stationary_moments} or \code{predictive_moments}
#' 
#' @param momentsObj an object of class \code{moments_hhh4} containing stationary or predictive moments, as returned by \code{stationary_moments} or \code{predictive_moments}
#' @param aggregation_matrix an aggregation matrix with either \code{momentsObj$n_units} columns
#' (for aggregation across units while keeping the temporal structure; set option \code{by_timepoint = TRUE} in this case) 
#' or \code{length(momentsObj$mu_vector)} (for aggregation that does not preserve the temporal structure; set option \code{by_timepoint = FALSE}).
#' @param by_timepoint Is aggregation only across units while preserving the temporal structure?
#' @export
aggregate_moments <- function(momentsObj, aggregation_matrix, by_timepoint = FALSE){
  
  # Sigma needs to be available in momentsObj
  if(is.null(momentsObj$Sigma)){
    "Aggregation requires that momentsObj$Sigma is available. Re-run calculation of momentsObj with return_Sigma = TRUE."
  }
  
  n_units <- ncol(momentsObj$mu_matrix)
  n_timepoints <- nrow(momentsObj$mu_matrix)
  if(is.null(rownames(aggregation_matrix))){
    rownames(aggregation_matrix) <- paste0("V", nrow(aggregation_matrix))
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
  ret$mu <- aggregation_matrix %*% momentsObj$mu_vector
  ret$Sigma <- aggregation_matrix%*%momentsObj$Sigma%*%t(aggregation_matrix)
  # names(ret$mu) <- colnames(ret$Sigma) <- rownames(ret$Sigma) <- rownames(aggregation_matrix)
  
  # some quantities can only be aggregated for by_timepoint
  if(by_timepoint){
    ret$mu_matrix <- matrix(ret$mu, nrow = n_timepoints, byrow = TRUE)
    ret$var_matrix <- matrix(diag(ret$Sigma), nrow = n_timepoints, byrow = TRUE)
    colnames(ret$mu_matrix) <- colnames(ret$var_matrix) <- names_new_units
    rownames(ret$mu_matrix) <- rownames(ret$var_matrix) <- rownames(momentsObj$mu_matrix)
    # realizations are only available in prediction case
    if(!is.null(momentsObj$realizations)){
      ret$realizations <- aggregation_matrix %*% momentsObj$realizations
      names(ret$realizations) <- names(ret$mu)
      rownames(ret$mu_matrix) <- rownames(momentsObj$mu_matrix)
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
    # M is not always availbable
    if(!is.null(momentsObj$M)){
      ret$M <- NULL
      warning("Aggregation of M currently not implemented.")
    }
    
    ret$start <- momentsObj$start
    ret$freq <- momentsObj$freq
    ret$n_seasons <- momentsObj$n_seasons
    ret$timepoints <- momentsObj$timepoints
    ret$condition <- momentsObj$condition
  }
  
  ret$type <- momentsObj$type
  ret$has_temporal_strucutre <- by_timepoint
  
  class(ret) <- class(momentsObj)

  return(ret)
}