#' Stationary distributions implied by a \code{hhh4}-model
#'
#' Returns the mean vector and covariance matrix of the stationary distribution implied by an hhh4 object.
#'
#' @param hhh4Obj \code{hhh4} object for which to calculate stationary moments
#' @param start start of the season
#' @param n_seasons number of
#' @param return_Sigma logical: return entire variance/covariance matrix of the prediction; can take a lot of storage
#' @param return_cov_array logical: return an array containing week-wise covariance matrices
#' @param return_mu_decomposed logical: return an array containing a decomposition of stationary means into the three
#' components \code{endemic}, \code{epi.own} and \code{epi.others}.
#' @param return_M logical: return the array M containing un-centered second moments (used internally for calculations)
#' @param max.iter maximum number of iterations before iterative algorithm stops
#' @param tolerance element-wise maximum tolerance (entering into temination criterion for the iterative calculation)
#'
#' @export
stationary_moments <- function(hhh4Obj, start = 1, n_seasons = 1,
                               return_Sigma = FALSE, return_cov_array = FALSE,
                               return_mu_decomposed = FALSE, return_M = FALSE,
                               max.iter = 10, tolerance = 1e-5){
  if(!("hhh4" %in% class(hhh4Obj))){
    stop("hhh4Obj neds to be of class hhh4 or an extension thereof.")
  }
  nu_lambda <- lambda_tilde_complex_neighbourhood(hhh4Obj, periodic = TRUE)

  psi <- if(hhh4Obj$control$family == "Poisson"){
    0
  }else{1/surveillance:::sizeHHH(hhh4Obj$coefficients, terms(hhh4Obj), subset = 1)}

  # calculate M in array form:
  M <- stationary_moments0(nu = nu_lambda$nu, phi = nu_lambda$lambda, psi = psi,
                           start = start, n_seasons = n_seasons,
                           return_Sigma = return_Sigma, return_cov_array = return_cov_array,
                           return_mu_decomposed = return_mu_decomposed, return_M = return_M,
                           max.iter = max.iter, tolerance = tolerance)

  # create, format and fill return objects:
  format_return(M = M, nu = nu_lambda$nu, phi = nu_lambda$lambda, n_seasons = n_seasons, start = start,
                return_Sigma = return_Sigma, return_cov_array = return_cov_array,
                return_mu_decomposed = return_mu_decomposed, return_M = return_M)
}

stationary_moments0 <- function(nu, phi, psi, start = 1, n_seasons = 1,
                                return_Sigma = FALSE, return_cov_array = FALSE,
                                return_mu_decomposed = FALSE, return_M = FALSE,
                                max.iter = 10, tolerance = 1e-5){

  # transform nu, phi to matrix and array if necessary
  if(is.matrix(phi) & is.vector(nu)){
    phi <- array(phi, dim = c(nrow(phi), ncol(phi), 1))
    nu <- matrix(nu, nrow = 1)
  }

  # extract dimensions:
  n_units <- ncol(nu)
  n_lags <- dim(phi)[2]/n_units
  length_of_period <- dim(phi)[3]
  names_units <- colnames(nu)

  # initialize matrix M contatining first and second moments:
  M <- array(1, dim = c(1 + n_units*n_lags, 1 + n_units*n_lags, length_of_period))
  names_M <- c("1", paste0(rep(names_units, n_lags), ";t-",
                           rep((n_lags - 1):0, each = n_units)))
  dimnames(M) <- list(names_M, names_M, paste0("t=", 1:length_of_period))

  # re-order if necessary:
  re_order <- c(start:length_of_period, if(start != 1) 1:(start - 1))
  M <- M[,,re_order, drop = FALSE]
  nu <- nu[re_order, , drop = FALSE]
  phi <- phi[,,re_order, drop = FALSE]

  # iterative calculation of M:
  for(j in 1:max.iter){
    temp_1 <- M[,,1]
    M[,,1] <- recursion_M(nu[1, ], phi[,,1],
                          psi, M[,,length_of_period])
    # stop when convergence is reached:
    # if(all.equal(M[,,1], temp_1, tolerance = tolerance) == TRUE){
    #   break()
    # }
    for(i in 2:length_of_period){
      M[,,i] <- recursion_M(nu[i, ], phi[,,i], psi, M[,,i - 1])
    }
  }

  return(M)
}

# helper function: takes phi_wide and sets all cross-region parameters to 0.
only_ar <- function(matr){
  n_units <- nrow(matr)
  max_lag <- ncol(matr)/n_units
  a <- c(rep(n_units + 1, n_units - 1), 1)
  b <- cumsum(c(1, rep(a, max_lag)))[1:(n_units*max_lag)]
  matr[-b] <- 0
  matr
}


# adding off-blockdiagonals
extend_M <- function(ana_mom, nu, phi, n_units, start, n_timepoints){
  length_of_period <- dim(ana_mom)[3]
  n_lags <- (ncol(ana_mom[,,1]) - 1)/n_units

  # re-order nu and phi if necessary:
  re_order <- c(start:length_of_period, if(start != 1) 1:(start - 1))
  ana_mom <- ana_mom[,,re_order]
  nu <- nu[re_order, ]
  phi <- phi[,,re_order]

  extended_M <- matrix(ncol = 1 + n_timepoints*n_units,
                       nrow = 1 + n_timepoints*n_units)
  extended_M[1, 1] <- 1
  for(i in n_lags:n_timepoints){
    # fill block diagonal:
    inds_blockdiag <- seq(to = i*n_units + 1, length.out = n_units*n_lags)
    ind_ana_mom <- ifelse((start + i - 1)%%length_of_period == 0,
                          length_of_period, (start + i - 1)%%length_of_period)
    extended_M[1, inds_blockdiag] <- extended_M[inds_blockdiag, 1] <- ana_mom[1, -1, ind_ana_mom]
    extended_M[inds_blockdiag, inds_blockdiag] <- ana_mom[-1, -1,ind_ana_mom]
    # fill remaining parts:
    if( i >= n_lags + 1){
      phi_star <- cbind(nu[ind_ana_mom, ], phi[, , ind_ana_mom])
      inds_t <- seq(to = i*n_units + 1, length.out = n_units)
      inds_off_blockdiag <- 2:((i - n_lags)*n_units + 1)
      inds_lags <- c(1, seq(to = (i - 1)*n_units + 1, length.out = n_lags*n_units))
      extended_M[inds_t, inds_off_blockdiag] <- phi_star %*% extended_M[inds_lags, inds_off_blockdiag]
      extended_M[inds_off_blockdiag, inds_t] <- t(extended_M[inds_t, inds_off_blockdiag])
      # extended_M[inds_off_blockdiag, inds_t] <- extended_M[inds_off_blockdiag, inds_lags]%*%t(phi_star)
      # extended_M[inds_t, inds_off_blockdiag] <- t(extended_M[inds_off_blockdiag, inds_t])
    }
  }
  extended_M
}


# format return object:
format_return <- function(M, nu, phi, n_seasons, start,
                          return_Sigma, return_cov_array,
                          return_mu_decomposed, return_M){
  n_units <- ncol(nu)
  n_lags <- dim(phi)[2]/n_units
  length_of_period <- dim(phi)[3]
  names_units <- colnames(nu)

  # little helper function to handle periodicity:
  ominus <- function(a, b, l){
    ifelse(a > b, a - b, l - ((b - a)%%l))
  }

  ret <- list(mu_matrix = NULL, var_matrix = NULL, cov_array = NULL,
              mu_vector = NULL, Sigma = NULL,
              mu_decomposed = NULL, M = NULL)

  # in matrix form
  mu_matrix0 <- matrix(nrow = length_of_period, ncol = n_units)
  var_matrix0 <- matrix(nrow = length_of_period, ncol = n_units)

  inds <- seq(to = dim(M)[2], length.out = n_units) # place where means/2nd moments of t are found in M
  for(i in 1:length_of_period){
    mu_matrix0[i, ] <- M[1, inds, i]
    mu2_temp <- diag(M[,,i])[inds]
    var_matrix0[i, ] <- mu2_temp - mu_matrix0[i, ]^2
  }
  # append n_season times:
  ret$mu_matrix <- matrix(rep(t(mu_matrix0), n_seasons), ncol = n_units, byrow = TRUE)
  ret$var_matrix <- matrix(rep(t(var_matrix0), n_seasons), ncol = n_units, byrow = TRUE)
  # naming:
  colnames(ret$mu_matrix) <- colnames(ret$var_matrix) <- names_units
  rownames(ret$mu_matrix) <- rownames(ret$var_matrix) <- paste0("t=", seq(start, length.out = n_seasons*length_of_period))

  # array containing covariances if needed:
  if(return_cov_array){
    cov_array0 <- array(dim = c(n_units, n_units, length_of_period))
    for(i in 1:length_of_period){
      cov_array0[,,i] <- M[inds, inds, i] - M[1, inds, i] %*% t(M[1, inds, i])
    }
    # append n_seasons times:
    ret$cov_array <- array(rep(cov_array0, n_seasons), dim = dim(cov_array0)*c(1, 1, n_seasons))
    dimnames(ret$cov_array) <- list(names_units, names_units, rownames(ret$mu_matrix))
  }
  rm(inds)

  # in vector form, if return_Sigma == TRUE
  if(return_Sigma){
    # get matrix of un-centered second moments:
    large_M <- extend_M(ana_mom = M, nu = nu, phi = phi, n_units = n_units,
                        start = start, n_timepoints = n_seasons*length_of_period)
    # calculate mean vector and covariance matrix from this:
    ret$mu_vector <- large_M[1, -1]
    ret$Sigma <- large_M[-1, -1] - large_M[1, -1] %*%t(large_M[1, -1])

    names_Sigma <- paste0(rep(names_units, n_seasons*length_of_period), ";t=",
                          rep(start + seq_len(n_seasons*length_of_period) - 1, each = n_units))
    names(ret$mu_vector) <- colnames(ret$Sigma) <- rownames(ret$Sigma) <- names_Sigma
  }

  if(return_mu_decomposed){# split up by component:
    mu_decomposed0 <- array(dim = c(dim(mu_matrix0), 3))
    dimnames(mu_decomposed0)[[3]] <- c("endemic", "epi.own", "epi.neighbours")
    mu_decomposed0[,,"endemic"] <- nu
    for(i in 1:length_of_period){
      mu_decomposed0[i,,"epi.own"] <-
        only_ar(phi[,,i])%*%as.vector(t(mu_matrix0[ominus(i, n_lags:1, length_of_period), ]))
    }
    mu_decomposed0[,,"epi.neighbours"] <- mu_matrix0 - mu_decomposed0[,,"endemic"] - mu_decomposed0[,,"epi.own"]
    # append n_seasons times:
    if(n_seasons > 1){
      ret$mu_decomposed <- array(dim = dim(mu_decomposed0)*c(n_seasons, 1, 1))
      for(i in 1:n_seasons){
        inds <- seq(to = i*length_of_period, length.out = length_of_period)
        ret$mu_decomposed[inds,,] <- mu_decomposed0
        dimnames(ret$mu_decomposed) <- list(rownames(ret$mu_matrix), colnames(ret$mu_matrix),
                                            c("endemic", "epi.own", "epi.neighbours"))
      }
    }else{
      ret$mu_decomposed <- mu_decomposed0
    }
  }

  if(return_M){
    ret$M <- M
  }

  ret$start <- start
  ret$freq <- length_of_period
  ret$n_seasons <- n_seasons
  ret$n_units <- n_units
  ret$timepoints <- seq(from = start, length.out = n_seasons*length_of_period)
  ret$condition <- NULL
  ret$type <- "stationary"
  ret$has_temporal_strucutre <- TRUE

  class(ret) <- c("stationary_moments_hhh4", "moments_hhh4")

  return(ret)
}
