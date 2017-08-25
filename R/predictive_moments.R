predictive_moments <- function(hhh4Obj, t_condition, lgt,
                               return_Sigma = FALSE, return_cov_array = FALSE,
                               return_mu_decomposed = FALSE, return_M = FALSE){
  nu_lambda <- lambda_tilde_complex_neighbourhood(hhh4Obj, subset = seq(from = t_condition + 1, length.out = lgt))

  psi <- if(hhh4Obj$control$family == "Poisson"){
    0
  }else{1/surveillance:::sizeHHH(hhh4Obj$coefficients, terms(hhh4Obj), subset = 1)}

  # calculate M in array form:
  M <- predictive_moments0(nu = nu_lambda$nu, phi = nu_lambda$lambda, psi = psi, stsObj = hhh4Obj$stsObj, t_condition = t_condition, lgt = lgt,
                           return_M = return_M, return_mu_split = return_mu_split)
  # format return object:
  format_return_prediction(M = M, nu = nu_lambda$nu, phi = nu_lambda$lambda,
                           stsObj = hhh4Obj$stsObj, t_condition = t_condition,
                           return_Sigma = return_Sigma, return_cov_array = return_cov_array,
                           return_mu_decomposed = return_mu_decomposed, return_M = return_M)
}



predictive_moments0 <- function(nu, phi, psi, stsObj, t_condition, lgt, M_start = NULL,
                                return_M = FALSE, return_mu_split = FALSE, return_cov_matr = FALSE){

  if(is.matrix(phi) & is.vector(nu)){
    phi <- array(phi, dim = c(nrow(phi), ncol(phi), 1))
    nu <- matrix(nu, nrow = 1)
  }

  n_units <- ncol(nu)
  n_lags <- dim(phi)[2]/n_units
  lgt_tilde <- lgt + n_lags

  M <- array(NA, dim = c(1 + n_units*n_lags, 1 + n_units*n_lags, lgt))
  M_0 <- c(1, as.vector(stsObj@observed[t_condition + 1 - n_lags:1, ]))%*%t(c(1, as.vector(stsObj@observed[t_condition + 1 - n_lags:1, ])))
  M[,,1] <- recursion_M(nu[1, ], phi[,,1], psi, M_0)

  inds_new <- seq(to = ncol(M), length.out = n_units)

  for(i in 2:lgt){
    M[,,i] <- recursion_M(nu[i, ], phi[,,i], psi, M[,,i - 1])
  }

  return(M)
}



# format return object:
format_return_prediction <- function(M, nu, phi, t_condition, stsObj,
                                     return_Sigma, return_cov_array,
                                     return_mu_decomposed, return_M){
  n_units <- ncol(nu)
  n_lags <- dim(phi)[2]/n_units
  lgt <- dim(M)[3]
  names_units <- colnames(nu)

  ret <- list(mu_matrix = NULL, var_matrix = NULL, cov_array = NULL,
              mu_vector = NULL, Sigma = NULL,
              mu_decomposed = NULL, M = NULL)

  # in matrix form
  ret$mu_matrix <- matrix(nrow = lgt, ncol = n_units)
  ret$var_matrix <- matrix(nrow = lgt, ncol = n_units)

  inds <- seq(to = dim(M)[2], length.out = n_units) # place where means/2nd moments of t are found in M
  for(i in 1:lgt){
    ret$mu_matrix[i, ] <- M[1, inds, i]
    mu2_temp <- diag(M[,,i])[inds]
    ret$var_matrix[i, ] <- mu2_temp - ret$mu_matrix[i, ]^2
  }
  # naming:
  colnames(ret$mu_matrix) <- colnames(ret$var_matrix) <- names_units
  rownames(ret$mu_matrix) <- rownames(ret$var_matrix) <- paste0("t=", seq(t_condition + 1, length.out = lgt))

  # array containing covariances if needed:
  if(return_cov_array){
    cov_array0 <- array(dim = c(n_units, n_units, lgt))
    for(i in 1:lgt){
      cov_array0[,,i] <- M[inds, inds, i] - M[1, inds, i] %*% t(M[1, inds, i])
    }
    dimnames(ret$cov_array) <- list(names_units, names_units, rownames(ret$mu_matrix))
  }
  rm(inds)

  # in vector form, if return_Sigma == TRUE
  if(return_Sigma){
    # get matrix of un-centered second moments:
    large_M <- extend_M(ana_mom = M, nu = nu, phi = phi, n_units = n_units,
                        start = 1, n_timepoints = lgt) # function is originally written for periodic case, but this works
    # calculate mean vector and covariance matrix from this:
    ret$mu_vector <- large_M[1, -1]
    ret$Sigma <- large_M[-1, -1] - large_M[1, -1] %*%t(large_M[1, -1])

    names_Sigma <- paste0(rep(names_units, lgt), ";t=",
                          rep(t_condition + seq_len(lgt), each = n_units))
    names(ret$mu_vector) <- colnames(ret$Sigma) <- rownames(ret$Sigma) <- names_Sigma
  }

  if(return_mu_decomposed){# split up by component:
    mu_decomposed0 <- array(dim = c(dim(mu_matrix0), 3))
    dimnames(mu_decomposed0)[[3]] <- c("endemic", "epi.own", "epi.neighbours")
    mu_decomposed0[,,"endemic"] <- nu
    for(i in 1:length_of_period){
      mu_decomposed0[i,,"epi.own"] <-
        only_ar(phi[,,i])%*%as.vector(t(mu_matrix0[i - n_lags:1, ]))
    }
    mu_decomposed0[,,"epi.neighbours"] <- mu_matrix0 - mu_decomposed0[,,"endemic"] - mu_decomposed0[,,"epi.own"]
  }

  if(return_M){
    ret$M <- M
  }

  ret$start <- stsObj@start
  ret$frequency <- stsObj@freq
  ret$t_condition <- t_condition
  ret$condition <- stsObj@observed[t_condition + 1 - n_lags:1, , drop = FALSE]
  ret$n_units <- n_units
  ret$timepoints <- t_condition+ 1:lgt
  ret$realizations_matrix <- stsObj@observed[t_condition + 1:lgt, , drop = FALSE]
  ret$type <- "predictive"
  ret$has_temporal_strucutre <- TRUE

  class(ret) <- c("stationary_moments_hhh4", "moments_hhh4")

  return(ret)
}
