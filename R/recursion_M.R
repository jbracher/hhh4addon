recursion_M <- function(nu_t, phi_t, psi, last_M){
  
  # extract features:
  n_units <- length(nu_t)
  n_lags <- if(is.matrix(phi_t)) ncol(phi_t)/n_units else 1 # handle univariate case with one lag
  inds_new <- seq(to = ncol(last_M), length.out = n_units)
  
  # extend phi to phi_tilde
  phi_tilde <- rbind(c(1, rep(0, n_units*n_lags)),
                     cbind(matrix(0, nrow = (n_lags - 1)*n_units, ncol = n_units + 1),
                           diag((n_lags - 1)*n_units)),
                     cbind(nu_t, phi_t)
  )

  # caclucation, step 1: quadratic form
  next_M <- phi_tilde%*%last_M%*%t(phi_tilde)
  # step 2: manipulation of the diagonal
  diag(next_M)[inds_new] <- (1 + psi)*diag(next_M)[inds_new] + next_M[1, inds_new]
  # step 3: setting first diagonal element back to 1:
  next_M[1, 1] <- 1
  # return
  next_M
}