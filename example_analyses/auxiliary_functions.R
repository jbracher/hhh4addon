# scp johannes@130.60.71.234:/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis/auxiliary_functions.R auxiliary_functions.R

csv_to_sts <- function(file, names, start, end, ...){
  # read data:
  dat <- read.csv2(file,
                   sep = ",", header = TRUE, stringsAsFactors = FALSE)
  dat$year <- dat$week <- NA
  # handle week 53:
  is_week53 <- which(grepl("w53", dat$time))
  dat <- dat[-is_week53, ]
  # handle time variable:
  for(i in 1:nrow(dat)){
    temp <- as.numeric(strsplit(dat$time[i], "w", fixed = TRUE)[[1]])
    dat$year[i] <- temp[1]
    dat$week[i] <- temp[2]
  }
  dat$time <- NULL
  colnames(dat)[1:length(names)] <- names
  # restrict to selected range
  if( (tail(dat$year, 1) < end[1]) ||
      (tail(dat$year, 1) == end[1]) & tail(dat$week, 1) < end[2]){
    stop("Either start or end is outside of the range of the provided data.")
  }

  dat <- subset(dat, year >= start[1])
  dat <- subset(dat, !(year == start[1] & week < start[2]))

  dat <- subset(dat, year <= end[1])
  dat <- subset(dat, !(year == end[1] & week > end[2]))

  dat$week <- NULL
  dat$year <- NULL

  # to sts object:
  stsObj <- new("sts", observed = dat, start = start,...)
  return(stsObj)
}

# obtain some quantities from the fits: endemic components:
get_nu_seas <- function(hhh4Obj, regions, christmas){
  pars <- hhh4Obj$coefficients
  nu <- matrix(ncol = 52, nrow = length(regions), dimnames = list(regions))
  for(region in regions){
    nu[region, ] <- exp(pars[paste0("end.1.", region)] +
                          pars[paste0("end.sin(2 * pi * t/52)")]*sin(2*pi/52*(1:52)) +
                          pars[paste0("end.cos(2 * pi * t/52)")]*cos(2*pi/52*(1:52)) +
                          pars[paste0("end.christmas")] * christmas)
  }

  return(nu)
}

# for comparison: empirical moments:
get_emp_moments <- function(sts, start = 1){
  obs <- sts@observed
  start_week <- sts@start[2]
  freq <- sts@freq
  end_week <- (start_week + nrow(obs) - 1) %% freq

  obs <- rbind(matrix(nrow = start_week - 1, ncol = ncol(obs)),
               obs,
               matrix(nrow = freq - end_week, ncol = ncol(obs)))
  weekwise_means <- weekwise_vars <- matrix(ncol = ncol(obs), nrow = freq)
  for(i in 1:ncol(sts)){
    matr_temp <- matrix(obs[, i], ncol = freq, byrow = TRUE)
    weekwise_means[, i] <- colMeans(matr_temp, na.rm = TRUE)
    weekwise_vars[, i] <- apply(matr_temp, 2, var, na.rm = TRUE)
  }
  colnames(weekwise_means) <- colnames(weekwise_vars) <- colnames(obs)
  reorder <- c(seq(from = start, to = freq),
               seq(to = start - 1, length.out = start - 1))
  weekwise_means <- weekwise_means[reorder, , drop = FALSE]
  weekwise_vars <- weekwise_vars[reorder, , drop = FALSE]
  return(list(mean = weekwise_means, var = weekwise_vars))
}

### THIS VERSION DOES NOT HANDLE THE END OF THE YEAR CORRECTLY!!!
# # compute the stationary moments of calendar week-wise averages
# compute_sm_of_means <- function(fit, n_seasons = 7, start = 1, return_Sigma = FALSE){
#   # obtain dimensions:
#   n_units <- ncol(fit$stsObj@observed)
#   freq <- fit$stsObj@freq
#   # get stationary moments for n_seasons seasons:
#   sm <- stationary_moments(fit, start = start, n_seasons = 3, return_Sigma = TRUE)
#   inds_first <- seq(from = 1, length.out = n_units*freq)
#   inds_last <- seq(to = 3*n_units*freq, length.out = n_units*freq)
#   inds_middle <- seq(to = 2*n_units*freq, length.out = n_units*freq)
#   # we avoid spanning up the matrix for the full n_seasons seasons by just computing
#   # everything for 3 seasons and re-using the middle part of the matrix n_seasons - 2 times.
#
#   new_mu_vector <- sm$mu_vector[inds_first]
#   new_Sigma <- (sm$Sigma[inds_first, inds_first] +
#                   (n_seasons - 2)*sm$Sigma[inds_middle, inds_middle] +
#                   sm$Sigma[inds_last, inds_last])/n_seasons^2
#
#   new_sm <- list()
#
#   if(return_Sigma){
#     new_sm$mu_vector <- new_mu_vector
#     new_sm$Sigma <- new_Sigma
#   }
#
#   new_sm$mu_matrix <- matrix(new_mu_vector, ncol = n_units, byrow = TRUE)
#   new_sm$var_matrix <- matrix(diag(new_Sigma), ncol = n_units, byrow = TRUE)
#
#   rownames(new_sm$mu_matrix) <- rownames(new_sm$var_matrix) <- 1:freq
#   colnames(new_sm$mu_matrix) <- colnames(new_sm$var_matrix) <- colnames(sm$mu_matrix)
#
#   return(new_sm)
# }

#########################################
# compute de-correlated Pearson residuals:
compute_decorr_pearson_residuals <- function(stat_mom_of_means, emp_mom){
  diffs <- (as.vector(t(emp_mom$mean)) - stat_mom_of_means$mu_vector)
  resid_uncorr <- t(chol(solve(stat_mom_of_means$Sigma)))%*%t(t(diffs))
  return(resid_uncorr)
}

# more computing intensive version, only kept for comparison:
compute_sm_of_means <- function(fit, n_seasons, start = 1, return_Sigma = FALSE){
  # obtain dimensions:
  n_units <- ncol(fit$stsObj@observed)
  freq <- fit$stsObj@freq
  # get stationary moments for n_seasons seasons:
  sm <- stationary_moments(fit, start = start, n_seasons = n_seasons, return_Sigma = TRUE)
  # transformation matrix to move to calendar week-wise means:
  trafo_matr <- matrix(0, nrow = n_units*freq, ncol = n_units*freq*n_seasons)
  for(i in 1:n_seasons){
    diag(trafo_matr[, ((i - 1)*(n_units*freq) + 1:(n_units*freq))]) <- 1/n_seasons
  }
  # aggregate moments:
  new_sm <- list()

  new_mu_vector <- trafo_matr%*%sm$mu_vector
  new_Sigma <- trafo_matr%*%sm$Sigma%*%t(trafo_matr)

  if(return_Sigma){
    new_sm$mu_vector <- new_mu_vector
    new_sm$Sigma <- new_Sigma
  }

  new_sm$mu_matrix <- matrix(new_mu_vector, ncol = n_units, byrow = TRUE)
  new_sm$var_matrix <- matrix(diag(new_Sigma), ncol = n_units, byrow = TRUE)

  rownames(new_sm$mu_matrix) <- rownames(new_sm$var_matrix) <- 1:freq
  colnames(new_sm$mu_matrix) <- colnames(new_sm$var_matrix) <- colnames(sm$mu_matrix)

  return(new_sm)
}

# idea: fill only block diagonal and
compute_sm_of_means2 <- function(fit, n_seasons, start = 1, return_Sigma = FALSE){
  # obtain dimensions:
  n_units <- ncol(fit$stsObj@observed)
  freq <- fit$stsObj@freq
  # get stationary moments for n_seasons seasons:
  sm2 <- stationary_moments(fit, start = start, n_seasons = 2, return_Sigma = TRUE)

  lgt_n_seas <- n_units*freq*n_seasons

  long_mu_vector <- rep(sm2$mu_vector, length.out = lgt_n_seas)
  large_Sigma <- matrix(0, nrow = lgt_n_seas, ncol = lgt_n_seas)
  for(i in 1:(n_seasons - 1)){
    large_Sigma[(i - 1)*n_units*freq + (1:(2*n_units*freq)),
          (i - 1)*n_units*freq + (1:(2*n_units*freq))] <-
      sm2$Sigma
  }

  # transformation matrix to move to calendar week-wise means:
  trafo_matr <- matrix(0, nrow = n_units*freq, ncol = n_units*freq*n_seasons)
  for(i in 1:n_seasons){
    diag(trafo_matr[, ((i - 1)*(n_units*freq) + 1:(n_units*freq))]) <- 1/n_seasons
  }
  # aggregate moments:
  new_sm <- list()

  new_mu_vector <- trafo_matr%*%long_mu_vector
  new_Sigma <- trafo_matr%*%large_Sigma%*%t(trafo_matr)

  if(return_Sigma){
    new_sm$mu_vector <- new_mu_vector
    new_sm$Sigma <- new_Sigma
  }

  new_sm$mu_matrix <- matrix(new_mu_vector, ncol = n_units, byrow = TRUE)
  new_sm$var_matrix <- matrix(diag(new_Sigma), ncol = n_units, byrow = TRUE)

  rownames(new_sm$mu_matrix) <- rownames(new_sm$var_matrix) <- 1:freq
  colnames(new_sm$mu_matrix) <- colnames(new_sm$var_matrix) <- colnames(sm$mu_matrix)

  return(new_sm)
}

##################
# p-values for region-wise versions:
get_pval <- function(stat_mom_of_means, emp_mom, correction_df = 0){
  diffs <- as.vector(t(emp_mom$mean)) - stat_mom_of_means$mu_vector
  test_stat <- (t(diffs)%*%solve(stat_mom_of_means$Sigma)%*%diffs)[1, 1]
  df <- length(diffs) - correction_df
  pval <- pchisq(test_stat, df, lower.tail = FALSE)
  return(list(test_stat = test_stat, df = df, pval = pval))
}

# p-values for region-wise versions with 2/3 correction:
get_pval2 <- function(stat_mom_of_means, emp_mom, correction_df = 0){
  h <- function(x) x^(2/3)
  diffs <- h(as.vector(t(emp_mom$mean))) - h(stat_mom_of_means$mu_vector)
  cov_diffs <- diag(2/3*as.vector(stat_mom_of_means$mu_vector)^(-1/3)) %*% stat_mom_of_means$Sigma %*% diag(2/3*as.vector(stat_mom_of_means$mu_vector)^(-1/3))
  test_stat <- t(diffs)%*%solve(cov_diffs)%*%diffs
  df <- length(diffs) - correction_df
  pval <- pchisq(test_stat, df, lower.tail = FALSE)
  return(list(test_stat = test_stat, df = df, pval = pval))
}

### one-step-ahead forecasts also updating the lag decay parameter
# arguments: fit: an hhh4_lag object from which to generate one-step-ahead-forecasts
#             tp: range of timepoints for which to generate forecasts, i.e. tp[1] + 1, ..., tp[2] + 1.
osa_refitting_par_lag <- function(fit, tp){
  # extract info from arguments:
  control <- fit$control
  sts <- fit$stsObj
  start_subset <- control$subset[1]
  tp_vect <- tp[1]:tp[2]

  # initialize osa list to store results:
  templ <- matrix(ncol = ncol(sts@observed), nrow = length(tp_vect),
                  dimnames = list(tp_vect + 1, colnames(sts@observed)))
  osa <- list(pred = templ, observed = templ, psi = templ, par_lag = numeric(length(tp_vect)))
  if(control$family == "Poisson") osa$psi <- NULL

  mod_temp <- fit

  # run through time points for which to obtain one-step-ahead-predictions
  for(i in seq_along(tp_vect)){
    # adapt subset:
    control$subset <- start_subset:(tp_vect[i])
    # adopt grid for lag decay parameter
    # new_grid <- get_grid_around_previous_par_lag(mod_temp)
    # update model:
    mod_temp <- profile_par_lag(sts, control = control, reltol_par_lag = 0.01, start_par_lag = mod_temp$par_lag)
    # get one-step-ahead prediction
    osa_temp <- suppressMessages(oneStepAhead_hhh4lag(mod_temp, tp = c(tp_vect[i], tp_vect[i]), type = "final"))
    # store result:
    osa$pred[i, ] <- osa_temp$pred
    osa$observed[i, ] <- osa_temp$observed
    if(control$family != "Poisson"){
      osa$psi[i, ] <- osa_temp$psi
    }
    osa$par_lag[i] <- mod_temp$par_lag
    print(i)
  }

  return(osa)
}

# Functions for decomposing the fitted values (used in plotting)
decompose_epidemic_component <- function(fit){
  # extract info:
  sts <- fit$stsObj
  max_lag <- if(class(fit)[1] == "hhh4lag") fit$max_lag else 1
  subset <- fit$control$subset
  n_units <- ncol(sts@observed)
  param <- hhh4addon:::lambda_tilde_complex_neighbourhood(fit, periodic = FALSE,
                                                          subset = 1:max(subset))

  # initialize:
  contributions <- array(dim = c(max(subset),
                                 n_units,
                                 n_units,
                                 max_lag),
                         dimnames = list(1:max(subset),
                                         paste0("from.", colnames(sts@observed)),
                                         paste0("to.", colnames(sts@observed)),
                                         paste0("lag", 1:max_lag)))
  # fill:
  for(t in subset){
    phi_temp <- param$lambda[,,t]
    obs_preceding <- sts@observed[t - max_lag:1, , drop = FALSE]
    for(lag in 1:max_lag){
      inds <- seq(to = n_units*(max_lag - lag + 1), length.out = n_units)
      phi_this_lag <- phi_temp[, inds]
      contributions[t, , , lag] <- t(phi_this_lag)*matrix(obs_preceding[max_lag - lag + 1, ], ncol = n_units, nrow = n_units)
    }
  }
  return(contributions)
}

decompose_coarse <- function(fit){
  sts <- fit$stsObj
  max_lag <- if(class(fit)[1] == "hhh4lag") fit$max_lag else 1
  subset <- fit$control$subset
  n_units <- ncol(sts@observed)
  param <- hhh4addon:::lambda_tilde_complex_neighbourhood(fit, periodic = FALSE,
                                                          subset = 1:max(subset))
  decomposition_epidemic <- decompose_epidemic_component(fit)
  contributions_coarse <- array(dim = c(max(subset),
                                        5,
                                        n_units),
                                dimnames = list(1:max(subset),
                                                c("endemic",
                                                  "epidemic.self.lag1", "epidemic.self.higher_lags",
                                                  "epidemic.other.lag1", "epidemic.other.higher_lags"),
                                                paste0("to.", colnames(sts@observed))))
  for(t in subset){
    contributions_coarse[t, "endemic", ] <- param$nu[t, ]
    contributions_coarse[t, "epidemic.self.lag1", ] <- diag(decomposition_epidemic[t, , , 1])
    contributions_coarse[t, "epidemic.other.lag1", ] <-
      colSums(decomposition_epidemic[t, , , 1]) - contributions_coarse[t, "epidemic.self.lag1", ]
    if(max_lag > 1){
      contributions_higher_lags_temp <- apply(decomposition_epidemic[t, , , 2:max_lag], 1:2, sum)
      contributions_coarse[t, "epidemic.self.higher_lags", ] <- diag(contributions_higher_lags_temp)
      contributions_coarse[t, "epidemic.other.higher_lags", ] <-
        colSums(contributions_higher_lags_temp) - contributions_coarse[t, "epidemic.self.higher_lags", ]
    }else{
      contributions_coarse[t, "epidemic.self.higher_lags", ] <- contributions_coarse[t, "epidemic.other.higher_lags", ] <- 0
    }
  }
  return(contributions_coarse)
}

# fancy plot of fitted values:
plot_fit_strat <- function(fit, unit = 1, col = c("lightgrey", brewer.pal(n = 4, name = 'RdBu')[c(2, 1, 3, 4)]),
                           ylim = NULL, cex.points = 0.7, draw_yax = TRUE,...){
  polygon_from_curve <- function(tp, coords, ...){
    tp <- tp[!is.na(coords)]; tp <- c(tp[1], tp, tail(tp, 1))
    coords <- coords[!is.na(coords)]; coords <- c(0, coords, 0)
    polygon(tp, coords, ...)
  }
  obs <- fit$stsObj@observed[, unit]
  timepoints_calendar <- seq(from = fit$stsObj@start[1] + fit$stsObj@start[2]/fit$stsObj@freq,
                             length.out = nrow(fit$stsObj@observed), by  =1/fit$stsObj@freq)
  fitted_vals <- fitted.values(fit)[, unit]
  decomp <- decompose_coarse(fit)
  cumsums <- decomp
  if(is.null(ylim)) ylim <- c(0, max(c(cumsums[, , unit], obs), na.rm = TRUE))
  for(i in 2:5) cumsums[, i, ] <- cumsums[, i, ] + cumsums[, i - 1, ]
  plot(NULL, xlim = range(timepoints_calendar), ylim = ylim,
       xlab = NA, ylab = "no. of reported cases", axes = FALSE, ...)
  axis(1, labels = NA)
  axis(1, lty = 0, at = 2011:2017 + 0.5, labels = 2011:2017)
  if(draw_yax) axis(2)
  box()
  polygon_from_curve(timepoints_calendar, cumsums[, "epidemic.other.higher_lags", unit], col = col[5], border = col[5])
  polygon_from_curve(timepoints_calendar, cumsums[, "epidemic.other.lag1", unit], col = col[4], border = col[4])
  polygon_from_curve(timepoints_calendar, cumsums[, "epidemic.self.higher_lags", unit], col = col[3], border = col[3])
  polygon_from_curve(timepoints_calendar, cumsums[, "epidemic.self.lag1", unit], col = col[2], border = col[2])
  polygon_from_curve(timepoints_calendar, cumsums[, "endemic", unit], col = col[1], border = col[1])
  lines(timepoints_calendar[fit$control$subset], fitted_vals)
  points(timepoints_calendar, obs, pch = 16, cex = cex.points)
}

# adding the (within-sample) predictive interval:
add_pi <- function(fit, unit, upper, lower,...){
  timepoints_calendar <- seq(from = fit$stsObj@start[1] + fit$stsObj@start[2]/fit$stsObj@freq,
                             length.out = nrow(fit$stsObj@observed), by  =1/fit$stsObj@freq)
  subset <- fit$control$subset
  fitted_vals <- fitted.values(fit)[, unit]
  size <- exp(fit$coefficients["-log(overdisp)"])
  upper <- qnbinom(upper, mu = fitted_vals, size = size)
  lower <- qnbinom(lower, mu = fitted_vals, size = size)

  lines(timepoints_calendar[subset], upper, lty = "dotted")
  lines(timepoints_calendar[subset], lower, lty = "dotted")
}

# Functions for residual analyses:
# function to compute within-regions ACFs:
my_acf <- function(pearson_resids){
  acf_matr <- matrix(ncol = ncol(pearson_resids), nrow = 52,
                     dimnames = list(NULL, colnames(pearson_resids)))
  for(i in 1:ncol(pearson_resids)){
    acf_matr[, i] <- acf(pearson_resids[, i], lag.max = 52, plot = FALSE)$acf[-1,,1]
  }
  return(acf_matr)
}

# a customized plotting function:
myplot_acf <- function(macf, unit = 1, nobs = 358, shift_x = 0, lwd = 2, ylim = c(-0.2, 0.2), add = FALSE, ...){
  if(!add){
    plot(c(1:5, 7) + shift_x,
         macf[c(1:5, 52), unit],
         xlim = c(0.5, 7.5), ylim = ylim, type = "h",
         axes = FALSE, xlab = "", ylab = "residual ACF",
         lwd = lwd, ...)
    axis(1, at = c(1:5, 7), labels = c(1:5, 52))
    mtext("lag", 1, line = 2.2)
    axis(2)
    box()
    abline(h = 0)
    x <- 6; d <- 0.03; y_ax <- 1.1*min(ylim)
    rect(x - d, y_ax -d, x + d, y_ax + d, col = "white", border = NA, xpd = TRUE)
    lines(c(x - 2*d, x), y_ax + c(d/2, -d/2), xpd = TRUE)
    lines(c(x, x + 2*d), y_ax + c(d/2, -d/2), xpd = TRUE)
    abline(h = c(-1.96, 1.96)/sqrt(nobs), lty = "dotted")
  }else{
    lines(c(1:5, 7) + shift_x,
          macf[c(1:5, 52), unit], type = "h", lwd = lwd, ...)
  }
}

# Plotting functions for periodically stationary moments:
# Means:
plot_stat_means <- function(stat_mom, emp_mom, disease, unit = 1,...){
  plot(stat_mom[[disease]]$full$geom$mu_matrix[, unit], type = "l",
       xlab = "calendar week", ylab = "mean", col = cols_model_versions[1],
       axes = FALSE,...)
  axis(1, at = c(1, 13, 26, 39, 52)); axis(2); box()
  lines(stat_mom[[disease]]$simple_seas$geom$mu_matrix[, unit], col = cols_model_versions[2])
  lines(stat_mom[[disease]]$no_cross$geom$mu_matrix[, unit], col = cols_model_versions[3])
  lines(stat_mom[[disease]]$neither$geom$mu_matrix[, unit], col = cols_model_versions[4])
  lines(stat_mom[[disease]]$full$ar1$mu_matrix[, unit], col = cols_model_versions[5])
  # lines(stat_mom_unstrat[[disease]]$geom$mu_matrix, col = cols_model_versions[6], lty = "dashed")
  points(emp_mom[[disease]]$mean[, unit], pch = 1, cex = 0.6)
}

# Standard deviations:
plot_stat_sds <- function(stat_mom, emp_mom, disease, unit = 1,...){
  plot(sqrt(stat_mom[[disease]]$full$geom$var_matrix[, unit]), type = "l",
       xlab = "calendar week", ylab = "standard deviation", col = cols_model_versions[1],
       axes = FALSE, ...)
  axis(1, at = c(1, 13, 26, 39, 52)); axis(2); box()
  lines(sqrt(stat_mom[[disease]]$simple_seas$geom$var_matrix[, unit]), col = cols_model_versions[2])
  lines(sqrt(stat_mom[[disease]]$no_cross$geom$var_matrix[, unit]), col = cols_model_versions[3])
  lines(sqrt(stat_mom[[disease]]$neither$geom$var_matrix[, unit]), col = cols_model_versions[4])
  # lines(sqrt(stat_mom[[disease]]$full$ar1$var_matrix[, unit]), col = cols_model_versions[5])
  # lines(sqrt(stat_mom_unstrat[[disease]]$geom$var_matrix), col = cols_model_versions[6], lty = "dashed")
  points(sqrt(emp_mom[[disease]]$var[, unit]), pch = 1, cex = 0.6)
}

# means, with "confidence" bands:
plot_sm_bands <- function(sm_of_means, emp_mom,
                          disease, model_version, lag_structure,
                          unit, add = FALSE, col = "black", ylim = c(0, 20), ...){
  means <- sm_of_means[[disease]][[model_version]][[lag_structure]]$mu_matrix[, unit]
  sds <- sqrt(sm_of_means[[disease]][[model_version]][[lag_structure]]$var_matrix[, unit])
  upper <- means + 1.96*sds
  lower <- means - 1.96*sds

  if(!add){
    par(mar = c(1, 4.2, 4, 1))
    plot(NULL, ylim = ylim, xlim = c(1, 52),
         cex = 0.7, xlab = "", ylab = expression(hat(mu)[it]),
         axes = FALSE, ...)
    axis(2)
    box()
  }
  lines(sm_of_means[[disease]][[model_version]][[lag_structure]]$mu_matrix[, unit],
        col = col)
  lines(upper, lty = "dotted", col = col)
  lines(lower, lty = "dotted", col = col)
  points(emp_mom[[disease]]$mean[, unit], cex = 0.6)
}

# residuals relative to periodically stationary moments:
plot_stat_resids <- function(sm_of_means, emp_mom,
                             disease,
                             model_version1, lag_structure1,
                             model_version2, lag_structure2,
                             unit, add = FALSE, col = "black", ylim = c(-5, 5), ...){

  my_lines <- function(v1, v2, col1 = cols_model_versions[1], col2 = cols_model_versions[2]){
    inds1larger <- which(abs(v1) > abs(v2))
    inds2larger <- which(abs(v2) > abs(v1))
    lines(inds1larger, v1[inds1larger], type = "h", col = col1)
    lines(inds2larger, v2[inds2larger], type = "h", col = col2, lwd = 1)
    lines(inds2larger, v1[inds2larger], type = "h", col = col1)
    lines(inds1larger, v2[inds1larger], type = "h", col = col2, lwd = 1)
  }

  emp_means <- emp_mom[[disease]]$mean[, unit]

  means1 <- sm_of_means[[disease]][[model_version1]][[lag_structure1]]$mu_matrix[, unit]
  sds1 <- sqrt(sm_of_means[[disease]][[model_version1]][[lag_structure1]]$var_matrix[, unit])
  pr1 <- (emp_means - means1)/sds1

  means2 <- sm_of_means[[disease]][[model_version2]][[lag_structure2]]$mu_matrix[, unit]
  sds2 <- sqrt(sm_of_means[[disease]][[model_version2]][[lag_structure2]]$var_matrix[, unit])
  pr2 <- (emp_means - means2)/sds2

  plot(NULL, xlim = c(1, 52), axes = FALSE,
       type = "h", ylim = ylim, col = cols_model_versions[2],
       xlab = "calendar week", ylab = "Pearson residual",...)
  my_lines(pr1, pr2)
  axis(1, at = c(1, 13, 26, 39, 52)); axis(2); box()
  abline(h = 0, col = "black")
  abline(h = c(-1.96, 1.96), col = "black", lty = "dotted")
}


# function to simulate p-values
sim_teststats <- function(fit, n_sim, n_seasons, seed = 123){
  set.seed(seed)
  teststats_sim <- teststats_trafo_sim <- numeric(n_sim)
  for(i in 1:n_sim){
    # generate data and fit model to simulated data::
    if(is.null(fit$control$funct_lag)){
      dat_sim0 <- simulate(fit, y.start = fit$stsObj@observed[1, ]) # run once to go through "burn-in" period (independence of starting values)
      dat_sim <- simulate(fit, y.start = dat_sim0@observed[52 + 1, ]) # feed values from the first run into a second as starting values
      fit_sim <- hhh4(dat_sim, fit$control)
    }else{
      dat_sim0 <- hhh4addon:::simulate.hhh4lag(fit, y.start = fit$stsObj@observed[1:5, ])
      dat_sim <- hhh4addon:::simulate.hhh4lag(fit, y.start = dat_sim0@observed[52 + 1:5, ])
      fit_sim <- hhh4addon::profile_par_lag(dat_sim, fit$control)
    }

    # compute stationary moments of week-wise averages:
    sm_of_means_sim <- compute_sm_of_means(fit_sim, n_seasons = n_seasons, return_Sigma = TRUE)
    # compute empirical moments:
    emp_mom_sim <- get_emp_moments(dat_sim)
    # get test statistic:
    teststats_sim[i] <- get_pval(stat_mom_of_means = sm_of_means_sim, emp_mom = emp_mom_sim)$test_stat
    teststats_trafo_sim[i] <- get_pval2(stat_mom_of_means = sm_of_means_sim, emp_mom = emp_mom_sim)$test_stat
    # print(i)
  }
  return(list(teststats_sim = teststats_sim, teststats_trafo_sim = teststats_trafo_sim))
}
