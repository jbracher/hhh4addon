detect_outbreak <- function(stsObj, control, range, alpha = 0.01, b = 3, pastWeeksNotIncluded = 26, type = c("seasonal", "osa")){
  n_units <- ncol(stsObj@observed)
  freq <- stsObj@freq
  mu0_matrix <- var0_matrix <- size0_matrix <- thresholds <-
    pvalue <- matrix(nrow = length(range), ncol = n_units)

    for(t in seq_along(range)){
    # get subset representing training data:
    control$subset <- seq(from = range[t] - b*freq, to = range[t] - pastWeeksNotIncluded)
    # fit model:
    fit_temp <- surveillance::hhh4(stsObj = stsObj, control = control)
    # if convergence fails: try Poisson model
    if(fit_temp$convergence == FALSE){
      print(t)
      control_new <- control
      control_new$start <- NULL
      fit_temp <- hhh4(stsObj, control = control_new)
    }
    # get the longterm prediction (at horizon pastWeeksNotIncluded):
    pred_moments_temp <- hhh4addon::predictive_moments(fit_temp,
                                                  t_condition = tail(control$subset, 1),
                                                  lgt = pastWeeksNotIncluded)
    if(fit_temp$convergence){
      mu0_matrix[t, ] <- pred_moments_temp$mu_matrix[pastWeeksNotIncluded, ]
      var0_matrix[t, ] <- pred_moments_temp$var_matrix[pastWeeksNotIncluded, ]
      size0_matrix[t, ] <- pmin(abs(mu0_matrix[t, ]^2/(var0_matrix[t, ] - mu0_matrix[t, ])), 10000)
      thresholds[t, ] <- qnbinom(1 - alpha, mu = mu0_matrix[t, ], size = size0_matrix[t, ])
      pvalue[t, ] <- pnbinom(pred_moments_temp$realizations_matrix[pastWeeksNotIncluded, ], mu = mu0_matrix[t, ], size = size0_matrix[t, ])
      # update start for next iteration:
      control$start$fixed <- fit_temp$coefficients
    }
  }
  # create return object:
  # make subset of sts:
  stsObj <- stsObj[range, ]
  # add slot "control" to sts object
  stsObj@control <- list()
  stsObj@control$alpha <- alpha
  stsObj@control$mu0_matrix <- stsObj@control$expected <- mu0_matrix
  stsObj@control$var0_matrix <- var0_matrix
  stsObj@control$size0_matrix <- size0_matrix
  stsObj@control$thresholds <- thresholds
  stsObj@control$pvalue <- pvalue
  stsObj@alarm <- stsObj@observed > thresholds
  return(stsObj)
}

plot_outbreak_detection <- function(stsObj, unit){
  timepoints <- seq(from = stsObj@start[1] + stsObj@start[2]/stsObj@freq, length.out = nrow(stsObj@observed),
           by = 1/stsObj@freq)
  ylim <- range(c(stsObj@control$thresholds[, unit], stsObj@observed[, unit]), na.rm = TRUE); ylim[1] <- 0
  plot(timepoints, stsObj@control$thresholds[, unit], ylim = ylim, type = "l")
  cols <- rep("black", nrow(stsObj@observed))
  cols[stsObj@observed[, unit] > stsObj@control$thresholds[, unit]] <- "red"
  cols[is.na(stsObj@control$thresholds[, unit])] <- "grey"
  lines(timepoints, stsObj@control$mu0_matrix[, unit], lty = "dashed")
  points(timepoints, stsObj@observed[, unit], col = cols, pch = 16)
  legend("topleft", legend = c("expectation",
                               paste(round(100*(1 - stsObj@control$alpha), 2), "% quantile"),
                               "flagged", "not converged"),
         lty = c("solid", "dashed", NA, NA), pch = c(NA, NA, 16, 16),
         col = c("black", "black", "red", "grey"))
}
