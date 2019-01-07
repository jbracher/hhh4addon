##### This file serves to reproduce parts of the analysis from
##### Bracher/Held (2018+): Multivariate endemic-epidemic models
##### with higher-order lags and an application to outbreak detection

# install the package hhh4addon:
# library(devtools)
# install_github("jbracher/hhh4addon", build_vignettes = TRUE)

#########################################
# get and format data:
library(surveillance)
library(hhh4addon)
library(RColorBrewer)

setwd("/home/johannes/Documents/hhh4predict/Theory/Article_Theory/data_analysis")
# The file auxiliary functions.R needs to be available (equally available in example_analyses.R)
source("auxiliary_functions.R")

library(hhh4contacts)
data("noroBE")
data("rotaBE")

# modify neighbourhood matrices to apply power law:
noroBE@neighbourhood <- noroBE@neighbourhood + 1
rotaBE@neighbourhood <- rotaBE@neighbourhood + 1

data <- list(noro = noroBE, rota = rotaBE)
names_diseases <- names(data)

library(RColorBrewer)
cols_model_versions <- brewer.pal(8, "Dark2")[1:4]
cols_lag_structures <- brewer.pal(8, "Dark2")[c(8, 5:7)]

#########################################
# fit models:

# specify controls:
names_model_versions <- c("full", "simple_seas", "no_cross", "neither")
names_lag_structures <- c("ar1", "ar2", "geom", "pois")
cols_model_versions <- brewer.pal(8, "Dark2")[1:4]
cols_lag_structures <- brewer.pal(8, "Dark2")[c(8, 5:7)]

# define variable for Christmas breaks:
christmas <- numeric(nrow(noroBE@observed))
christmas[(seq_along(christmas) %% 52) %in% c(0, 1)] <- 1

# number of sine/cosine waves to include in the endemic and epidemic components:
S_end <- 1
S_epid <- 1

# full model (a)
ctrl_full <- list(end = list(f = addSeason2formula(~0 + fe(1, unitSpecific = TRUE) + christmas, S = S_end)),
                  ar = list(f = ~-1),
                  ne = list(f = addSeason2formula(~0 + fe(1, unitSpecific = TRUE) + christmas, S = S_epid),
                            weights = W_powerlaw(maxlag=5, normalize=TRUE, log=TRUE)),
                  family = "NegBin1",
                  subset = 6:nrow(noroBE@observed),
                  data = list(christmas = christmas))

# model (b) with simplified seasonality (no seasonality in epidemic parameters)
ctrl_simple_seas <- list(end = list(f = addSeason2formula(~0 + fe(1, unitSpecific = TRUE) + christmas, S = S_end)),
                         ar = list(f = ~-1),
                         ne = list(f = ~0 + fe(1, unitSpecific = TRUE),
                                   weights = W_powerlaw(maxlag=5, normalize=TRUE, log=TRUE)),
                         family = "NegBin1",
                         subset = 6:nrow(noroBE@observed),
                         data = list(christmas = christmas))

# model (c) without cross lags
ctrl_no_cross <- list(end = list(f = addSeason2formula(~0 + fe(1, unitSpecific = TRUE) + christmas, S = S_end)),
                      ne = list(f = ~-1),
                      ar = list(f = addSeason2formula(~0 + fe(1, unitSpecific = TRUE) + christmas, S = S_epid)),
                      family = "NegBin1",
                      subset = 6:nrow(noroBE@observed),
                      data = list(christmas = christmas))

# model (d) with neither cross lags nor seasonality in epidemic parameters
ctrl_neither <- list(end = list(f = addSeason2formula(~0 + fe(1, unitSpecific = TRUE) + christmas, S = S_end)),
                     ne = list(f = ~-1),
                     ar = list(f = ~0 + fe(1, unitSpecific = TRUE)),
                     family = "NegBin1",
                     subset = 6:nrow(noroBE@observed),
                     data = list(christmas = christmas))

controls <- list(full = ctrl_full,
                 simple_seas = ctrl_simple_seas,
                 no_cross = ctrl_no_cross,
                 neither = ctrl_neither)

# fit all models (this takes a while):
return_full_cov <- FALSE
fits <- full_cov_matrices <- list()
for(disease in names_diseases){
  print(paste("--- Started", disease, "---"))
  for(model_version in names_model_versions){
    # ar1 lag structure:
    print(paste(model_version, "ar1"))
    fits[[disease]][[model_version]]$ar1 <- hhh4(data[[disease]], controls[[model_version]])
    # ar2 lag structure:
    print(paste(model_version, "ar2"))
    control_ar2_temp <- controls[[model_version]]; control_ar2_temp$funct_lag <- ar2_lag
    fit_temp <- profile_par_lag(data[[disease]], control_ar2_temp, return_full_cov = return_full_cov)
    if(return_full_cov){
      fits[[disease]][[model_version]]$ar2 <- fit_temp$best_mod
      full_cov_matrices[[disease]][[model_version]]$ar2 <- fit_temp$cov
    }else{
      fits[[disease]][[model_version]]$ar2 <- fit_temp
    }

    # geometric lag structure:
    print(paste(model_version, "geom"))
    control_geom_temp <- controls[[model_version]]; control_geom_temp$funct_lag <- geometric_lag
    fit_temp <- profile_par_lag(data[[disease]], control_geom_temp, return_full_cov = return_full_cov)
    if(return_full_cov){
      fits[[disease]][[model_version]]$geom <- fit_temp$best_mod
      full_cov_matrices[[disease]][[model_version]]$geom <- fit_temp$cov
    }else{
      fits[[disease]][[model_version]]$geom <- fit_temp
    }


    # Poisson lag structure:
    print(paste(model_version, "pois"))
    control_pois_temp <- controls[[model_version]]; control_pois_temp$funct_lag <- poisson_lag
    fit_temp <- profile_par_lag(data[[disease]], control_pois_temp, return_full_cov = return_full_cov)
    if(return_full_cov){
      fits[[disease]][[model_version]]$pois <- fit_temp$best_mod
      full_cov_matrices[[disease]][[model_version]]$pois <- fit_temp$cov
    }else{
      fits[[disease]][[model_version]]$pois <- fit_temp
    }
  }
}

#########################################
# Assessing models fits etc.

# obtain AICs, lag weights, neighbourhood weights:
AICs <- overdisp <- neighbourhood_weights <- lag_weights <- list()

for(disease in names_diseases){
  AICs[[disease]] <- overdisp[[disease]] <-
    matrix(nrow = length(names_lag_structures), ncol = length(names_model_versions),
           dimnames = list(names_lag_structures, names_model_versions))
  neighbourhood_weights[[disease]] <- lag_weights[[disease]] <-
    array(dim = c(length(names_lag_structures), length(names_model_versions), 5),
          dimnames = list(names_lag_structures, names_model_versions, paste0("lag.", 1:5)))
  for(model_version in names_model_versions){
    for(lag_structure in names_lag_structures){
      AICs[[disease]][lag_structure, model_version] <- AIC(fits[[disease]][[model_version]][[lag_structure]])

      overdisp[[disease]][lag_structure, model_version] <-
        fits[[disease]][[model_version]][[lag_structure]]$coefficients["-log(overdisp)"]

      neighbourhood_weights[[disease]][lag_structure, model_version, ] <-
        (1:5)^-exp(fits[[disease]][[model_version]][[lag_structure]]$coefficients["neweights.d"])

      lag_weights[[disease]][lag_structure, model_version, ] <-  if(lag_structure == "ar1"){
        c(1, 0, 0, 0, 0)
      }else{
        fits[[disease]][[model_version]][[lag_structure]]$distr_lag
      }

    }
  }
}

lag_weights

Delta_AICs_noro <- AICs$noro - max(AICs$noro)
Delta_AICs_rota <- AICs$rota - max(AICs$rota)

# Observations:
# - more weight on higher lags
# - different weighting schemes now differ a bit more
# - weighted lag models borrow less strength across regions
# - slightly decreased overdispersion, i.e. less unexplained variation

#########################################
# Analyse periodically stationary moments:

# compute empirical moments:
emp_mom <- list()
emp_mom$noro <- get_emp_moments(noroBE, start = 1)
emp_mom$rota <- get_emp_moments(rotaBE, start = 1)

# compute model-based moments (both for single observations and calendar-week-wise means):
sm <- sm_of_means <- list()
for(disease in names_diseases){
  print(disease)
  for(model_version in names_model_versions){
    print(model_version)
    for(lag_structure in names_lag_structures){
      print(lag_structure)
      sm[[disease]][[model_version]][[lag_structure]] <-
        stationary_moments(fits[[disease]][[model_version]][[lag_structure]], return_Sigma = TRUE)
      sm_of_means[[disease]][[model_version]][[lag_structure]] <-
        compute_sm_of_means(fits[[disease]][[model_version]][[lag_structure]], n_seasons = 7, return_Sigma = TRUE)
    }
  }
}
# save(sm_of_means, file = "../model_fits/sm_of_means_nb1.rda") # with faster version actually no need to store this separately.
# load("sm_of_means_nb1.rda")

# plots for comparison:
disease <- "rota"
# means:
par(mfrow = 3:4)
for(unit in 1:12){
  plot_stat_means(sm, emp_mom, "rota", unit)
}

par(mfrow = 3:4)
par(mfrow = 3:4)
for(unit in 1:12){
  plot_stat_sds(sm, emp_mom, "rota", unit)
}

#########################################
# p-values

# Compute approximative p-values as given in the Supplementary Material
# (the simulation-based ones from the article were computed
# using some custom functions not included in the package and computations are quite long)
teststats <- pvals_approx_liberal <- pvals_approx_conservative <-
  teststats_trafo <- pvals_trafo_liberal <- pvals_trafo_conservative <- list()
for(disease in names_diseases){
  pvals_approx_liberal[[disease]] <- teststats[[disease]] <-
    matrix(nrow = length(names_lag_structures), ncol = length(names_model_versions),
           dimnames = list(names_lag_structures, names_model_versions))
  for(model_version in names_model_versions){
    for(lag_structure in names_lag_structures){
      testres_liberal_temp <- get_pval(stat_mom_of_means = sm_of_means[[disease]][[model_version]][[lag_structure]],
                                       emp_mom = emp_mom[[disease]],
                                       correction_df = fits[[disease]][[model_version]][[lag_structure]]$dim[1])
      pvals_approx_liberal[[disease]][lag_structure, model_version] <- testres_liberal_temp$pval
      teststats[[disease]][lag_structure, model_version] <- testres_liberal_temp$test_stat
    }
  }
}


########################################
# Analysis of conventional Pearson residuals

pearson_resids <- acfs <- list()
for(disease in names_diseases){
  for(model_version in names_model_versions){
    for(lag_structure in names_lag_structures){
      size_temp <- exp(-fits[[disease]][[model_version]][[lag_structure]]$coefficients["-log(overdisp)"])
      fitted_values_temp <- fits[[disease]][[model_version]][[lag_structure]]$fitted.values
      cond_sds_temp <- sqrt(fitted_values_temp + size_temp*fitted_values_temp^2)
      obs_values_temp <- data[[disease]]@observed[-(1:5), ]
      pearson_resids_temp <- (fitted_values_temp - obs_values_temp)/cond_sds_temp
      pearson_resids[[disease]][[lag_structure]][[model_version]] <- pearson_resids_temp
      pearson_resids[[disease]][[lag_structure]][[model_version]] <- my_acf(pearson_resids_temp)
    }
  }
}

# choose which model to display:
model_version <- "full"
cols_lag_structures <- c("black", "orange", "purple", "chartreuse3")
lwd_weights <- 2
unit1 <- "pank"
unit2 <- "span"

par(mfrow = c(2, 2))
myplot_acf(pearson_resids$noro$ar1[[model_version]], unit = unit1, main = paste("(a) Norovirus,", unit1),
           shift_x = -0.15, col = cols_lag_structures[1])
myplot_acf(pearson_resids$noro$ar2[[model_version]], unit = unit1,
           add = TRUE, shift_x = -0.05, col = cols_lag_structures[2])
myplot_acf(pearson_resids$noro$geom[[model_version]], unit = unit1,
           add = TRUE, shift_x = 0.05, col = cols_lag_structures[3])
myplot_acf(pearson_resids$noro$pois[[model_version]], unit = unit1,
           add = TRUE, shift_x = 0.15, col = cols_lag_structures[4])

myplot_acf(pearson_resids$noro$ar1[[model_version]], unit = "span", main = paste("(b) Norovirus,", unit2),
           shift_x = -0.15, col = cols_lag_structures[1])
myplot_acf(pearson_resids$noro$ar2[[model_version]], unit = "span",
           add = TRUE, shift_x = -0.05, col = cols_lag_structures[2])
myplot_acf(pearson_resids$noro$geom[[model_version]], unit = "span",
           add = TRUE, shift_x = 0.05, col = cols_lag_structures[3])
myplot_acf(pearson_resids$noro$pois[[model_version]], unit = "span",
           add = TRUE, shift_x = 0.15, col = cols_lag_structures[4])

myplot_acf(pearson_resids$rota$ar1[[model_version]], unit = unit1, main = paste("(c) Rotavirus,", unit1),
           shift_x = -0.15, col = cols_lag_structures[1])
myplot_acf(pearson_resids$rota$ar2[[model_version]], unit = unit1,
           add = TRUE, shift_x = -0.05, col = cols_lag_structures[2])
myplot_acf(pearson_resids$rota$geom[[model_version]], unit = unit1,
           add = TRUE, shift_x = 0.05, col = cols_lag_structures[3])
myplot_acf(pearson_resids$rota$pois[[model_version]], unit = unit1,
           add = TRUE, shift_x = 0.15, col = cols_lag_structures[4])

myplot_acf(pearson_resids$rota$ar1[[model_version]], unit = unit2, main = paste("(d) Rotavirus,", unit2),
           shift_x = -0.15, col = cols_lag_structures[1])
myplot_acf(pearson_resids$rota$ar2[[model_version]], unit = unit2,
           add = TRUE, shift_x = -0.05, col = cols_lag_structures[2])
myplot_acf(pearson_resids$rota$geom[[model_version]], unit = unit2,
           add = TRUE, shift_x = 0.05, col = cols_lag_structures[3])
myplot_acf(pearson_resids$rota$pois[[model_version]], unit = unit2,
           add = TRUE, shift_x = 0.15, col = cols_lag_structures[4])


##########
# New Pearson residuals relative to stationary moments
library(RColorBrewer)

# noro

par(mar = c(1, 4.5, 2, 1), font.main = 1, family = "serif", las = 1)
layout_matr <- matrix(c(1, 2,
                        1, 2,
                        1, 2,
                        3, 4,
                        3, 4,
                        3, 4,
                        5, 6,
                        5, 6,
                        5, 6,
                        5, 6#,
                        # 7, 7,
                        # 7, 7
), byrow = TRUE, ncol = 2)
layout(layout_matr)
ylim_mu <- c(0, 25)
ylim_sd <- c(0, 15)

plot_stat_means(stat_mom = sm, emp_mom = emp_mom, "noro", unit = unit1, ylim = ylim_mu)
title("(a) Pankow")
legend("top", col = c(cols_model_versions[1:4], "black"),
       legend = c("Model (a)   ", "Model (b)",
                  "Model (c)   ", "Model (d)"),
       lty = c(1, 1, 1, 1), pch = c(NA, NA, NA, NA), ncol = 2, cex = 1, bty = "n")

plot_stat_means(stat_mom = sm, emp_mom = emp_mom, "noro", unit = unit2, ylim = ylim_mu)
title("(b) Spandau")

plot_stat_sds(stat_mom = sm, emp_mom = emp_mom, "noro", unit = unit1, ylim = ylim_sd)
# mtext("(a) Pankow", line = 1.3, at = -10, cex = 0.8)
plot_stat_sds(stat_mom = sm, emp_mom = emp_mom, "noro", unit = unit2, ylim = ylim_sd)
# mtext("(b) Spandau", line = 1.3, at = -10, cex = 0.8)

par(mar = c(4.3, 4.5, 2, 1), font.main = 1, family = "serif", las = 1)
plot_stat_resids(sm_of_means, emp_mom, "noro",
                 "full", "geom",
                 "simple_seas", "geom",
                 unit1)

plot_stat_resids(sm_of_means, emp_mom, "noro",
                 "full", "geom",
                 "simple_seas", "geom",
                 unit2)

# rota:

plot_stat_means(stat_mom = sm, emp_mom = emp_mom, "rota", unit = unit1, ylim = ylim_mu)
title("(a) Pankow")
legend("topright", col = c(cols_model_versions[1:4], "black"),
       legend = c("Model (a)   ", "Model (b)",
                  "Model (c)   ", "Model (d)"),
       lty = c(1, 1, 1, 1), pch = c(NA, NA, NA, NA), ncol = 2, cex = 1, bty = "n")

plot_stat_means(stat_mom = sm, emp_mom = emp_mom, "rota", unit = unit2, ylim = ylim_mu)
title("(b) Spandau")


plot_stat_sds(stat_mom = sm, emp_mom = emp_mom, "rota", unit = unit1, ylim = ylim_sd)
# mtext("(a) Pankow", line = 1.3, at = -10, cex = 0.8)
plot_stat_sds(stat_mom = sm, emp_mom = emp_mom, "rota", unit = unit2, ylim = ylim_sd)
# mtext("(b) Spandau", line = 1.3, at = -10, cex = 0.8)

par(mar = c(4.3, 4.5, 2, 1), font.main = 1, family = "serif", las = 1)
plot_stat_resids(sm_of_means, emp_mom, "rota",
                 "full", "geom",
                 "simple_seas", "geom",
                 unit1)

plot_stat_resids(sm_of_means, emp_mom, "rota",
                 "full", "geom",
                 "simple_seas", "geom",
                 unit2)

