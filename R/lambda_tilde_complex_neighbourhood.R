# get the indices corresponding to the first full period
get_inds_period <- function(hhh4Obj){
  length_of_period <- hhh4Obj$stsObj@freq
  inYear <- surveillance::epochInYear(hhh4Obj$stsObj)
  ind_first_week1 <- min(which(inYear == 1))
  ind_first_full_year <- ind_first_week1 + seq(from = 0, length.out = length_of_period)
  return(ind_first_full_year)
}

# function to check whether lags are handled in a uniform fashion
check_lags <- function(hhh4Obj){
  ar <- hhh4Obj$control$ar
  ne <- hhh4Obj$control$ne

  # check compatibility with lag specification:
  # for hhh4 class lags need to be 1 (could be generalized, but nobody ever uses the others)
  if(class(hhh4Obj)[1] == "hhh4"){
    if((ne$lag != ar$lag) & (ne$f != ~-1) & (ar$f != ~-1)){
      stop("Lags in all components need to be the same for these algorithms to work.")
    }
  }
}

### Change this in a way that it returns a list of arrays nu and lambda.

#' Extracting Lambda_Tilde from an hhh4 object with complex neighbourhood structure
#'
#' Extracting Lambda_Tilde from an hhh4 object with complex neighbourhood structure. Used for calculations of
#' longterm predictions and stationary distributions.
#'
#' @param hhh4Obj a hhh4 object for which to extract Lambda_tilde
#' @param subset a subset (in time); only required when periodic == FALSE
#' @param periodic choose subset to correspond to one full cycle
lambda_tilde_complex_neighbourhood <- function(hhh4Obj, subset = NULL, periodic = FALSE){

  check_lags(hhh4Obj)

  # extract information from model:
  n_units <- ncol(hhh4Obj$stsObj@observed)
  n_timepoints <- nrow(hhh4Obj$stsObj@observed)

  # means etc:
  meanHHH_temp <- surveillance:::meanHHH(hhh4Obj$coefficients,
                                         terms(hhh4Obj),
                                         subset = 1:n_timepoints)
  # lagged version:
  Ylagged <- terms(hhh4Obj)$offset$ar/hhh4Obj$control$ar$offset # the offset in terms also contains the offset from control
  # maximium lag:
  max_lag <- max(c(if(class(hhh4Obj)[1] == "hhh4lag")hhh4Obj$control$max_lag,
                   hhh4Obj$control$ar$lag,
                   hhh4Obj$control$ne$lag), na.rm = TRUE) #BJ

  # weights of lags:
  weights_lag <- if(!is.null(hhh4Obj$distr_lag)){
    hhh4Obj$distr_lag # have to be identical in ar and ne, checked by check_lags
  }else{c(rep(0, max_lag - 1), 1)} # all weight to one lag if no distributed lags are used.

  # (potentially time-varying) nu:
  nu <- meanHHH_temp$endemic
  # (potentially time-varying) lambda:
  lambda <- meanHHH_temp$ar.exppred*hhh4Obj$control$ar$offset
  # (potentially time-varying) Phi:
  Phi <- meanHHH_temp$ne.exppred*hhh4Obj$control$ne$offset
  # (potentially normalized, potentially scaled) weights
  scaled_w <- surveillance:::getNEweights(hhh4Obj)
  if(is.null(scaled_w)){ # if nothing provided: can be set to zero matrix
    scaled_w <- diag(0, n_units)
  }

  # if periodic: choose first full period as subset
  if(periodic){
    # extract additional info:
    subset <- get_inds_period(hhh4Obj = hhh4Obj)
    # check whether model is actually periodic:
    if(!matrix_is_cyclic(matr = nu, length_of_period = length(subset)) |
       !matrix_is_cyclic(matr = Phi, length_of_period = length(subset))){
      stop("hhh4Obj does not seem to have a periodic structure, funtions not applicable.")
    }
  }else{
    if(is.null(subset)){
      stop("When periodic == FALSE a value for subset must be provided.")
    }
  }

  # length of subset:
  n_subset <- length(subset)

  # if subset is outside of time range of hhh4 model:
  if(max(subset)>n_timepoints){
    # add NAs to observed
    hhh4Obj$stsObj@observed <- rbind(hhh4Obj$stsObj@observed, matrix(ncol = 3, nrow = max(subset) - n_timepoints))
    # and extend time
    hhh4Obj$control$data$t <- c(hhh4Obj$control$data$t,
                                seq(from = max(hhh4Obj$control$data$t) + 1, to = max(subset) - 1, by = 1))
    n_timepoints <- nrow(hhh4Obj$stsObj@observed)
  }

  # initialise Lambda:
  Lambda_wide <- array(NA, dim = c(n_units, n_units*max_lag, n_subset))
  # to check whether Lambda_tilde is in line with the fitted means of the model: calculate implied epidemic components
  means_from_lambda_tilde <- matrix(NA, nrow = n_subset, ncol = n_units)

  # Fill Lambda:
  for(i in 1:n_subset){
    t <- subset[i]
    # accout for potentially time-varying scaled_w
    if(length(dim(scaled_w)) == 3){
      scaled_wt <- scaled_w[,,t]
    }else{
      scaled_wt <- scaled_w
    }
    # add phi
    lambda_unweighted_temp <- Phi[t, ]*t(scaled_wt)
    # add lambda (in case that three-component-version is used)
    diag(lambda_unweighted_temp) <- diag(lambda_unweighted_temp) + lambda[t, ]

    for(j in 1:max_lag){
      inds <- seq(to = (max_lag - j + 1)*n_units, length.out = n_units)
      Lambda_wide[,inds,i] <- weights_lag[j]*lambda_unweighted_temp
    }

    # to check correctness:
    if(i > max_lag){
      means_from_lambda_tilde[i, ] <- nu[t, ] + lambda_unweighted_temp %*% Ylagged[t, ]
    }
  }

  # restrict nu to subset:
  nu <- nu[subset, , drop = FALSE]

  # name elements:
  labels_t <- paste0("-", (max_lag - 1):0); labels_t[length(labels_t)] <- ""
  colnames(Lambda_wide) <- paste0(rep(colnames(hhh4Obj$stsObj@observed)), rep(labels_t, each = n_units))
  rownames(Lambda_wide) <- colnames(nu) <- colnames(hhh4Obj$stsObj)

  if(periodic){
    dimnames(Lambda_wide)[[3]] <- rownames(nu) <- paste0("m=", 1:n_subset)
  }else{
    dimnames(Lambda_wide)[[3]] <- rownames(nu) <- paste0("t=", subset)
  }

  # check whether everything went right:
  if(any(abs(as.vector(meanHHH_temp$mean[subset[-(1:max_lag)], , drop = FALSE]) -
             as.vector(means_from_lambda_tilde[-(1:max_lag), , drop = FALSE])) > 0.00001, na.rm = TRUE)){
    # added na.rm = TRUE so that this works for models where the last observations are not yet available.
    stop("Extracted Lambda is not in agreement with fitted values returned by surveillance:::meanHHH.
         Model does not seem to be covered by extraction algorithm")
  }# else{print("Equality check in lambda_tilde_complex_neighbourhood passed")}

  return(list(nu = nu, lambda = Lambda_wide))
}
