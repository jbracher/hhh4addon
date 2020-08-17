#' Data set on norovirus gastroenteritis in Bremen and Lower Saxony
#'
#' Case counts of norovirus gastroenteritis in the German states of Bremen and Lower Saxony, 2011-2017; stored as an sts object
#'
#' @name noroBL
#' @docType data
#' @author Johannes Bracher
#' @source Surveillance counts retrieved from SurvStat@RKI 2.0 service (https://survstat.rki.de), Robert Koch Institute, Berlin as of 30 May 2017.
#' @keywords data
NULL

#' Data set on rotavirus gastroenteritis in Bremen and Lower Saxony
#'
#' Case counts of rotavirus gastroenteritis in the German states of Bremen and Lower Saxony, 2011-2017; stored as an sts object
#'
#' @name noroBL
#' @docType data
#' @author Johannes Bracher
#' @source Surveillance counts retrieved from SurvStat@RKI 2.0 service (https://survstat.rki.de), Robert Koch Institute, Berlin as of 30 May 2017.
#' @keywords data
NULL

#' Data set on campylobacteriosis in Bremen and Lower Saxony
#'
#' Case counts of campylobacteriosis in the German states of Bremen and Lower Saxony, 2011-2017; stored as an sts object
#'
#' @name noroBL
#' @docType data
#' @author Johannes Bracher
#' @source Surveillance counts retrieved from SurvStat@RKI 2.0 service (https://survstat.rki.de), Robert Koch Institute, Berlin as of 30 May 2017.
#' @keywords data
NULL


#' Data set on norovirus in Berlin
#'
#' Case counts of norovirus gastroenteritis in the twelve districts of Berlin, Germany, 2011-2017; stored as an sts object
#'
#' @name noroBL
#' @docType data
#' @author Johannes Bracher
#' @source Surveillance counts retrieved from SurvStat@RKI 2.0 service (https://survstat.rki.de), Robert Koch Institute, Berlin.
#' @keywords data
NULL

#' Data set on rotavirus in Berlin
#'
#' Case counts of rotavirus gastroenteritis in the twelve districts of Berlin, Germany, 2011-2017; stored as an sts object
#'
#' @name noroBL
#' @docType data
#' @author Johannes Bracher
#' @source Surveillance counts retrieved from SurvStat@RKI 2.0 service (https://survstat.rki.de), Robert Koch Institute, Berlin.
#' @keywords data
NULL

#' Data set on dengue in San Juan, Puerto Rico
#'
#' Case counts of dengue in San Juan, Puerto Rico, 1990-2013; stored as an sts object
#'
#' @name dengueSJ
#' @docType data
#' @author Johannes Bracher
#' @source Counts retrieved from the supplement of Ray et al (2017): Infectious disease prediction with kernel conditional density estimation, Statistics in Medicine 36(30):4908-4929.
#' These data originally stem from a forecasting competition organized by the US federal government: http://dengueforecasting.noaa.gov/
#' @keywords data
NULL

#' Check whether the rows of a matrix show a cyclic pattern
#'
#' Needed to determine whether \code{stationary_moments} is applicable (works only for
#' models with periodic parameter structure)
#'
#' @param matr The parameter matrix to check.
#' @param length_of_period Usually 52 (52 weeks per year).
#' @return logical: does the matrix show a periodic pattern?
matrix_is_cyclic <- function(matr, length_of_period){
  n_timepoints <- nrow(matr)
  cyclic_for_ith_timepoint <- rep(NA, length_of_period)
  for(i in 1:length_of_period){
    rows_i <- seq(from = i, to = n_timepoints, by = length_of_period)
    cyclic_for_ith_timepoint[i] <- all(apply(matr[rows_i, , drop = FALSE], 2, function(x) all.equal(min(x), max(x)) == TRUE))
  }
  return(all(cyclic_for_ith_timepoint))
}

#' Get diagnoal elements of all slices of an array
#'
#' Extracts diagonals of all slices of an array (i.e. of \code{arr[,,1], arr[,,2], ...} and stacks them in one vector.)
#'
#' @param arr An array.
#'
get_diags_of_array <- function(arr){
  n_units <- dim(arr)[1]
  lgt <- dim(arr)[3]
  inds_one_slice <- seq(from = 1, to = n_units^2, length.out = n_units)
  inds_in_array <- rep(0:(lgt - 1)*n_units^2, each = n_units) + inds_one_slice
  return(arr[inds_in_array])
}

#' Check if the par_lag parameter was fitted
is_fitted_par_lag <- function(object){
  return(sum(object$dim) > length(object$coefficients))
}
