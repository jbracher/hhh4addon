#' Determine whether an hhh4 object was fitted using one of the more complex techniques for handling neighbourhoods
#'
#' @param hhh4Obj an hhh4 object
is_complex_neighbourhood <- function(hhh4Obj){
  if(!is.null(hhh4Obj$logpower) |
     is.list(hhh4Obj$control$ne$weights) |
     any(grepl("neweights", names(hhh4Obj$coefficients)))){
    return(TRUE)
  }
  if(is.matrix(hhh4Obj$control$ne$weights)){
    if(any(diag(hhh4Obj$control$ne$weights) > 0)){
      return(TRUE)
    }
  }
  return(FALSE)
}

#' Extracting Lambda_Tilde from an hhh4 object with complex neighbourhood structure
#'
#' A wrapper around \code{lambda_tilde_complex_neighbourhood} and \code{lambda_tilde_simple_neighbourhood}.
#'
#' @param hhh4Obj a hhh4 object for which to extract Lambda_tilde
#' @param subset a subset (in time); only required when periodic == FALSE
#' @param periodic choose subset to correspond to one full cycle

lambda_tilde <- function(hhh4Obj, subset = NULL, periodic = FALSE){
  # check which case to handle:
  # complex_neighbourhood <- is_complex_neighbourhood(hhh4Obj = hhh4Obj)
  # more complex case of Sebastian's power law and social contact models
  # if(complex_neighbourhood){
    lambda_tilde_complex_neighbourhood(hhh4Obj = hhh4Obj, subset = subset, periodic = periodic)
  # }else{# simpler case with binary neighbourhood matrix and no social contact matrix:
  #   lambda_tilde_simple_neighbourhood(hhh4Obj = hhh4Obj, subset = subset, periodic = periodic)
  # }
}
