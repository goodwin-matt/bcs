#' Finds the discrete wavelet transform matrix.
#'
#' @param PHI Measurement matrix.
#' @param y Sample from original signal.
#' @param eta Tolerance level in determining convergence of marginal likelihood.
#' @param verbose Print to screen which basis are added, re-estimated, or deleted.
#' @return The sparse signal as found by the fast Laplace algorithm.
#' @export
WaveletBasis <- function(train,N,signal,wavelet){
  # Function for finding the DWT matrix
  # Args:
  #   train: the time indice points to evaluate the basis functions at
  #   N: the number of basis functions to find
  #   signal: the original signal be decomposed into basis functions.
  #   wavelet: the type of wavelet for the decomposition.
  # Returns:
  #   The DWT matrix evaluated at the training point indices
  w <- wavDWT(signal,wavelet=wavelet)
  M <- length(train)
  lvl <- w$dictionary$n.levels
  #basis <- matrix(NA,2^lvl,2^lvl)
  basis <- matrix(NA,M,N)
  compt <- 1
  w$data[lvl+1][[1]] <- 0
  for (i in 1:lvl){
    w$data[[i]]<-rep(0,2^(lvl-i))
  }
  for (i in 1:(lvl+1)){
    for (j in 1:2^(lvl-i)){
      w$data[[i]][j] <- 1
      basis[,compt] <- reconstruct(w)[train]
      w$data[[i]][j] <- 0
      compt <- compt + 1
    }
  }
  return(basis)
}


#' Finds the discrete Fourier transformation matrix.
#'
#' Given samples from some signal and given a measurement matrix, find the
#' sparse respresentation of a signal using compressive sensing.
#'
#' @param PHI Measurement matrix.
#' @param y Sample from original signal.
#' @param eta Tolerance level in determining convergence of marginal likelihood.
#' @param verbose Print to screen which basis are added, re-estimated, or deleted.
#' @return The sparse signal as found by the fast Laplace algorithm.
#' @export
FourierBasis <- function(train,tlist,M,N){
  # Function for finding the DFT matrix
  # Args:
  #   train: the time indice points to evaluate the basis functions at
  #   tlist: an array of the specific time points where the signal is evaluated at
  #   M: how many time points there are, the length of tlist
  #   N: how many basis to keep (how many columns in the matrix)
  # Returns:
  #   The DFT matrix evaluated at the training point indices
  if(N%%2==0)
    bFor <- create.fourier.basis(c(tlist[1],tlist[M]),N,dropind=N-1)
  else
    bFor <- create.fourier.basis(c(tlist[1],tlist[M]),N)
  return(eval.basis(tlist[train],bFor))
}

#' Finds the transformation matrix for b-spline basis..
#'
#' Given samples from some signal and given a measurement matrix, find the
#' sparse respresentation of a signal using compressive sensing.
#'
#' @param PHI Measurement matrix.
#' @param y Sample from original signal.
#' @param eta Tolerance level in determining convergence of marginal likelihood.
#' @param verbose Print to screen which basis are added, re-estimated, or deleted.
#' @return The sparse signal as found by the fast Laplace algorithm.
#' @export
BSplineBasis <- function(train,tlist,M,N){
  # Function for finding the DFT matrix
  # Args:
  #   train: the time indice points to evaluate the basis functions at
  #   tlist: an array of the specific time points where the signal is evaluated at
  #   M: how many time points there are (how many rows in the matrix)
  #   N: how many basis to keep (how many columns in the matrix)
  # Returns:
  #   The matrix representing the b-spline basis evaluated at the training point indices
  # The first two lines below shift the time scale so that the difference
  # between knots is roughly equal to 1. To see where the knots are placed
  # take bSpl$params
  tlist.diff <- tlist[2] - tlist[1]
  tlist <- tlist/tlist.diff
  bSpl <- create.bspline.basis(c(tlist[1],tlist[M]),N)
  return(eval.basis(tlist[train],bSpl))
}

