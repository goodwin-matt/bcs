#' Implements the fast Laplace algorithm.
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
FindSparse <- function(PHI,y,eta,verbose=FALSE){
  M <- dim(PHI)[1]
  N <- dim(PHI)[2]
  fl <- FastLaplace(PHI, y, sd(y)^2/M, eta, verbose = verbose)
  x.lap <- rep(0,N)
  x.lap[fl[[2]]] <- fl[[1]]
  return(x.lap)
}



