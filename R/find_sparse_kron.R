#' Implements the fast Laplace Kronecker algorithm.
#'
#' Given samples from some signal and a measurement matrix with Kronecker
#' structure, find the sparse respresentation of the signal using Kronecker
#' compressive sensing. Note that this function is a convenient wrapper for the
#' function \code{\link{FastLaplaceKron}}.
#'
#' The fast Laplace Kronecker algorithm is a method used to solve the compressive
#' sensing problem when the signal is of a large dimension. This typically occurs
#' when working with images where even a relatively small image can require a large
#' measurement matrix when the image is flattened. For example, take the
#' \eqn{M x M} image \eqn{Y}. Let \eqn{y} equal the flattened \eqn{Y} image, with
#' dimensions \eqn{M^2 x 1}. Using the original fast Laplace algorithm would
#' require a \eqn{N x M^2} measurement matrix PHI. However, if it is
#' assumed that the measurement matrix PHI has a Kronecker structure made up of
#' matrices A, B then we can write:
#' \deqn{PHI y = kron(B,A) y = AYB'}
#'
#' which provides a more simple calculation that saves time and memory.
#' See [1] for details.
#'
#' @param A left measurement matrix.
#' @param B right measurement matrix.
#' @param y sample from original signal.
#' @param matrix.return whether or not to return sparse signal in matrix form.
#' @param eta tolerance level in determining convergence of marginal likelihood.
#' @param roundit whether or not to round the marginal likelihood, in order to
#'       avoid machine precision error when comparing across platforms.
#' @param verbose print to screen which basis are added, re-estimated, or deleted.
#' @return The sparse signal as found by the fast Laplace algorithm.
#' @references [1] Cesar F. Caiafa and Andrzej Cichocki, "Computing Sparse
#' Representations of Multidimensional Signals Using Kronecker Bases," in Neural
#' Computation, vol. 25, no. 1, pp. 186-220, 2013.
#' @export
FindSparseKron <- function(A, B, y, matrix.return = T, eta = 1e-8, roundit = FALSE,
                           verbose = FALSE){
  M <- length(y)
  N <- dim(B)[2] * dim(A)[2]
  fl <- FastLaplaceKron(A, B, y, stats::sd(y) ^ 2 / M, eta, roundit = roundit,
                    verbose = verbose)
  x.lap <- rep(0, N)
  x.lap[fl[[2]]] <- fl[[1]]
  if(matrix.return == T){
    x.lap <- matrix(x.lap, dim(A)[2], dim(B)[2])
  }
  return(x.lap)
}
