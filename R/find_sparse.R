#' Implements the fast Laplace algorithm.
#'
#' Given samples from some signal and a measurement matrix, find the sparse
#' respresentation of the signal using compressive sensing. Note that this
#' function is a convenient wrapper for the function \code{\link{FastLaplace}}.
#'
#' The fast Laplace algorithm is a method used to solve the compressive sensing
#' problem, or in general, a highly underdetermined system of equations. This
#' system can be written out as:
#'   \deqn{y = \Phiw + n}
#' where \eqn{w} is the vector of unknown coefficients to solve for and \eqn{n}
#' is random noise. The method uses a Bayesian framework, and in particular, uses
#' a Laplace prior to incorporate the information that most of the coefficients
#' in the solution vector are zero or close to zero. See [1] for details.
#'
#' @param PHI measurement matrix.
#' @param y sample from original signal.
#' @param eta tolerance level in determining convergence of marginal likelihood.
#' @param roundit whether or not to round the marginal likelihood, in order to
#'       avoid machine precision error when comparing across platforms.
#' @param verbose print to screen which basis are added, re-estimated, or deleted.
#' @return The sparse signal as found by the fast Laplace algorithm.
#' @references [1] S. D. Babacan, R. Molina and A. K. Katsaggelos, "Bayesian
#' Compressive Sensing Using Laplace Priors," in IEEE Transactions on Image
#' Processing, vol. 19, no. 1, pp. 53-63, Jan. 2010.
#' @export
FindSparse <- function(PHI, y, eta = 1e-8, roundit = FALSE, verbose=FALSE){
  M <- dim(PHI)[1]
  N <- dim(PHI)[2]
  fl <- FastLaplace(PHI, y, stats::sd(y) ^ 2 / M, eta, roundit = roundit,
                    verbose = verbose)
  x.lap <- rep(0, N)
  x.lap[fl[[2]]] <- fl[[1]]
  return(x.lap)
}
