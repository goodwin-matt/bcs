#' Implements the Fast Laplace Algorithm
#'
#' Implements the fast Laplace algorithm given a measurement matrix and samples
#' from a signal or function.  Note that this function is a convenient wrapper
#' for the function \code{\link{FastLaplace}}.
#'
#' This code implements the fast Laplace algorithm.
#' The fast Laplace algorithm is a method
#' used to solve the compressive sensing problem, or in general, a highly
#' underdetermined system of equations. It does this by taking the
#' system of equations
#' \deqn{y = \Phi w + n}
#' and converting it into a minimization problem
#' where we minimize the error with a constraint on \eqn{w}
#' (the vector we are solving for) that enforces
#' sparsity. The fast Laplace method uses a Bayesian framework, and in
#' particular, uses a Laplace prior to enforce sparsity on \eqn{w}.
#' See [1] for more information.
#'
#' @param PHI typically equals the product of a measurment matrix and basis
#' representation matrix, such as the wavelet basis.
#' The solution vector \eqn{w} (see Details below) is assumed to be sparse in the chosen basis.
#' @param y CS measurements, samples from the signal or function.
#' @param eta tolerance level in determining convergence of marginal likelihood.
#' @param roundit whether or not to round the marginal likelihood, in order to
#'       avoid machine precision error when comparing across platforms.
#' @param verbose print to screen which basis are added, re-estimated, or deleted.
#' @return The sparse signal \eqn{w} as found by the fast Laplace algorithm.
#' @references [1] S. D. Babacan, R. Molina and A. K. Katsaggelos, "Bayesian
#' Compressive Sensing Using Laplace Priors," in IEEE Transactions on Image
#' Processing, vol. 19, no. 1, pp. 53-63, Jan. 2010.
#' @examples
#' # size of the basis function expansion
#' N <- 64
#'
#' # generate sparse coefficient vector
#' w <- rep(0,N)
#' w[sample(1:N,10)] <- runif(10,-1,1)
#'
#' # create wavelet basis trasform matrix
#' wavelet.basis <- WaveletBasis(N)
#' # generate actual signal
#' signal <- wavelet.basis%*%w
#'
#' # now we try and recover 'w' and 'signal' from samples
#' num_samps <- 25
#'
#' # create random measurement matrix
#' measure.mat <- matrix(runif(num_samps*N),num_samps,N)
#' measure.mat <- measure.mat/matrix(rep(sqrt(apply(measure.mat^2,2,sum)),
#'                                       num_samps),num_samps,N,byrow=TRUE);
#'
#' PHI <- measure.mat%*%wavelet.basis
#' # actual samples we see
#' y <- measure.mat%*%signal
#'
#' # use fast Laplace algorithm
#' w_est <- FindSparse(PHI, y)
#'
#' # compare plots of the sparse vector and the estimated sparse vector
#' plot(w,type='h')
#' lines(w_est,type='h',col='red')
#'
#' # estimate signal
#' signal_est <- wavelet.basis%*%w_est
#'
#' # Root mean squared error of estimate
#' error <- sqrt(mean((signal - signal_est)^2))
#' @export
FindSparse <- function(PHI, y, eta = 1e-8, roundit = FALSE, verbose=FALSE){
  M <- dim(PHI)[1]
  N <- dim(PHI)[2]
  fl <- bcs::FastLaplace(PHI, y, stats::sd(y) ^ 2 / M, eta, roundit = roundit,
                    verbose = verbose)
  x.lap <- rep(0, N)
  x.lap[fl[[2]]] <- fl[[1]]
  return(x.lap)
}
