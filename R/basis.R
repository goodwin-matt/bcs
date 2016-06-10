#' Finds the discrete wavelet transform matrix.
#'
#' Uses the package \code{wmtsa} to find the matrix of a wavelet basis. Each
#' column is a basis function from the given basis, evaluated at particular points
#' along the rows of the matrix.
#'
#' @param signal original signal to find the wavelet transform for. Must be the
#'        same length as the number of basis.
#' @param N number of wavelet basis to keep.
#' @param train the rows to keep of the matrix.
#' @param wavelet the type of wavelet basis to use. See the documentation for the
#'        package \code{wmtsa} for more information.
#' @return A PxN discrete wavelet transform matrix, where P is equal to the
#'        length of \code{train} and N is the number of basis.
#' @export
WaveletBasis <- function(signal, N, train = NULL, wavelet = 'Haar'){
  P <- length(signal)
  # Default for training points is to use all points
  if(is.null(train)){
    train <- 1:P
  }
  M <- length(train)
  w <- wmtsa::wavDWT(signal,wavelet=wavelet)
  lvl <- w$dictionary$n.levels
  basis <- matrix(NA,M,N)
  compt <- 1
  w$data[lvl+1][[1]] <- 0
  for (i in 1:lvl){
    w$data[[i]]<-rep(0,2^(lvl-i))
  }
  for (i in 1:(lvl+1)){
    for (j in 1:2^(lvl-i)){
      w$data[[i]][j] <- 1
      basis[,compt] <- wmtsa::reconstruct(w)[train]
      w$data[[i]][j] <- 0
      compt <- compt + 1
    }
  }
  return(basis)
}

#' Finds the discrete Fourier transformation matrix.
#'
#' Uses the package \code{fda} to find a matrix where the columns are the different
#' basis and the rows are where the basis are evaluated at.
#'
#' @param tlist an array of the specific points where the basis are evaluated at.
#' @param N number of basis in the matrix.
#' @param train indices corresponding to which rows of the matrix to keep.
#' @return A PxN discrete Fourier transformation matrix where P is equal to the
#'        length of \code{train} and N is the number of basis.
#' @export
FourierBasis <- function(tlist, N, train = NULL){
  M <- length(tlist)
  if(is.null(train)){
    train <- 1:M
  }
  # Forces the basis to be even
  if(N%%2==0)
    suppressWarnings(bFor <- fda::create.fourier.basis(c(tlist[1],tlist[M]),N,dropind=N-1))
  else
    bFor <- fda::create.fourier.basis(c(tlist[1],tlist[M]),N)
  return(fda::eval.basis(tlist[train],bFor))
}

#' Finds the transformation matrix for the B-spline basis.
#'
#' Uses the package \code{fda} to find a matrix where each column is a basis
#' function from the B-spline basis, evaluated at particular points along the
#' rows of the matrix.
#'
#' @param tlist an array of the specific points where the basis are evaluated at.
#' @param N number of basis in the matrix.
#' @param train indices corresponding to which rows of the matrix to keep.
#' @return A PxN matrix where P is equal to the length of \code{train} and N is
#'        the number of basis.
#' @export
BSplineBasis <- function(tlist, N, train = NULL){
  M <- length(tlist)
  if(is.null(train)){
    train <- 1:M
  }
  # The first two lines below shift the time scale so that the difference between
  # knots is roughly equal to 1. To see where the knots are placed take
  # bSpl$params
  tlist.diff <- tlist[2] - tlist[1]
  tlist <- tlist/tlist.diff
  bSpl <- fda::create.bspline.basis(c(tlist[1],tlist[M]),N)
  return(fda::eval.basis(tlist[train],bSpl))
}

