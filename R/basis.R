#' Finds the Discrete Wavelet Transform Matrix
#'
#' Uses the functions \code{\link[wmtsa]{wavDWT}} and
#' \code{\link[wmtsa]{reconstruct}} from the \code{wmtsa} package to find the
#' transformation matrix of the given wavelet basis type. Each column of the
#' matrix is a wavelet basis function.
#'
#' @param N number of wavelet basis functions to include in matrix. Note that
#'  N must be a
#' power of 2, otherwise the matrix will include NA's. The reason for this has
#' to do with how the wavelet basis is defined.
#' @param train indices corresponding to which rows of the matrix to keep.
#'        Default is to keep all rows.
#' @param wavelet the type of wavelet basis to use. See
#'        \code{\link[wmtsa]{wavDaubechies}} from \code{wmtsa} for types.
#' @return A P x N discrete wavelet transform matrix, where P is equal to the
#'        length of \code{train} and N is the number of basis. If \code{train}
#'        is \code{NULL} then P equals N.
#' @examples
#' # Find first 8 basis functions of the Haar wavelet type
#' w.Haar <- WaveletBasis(8)
#'
#' # Find first 8 basis functions of the d4 wavelet type, keeping the first
#' # half of the rows
#' w.d4 <- WaveletBasis(8, 1:4, wavelet='d4')
#'
#' @export
WaveletBasis <- function(N, train = NULL, wavelet = "Haar"){
  signal <- 1:N
  P <- length(signal)
  # Default for training points is to use all points
  if (is.null(train)){
    train <- 1:P
  }
  M <- length(train)
  w <- wmtsa::wavDWT(signal, wavelet = wavelet)
  lvl <- w$dictionary$n.levels
  basis <- matrix(NA, M, N)
  compt <- 1
  w$data[lvl + 1][[1]] <- 0
  for (i in 1:lvl){
    w$data[[i]] <- rep(0, 2 ^ (lvl - i))
  }
  for (i in 1:(lvl + 1)){
    for (j in 1:2 ^ (lvl - i)){
      w$data[[i]][j] <- 1
      basis[, compt] <- wmtsa::reconstruct(w)[train]
      w$data[[i]][j] <- 0
      compt <- compt + 1
    }
  }
  return(basis)
}

#' Finds the Discrete Fourier Transformation Matrix
#'
#' Uses the package \code{\link{fda}} to find a transformation matrix where the
#' columns are the different Fourier basis functions, evaluated at specified
#' points along the rows of the matrix.
#'
#' @param tlist an array of the specific points where the basis functions are
#' evaluated at.
#' @param N number of basis functions in the matrix.
#' @param train indices corresponding to which rows of the matrix to keep.
#'        Default is to keep all rows.
#' @return A P x N discrete Fourier transformation matrix where P is equal to the
#'        length of \code{train} and N is the number of basis. If \code{train}
#'        is \code{NULL} then P equals the length of \code{tlist}.
#' @examples
#' # Points to evaluate the basis functions at
#' points <- c(0,1,2,3,4)
#' # Find first 8 Fourier basis functions evaluted at "points"
#' f <- FourierBasis(points, 8)
#'
#' # Find first 8 Fourier basis functions evaluated at "points" but only keep
#' # last two rows
#' f <- FourierBasis(points, 8, train = c(4,5))
#' @export
FourierBasis <- function(tlist, N, train = NULL){
  P <- length(tlist)
  # Default for training points is to use all points
  if (is.null(train)){
    train <- 1:P
  }
  # Forces fda to keep an even number of basis
  if (N %% 2 == 0)
    suppressWarnings(b.for <- fda::create.fourier.basis(c(tlist[1], tlist[P]),
                                                        N, dropind = N - 1))
  else
    b.for <- fda::create.fourier.basis(c(tlist[1], tlist[P]), N)
  return(fda::eval.basis(tlist[train], b.for))
}

#' Finds the Transformation Matrix for the B-spline Basis
#'
#' Uses the package \code{\link{fda}} to find a transformation matrix where each
#' column is a basis function from the B-spline basis, evaluated at specified
#' points along the rows of the matrix.
#'
#' @param tlist an array of the specific points where the basis functions are
#' evaluated at.
#' @param N number of basis functions in the matrix.
#' @param train indices corresponding to which rows of the matrix to keep.
#'        Default is to keep all rows.
#' @return A P x N discrete Fourier transformation matrix where P is equal to the
#'        length of \code{train} and N is the number of basis. If \code{train}
#'        is \code{NULL} then P equals the length of \code{tlist}.
#' @examples
#' # Points to evaluate the basis functions at
#' points <- c(0,1,2,3,4)
#' # Find first 8 B-spline basis functions evaluted at "points"
#' b <- BSplineBasis(points, 8)
#'
#' # Find first 8 B-spline basis functions evaluated at "points" but only keep
#' # last two rows
#' b <- BSplineBasis(points, 8, train = c(4,5))
#' @export
BSplineBasis <- function(tlist, N, train = NULL){
  P <- length(tlist)
  # Default for training points is to use all points
  if (is.null(train)){
    train <- 1:P
  }
  # The first two lines below shift the time scale so that the difference
  # between knots is roughly equal to 1. To see where the knots are placed take
  # bSpl$params
  tlist.diff <- tlist[2] - tlist[1]
  tlist <- tlist / tlist.diff
  b.spl <- fda::create.bspline.basis(c(tlist[1], tlist[P]), N)
  return(fda::eval.basis(tlist[train], b.spl))
}
