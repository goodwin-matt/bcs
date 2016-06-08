context("Testing basis")

test_that("wavelets, Fourier, B-splines", {
  bWavTest <- data.matrix(read.csv('bWavTest.csv',header=TRUE))
  bWavNULLTest <- data.matrix(read.csv('bWavNULLTest.csv',header=TRUE))
  bForTest <- data.matrix(read.csv('bForTest.csv',header=TRUE))
  bForNULLTest <- data.matrix(read.csv('bForNULLTest.csv',header=TRUE))
  bSplineTest <- data.matrix(read.csv('bSplineTest.csv',header=TRUE))
  bSplineNULLTest <- data.matrix(read.csv('bSplineNULLTest.csv',header=TRUE))
  dimnames(bWavTest) <-  NULL
  dimnames(bWavNULLTest) <- NULL
  expect_equal(WaveletBasis(1:8,8,seq(1:5)), bWavTest)
  expect_equal(WaveletBasis(1:8,8), bWavNULLTest)
  expect_equal(FourierBasis(1:10,11,1:5), bForTest)
  expect_equal(FourierBasis(1:10,11), bForNULLTest)
  expect_equal(BSplineBasis(1:10,11,1:5), bSplineTest)
  expect_equal(BSplineBasis(1:10,11), bSplineNULLTest)

})
