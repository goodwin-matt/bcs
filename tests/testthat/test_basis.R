context("Testing basis functions")

test_that("the WaveletBasis, FourierBasis, and BSplineBasis functions produce
          the correct matrices", {
  b.wav.test <- data.matrix(read.csv("bWavTest.csv", header = TRUE))
  b.wav.null.test <- data.matrix(read.csv("bWavNULLTest.csv", header = TRUE))
  b.for.test <- data.matrix(read.csv("bForTest.csv", header = TRUE))
  b.for.null.test <- data.matrix(read.csv("bForNULLTest.csv", header = TRUE))
  b.spline.test <- data.matrix(read.csv("bSplineTest.csv", header = TRUE))
  b.spline.null.test <- data.matrix(read.csv("bSplineNULLTest.csv",
                                             header = TRUE))
  dimnames(b.wav.test) <-  NULL
  dimnames(b.wav.null.test) <- NULL
  expect_equal(WaveletBasis(8, seq(1:5)), b.wav.test)
  expect_equal(WaveletBasis(8), b.wav.null.test)
  expect_equal(FourierBasis(1:10, 11, 1:5), b.for.test)
  expect_equal(FourierBasis(1:10, 11), b.for.null.test)
  expect_equal(BSplineBasis(1:10, 11, 1:5), b.spline.test)
  expect_equal(BSplineBasis(1:10, 11), b.spline.null.test)

})
