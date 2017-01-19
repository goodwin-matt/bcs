context("Testing support")

test_that("intersection of two sets", {
  #TODO: test and optimize intersection better
  expect_equal(intersect(matrix(c(2, 3), 2, 1), matrix(c(1, 2, 3), 3, 1)),
               matrix(c(2, 1, 3, 2), 2, 2, byrow = TRUE))
  expect_equal(intersect(matrix(c(1, 2, 3), 3, 1), matrix(c(2, 3), 2, 1)),
               matrix(c(2, 0, 3, 1), 2, 2, byrow = TRUE))
  expect_equal(intersect(matrix(c(1, 2, 3), 3, 1), matrix(c(1, 2, 3), 3, 1)),
               matrix(c(1, 0, 2, 1, 3, 2), 3, 2, byrow = TRUE))
  expect_equal(intersect(matrix(c(1, 3, 2), 3, 1), matrix(c(1, 2, 3), 3, 1)),
               matrix(c(1, 0, 2, 1, 3, 2), 3, 2, byrow = TRUE))
})

test_that("set difference of two sets", {
  expect_equal(setdiff(matrix(c(3), 1, 1), matrix(c(2, 3), 2, 1)),
               matrix(0, 0, 0))
  expect_equal(setdiff(matrix(c(2, 3), 2, 1), matrix(c(2, 3), 2, 1)),
               matrix(0, 0, 0))
  expect_equal(setdiff(matrix(c(2, 3), 2, 1), matrix(c(4, 5), 2, 1)),
               matrix(c(2, 3), 2, 1))
  expect_equal(setdiff(matrix(c(2, 3), 2, 1), matrix(2, 1, 1)), matrix(3, 1, 1))
})

test_that("get column of kronecker structure matrix",{
  A <- matrix(1:16, 4, 4)
  B <- matrix(1:100,10,10)
  expect_equal(GetCol(6, A, B), matrix(kronecker(A,B)[,7]))
  expect_equal(GetCol(9, A, B), matrix(kronecker(A,B)[,10]))
  expect_equal(GetCol(15, A, B), matrix(kronecker(A,B)[,16]))
})

test_that("get column sum of squared values of kronecker structure matrix",{
  A <- matrix(1:16, 4, 4)
  B <- matrix(1:100,10,10)
  expect_equal(GetColSumSquared(A, B), matrix(colSums(kronecker(A,B)^2),1,40))
})

test_that("get matrix multiplication of matrix with kronecker structure",{
  A <- matrix(1:100,10,10)[1:4,]
  B <- matrix(sample(1:100,100),10,10)[1:4,]
  x <- matrix(1:100,100,1)
  x1 <- matrix(1:16,16,1)
  expect_equal(MultMatrix(A, B, x, 1), kronecker(B,A)%*%x)
  expect_equal(MultMatrix(A, B, x1, 2), kronecker(t(B),t(A))%*%x1)
})
