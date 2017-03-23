context("Testing support")

test_that("the intersection support function is correct", {
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

test_that("the set difference support function is correct", {
  expect_equal(setdiff(matrix(c(3), 1, 1), matrix(c(2, 3), 2, 1)),
               matrix(0, 0, 0))
  expect_equal(setdiff(matrix(c(2, 3), 2, 1), matrix(c(2, 3), 2, 1)),
               matrix(0, 0, 0))
  expect_equal(setdiff(matrix(c(2, 3), 2, 1), matrix(c(4, 5), 2, 1)),
               matrix(c(2, 3), 2, 1))
  expect_equal(setdiff(matrix(c(2, 3), 2, 1), matrix(2, 1, 1)), matrix(3, 1, 1))
})

