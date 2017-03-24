context("Testing the main fast Laplace algorithm")

test_that("the output of the fast Laplace algorithm is correct for small signal"
          , {
  b.for <- data.matrix(read.csv("bForEZ.csv", header = TRUE)[,-1])
  signal <- matrix(read.csv("signalEZ.csv", header = TRUE)[,-1], 100, 1)
  # testing when we round the maximum likelihood to 7 digits
  test.e <- FastLaplace(b.for, signal, sd(signal) ^ 2 / 100, 1e-8, TRUE)
  expect_equal(sum(test.e[[1]]), -0.142390884766597281)
  expect_equal(sum(test.e[[2]]), 4430)
  expect_equal(test.e[[3]][1], 0.40479148230739425696)
  expect_equal(sum(test.e[[4]]), 0.23274684494650607625)
  expect_equal(sum(test.e[[5]]), 1819749.4156093874481)
  expect_equal(length(test.e[[1]]), 87)
  x.lap <- rep(0, 150)
  x.lap[test.e[[2]]] <- test.e[[1]]
  expect_equal(FindSparse(b.for, signal, roundit = TRUE), x.lap)
})

# Using the bForEZ for idx==132 for Matlab and idx==126 and for R count==81
# and idx==65 is where the rounding difference is between R and matlab.

# skip_on_cran()
# test_that("the output of the fast Laplace algorithm is correct for large signal"
#           , {
#   bFor <- data.matrix(read.csv(paste('/Users/goodwinm/MyStuff/Research/Code/',
#                                      'bcs_kron/Test_Files_To_Stay_Local/',
#                                      'bFor2.csv',sep=""),
#                                header=TRUE)[,-1])
#   signal <- matrix(read.csv(paste('/Users/goodwinm/MyStuff/Research/Code/',
#                                   'bcs_kron/Test_Files_To_Stay_Local/',
#                                   'signal2.csv',sep=""),
#                             header=TRUE)[,-1],800,1)
#   # testing when we round the maximum likelihood to 7 digits
#   test.e <- FastLaplace(bFor, signal, sd(signal)^2/800, 1e-8, TRUE)
#   expect_equal(sum(test.e[[1]]),0.095116009295516487643)
#   expect_equal(sum(test.e[[2]]),298117)
#   expect_equal(test.e[[3]][1], 0.011482830173325)
#   expect_equal(sum(test.e[[4]]),0.75847148556943433384)
#   expect_equal(sum(test.e[[5]]),18815188.135733500123)
#   expect_equal(length(test.e[[1]]), 551)
# })



