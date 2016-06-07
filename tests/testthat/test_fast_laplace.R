context("Testing Fast Laplace")

test_that("General output", {
  bFor <- data.matrix(read.csv('bForEZ.csv',header=TRUE)[,-1])
  signal <- matrix(read.csv('signalEZ.csv',header=TRUE)[,-1],100,1)
  # testing when we round the maximum likelihood to 7 digits
  test.e <- FastLaplace(bFor,signal,sd(signal)^2/100,1e-8,FALSE)
  expect_equal(sum(test.e[[1]]),-0.142390884766597281)
  expect_equal(sum(test.e[[2]]),4430)
  expect_equal(test.e[[3]][1],0.40479148230739425696)
  expect_equal(sum(test.e[[4]]),0.23274684494650607625)
  expect_equal(sum(test.e[[5]]),1819749.4156093874481)
  expect_equal(length(test.e[[1]]), 87)
})


# # Using the bForEZ for idx==132 for Matlab and idx==126 and for R count==81 and idx==65
# # is where the rounding difference is between R and matlab.
