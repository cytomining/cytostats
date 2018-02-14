context("`parallel_cov` computes covariance")

test_that("The result of `parallel_cov` is nearly the same as `cov` ", {
  x <- matrix(rnorm(1000), 100, 10)

  expect_equal(cytostats::parallel_cov(x),
               cov(x))

})
