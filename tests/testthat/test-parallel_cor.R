context("`parallel_cor` computes correlation")

test_that("The result of `parallel_cor` is nearly the same as `cor` ", {
  x <- matrix(rnorm(1000), 100, 10)

  expect_equal(cytostats::parallel_cor(x),
               cor(x))

})
