#' Compute correlation function in batches, in parallel
#'
#' @param x data matrix
#' @param splits shows the number of splits of data that we want to combine  (should be of 2^n form)
#' @param cores number of cores that will be used for the parallel computation
#' @param cov_fun covariance function used on each split
#'
#' @export
#'
#' @importFrom magrittr %>%
#' @importFrom magrittr %<>%
#' @importFrom foreach %dopar%
#'
parallel_cor <- function(x, splits = 2, cores = 2, cov_fun = "two_pass_multi_covar") {

  cov_mat <-as.matrix(parallel_cov(x, splits = splits, cores = cores, cov_fun = cov_fun))

  # \cor(X, Y) = \frac{\cov(X, Y)}{\sd(X) \sd(Y)}

  inv_sd <- as.matrix(diag(diag(cov_mat)^-0.5))

  result <- inv_sd %*% cov_mat %*% inv_sd

  rownames(result) <- colnames(x)

  colnames(result) <- colnames(x)

  result
}