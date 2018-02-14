#' Compute covariance function in batches, in parallel
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
parallel_cov <- function(x, splits = 2, cores = 2, cov_fun = "two_pass_multi_covar") {

  x <- as.matrix(x)

  n <- NROW(x)

  doParallel::registerDoParallel(cores = cores)

  batches <- rep(round(n/splits), splits - 1)

  batches <- c(batches, n - sum(batches))

  j <- 0 # to avoid warning no visible binding for global variable j

  result <- foreach::foreach(j = 1:length(batches)) %dopar% {
    if (j == 1) {
      i <- 1

    } else {
      i <- sum(batches[1:(j-1)]) + 1

    }

    k <- sum(batches[1:j])

    batch_mean <- apply(x[i:k,], 2, mean)

    batch_cov <- do.call(cov_fun, list(x[i:k, ]))

    rbind(batch_mean, batch_cov)
  }

  batch_mean_cov <- do.call(cbind, result)

  mean_cov <- combine_cov_estimates(batch_mean_cov, batches)

  result <- mean_cov[2:NROW(mean_cov), ]

  rownames(result) <- colnames(x)

  colnames(result) <- colnames(x)

  result
}