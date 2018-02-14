#include <Rcpp.h>

using namespace Rcpp;

double online_covar(double *x1, double *x2, int size) {
    int n = 0;
    double mean_1 = 0;
    double mean_2 = 0;
    double mean_12 = 0;

    for (int i = 0; i < size; i++) {
        double xi1 = x1[i];
        double xi2 = x2[i];
        n++;
        double delta_1 = (xi1 - mean_1) / n;
        mean_1 += delta_1;
        double delta_2 = (xi2 - mean_2) / n;
        mean_2 += delta_2;
        mean_12 += (n - 1) * delta_1 * delta_2 - mean_12 / n;
    }

    if (n < 2) {
        return 0;
    } else {
        return (mean_12 * n / (n - 1));
    }
}

double two_pass_covar(double *x1, double *x2, double mean_1, double mean_2, int size) {
    int n = size;
    double mean_12 = 0;

    for (int i = 0; i < size; i++) {
        double delta_1 = (x1[i] - mean_1);
        double delta_2 = (x2[i] - mean_2);
        mean_12 += (delta_1 * delta_2) / (n - 1);
    }

    if (n < 2) {
        return 0;
    } else {
        return (mean_12);
    }
}

NumericMatrix combine_covs_base(NumericMatrix mean_cov_1, NumericMatrix mean_cov_2, int b1, int b2) {
    NumericMatrix result(mean_cov_1.nrow(), mean_cov_1.nrow() - 1);
    for (int i = 1; i < mean_cov_1.nrow(); i++) {
        result(0, i - 1) = (mean_cov_1(0, i - 1) * b1 + mean_cov_2(0, i - 1) * b2) / (b1 + b2);
        for (int j = 0; j < mean_cov_1.nrow() - 1; j++) {
            double cx = ((b1 - 1) * mean_cov_1(i, j) + (b2 - 1) * mean_cov_2(i, j) +
                         (mean_cov_1(0, i - 1) - mean_cov_2(0, i - 1)) * (mean_cov_1(0, j) - mean_cov_2(0, j)) * (b1) *
                         (b2) / (b1 + b2)) / (b1 + b2 - 1);
            result(i, j) = cx;
        }
    }

    return result;
}

//' online_covar
//'
//' @param x1 a numeric vector containing samples of the first random variable
//' @param x2 a numeric vector containing corresponding samples of the second random variable
//'
//' @export
//'
// [[Rcpp::export]]
double online_covar(NumericVector x1, NumericVector x2) {
    return (online_covar(&x1[0], &x2[0], x1.size()));

}

//' two_pass_multi_covar
//'
//' @param s a data matrix whose column covariances are sought
//'
//' @export
//'
// [[Rcpp::export]]
NumericMatrix two_pass_multi_covar(NumericMatrix s) {
    NumericMatrix result(s.ncol(), s.ncol());
    double *Means = new double[s.ncol()];
    double *sx = &s[0];
    for (int i = 0; i < s.ncol(); i++) {
        double sm = 0;
        for (int j = 0; j < s.nrow(); j++) {
            sm += s(j, i);
        }
        sm /= s.nrow();
        Means[i] = sm;
    }

    int n = s.nrow();
    for (int i = 0; i < s.ncol(); i++) {
        for (int j = 0; j <= i; j++) {
            double mean_12 = 0;
            double mean_1 = Means[i];
            double mean_2 = Means[j];
            for (int k = 0; k < n; k++) {
                mean_12 += (sx[i * n + k] - mean_1) * (sx[j * n + k] - mean_2) / (n - 1);   // division inside
                // the loop makes it slower
                // but numerically stabler
            }
            result(i, j) = mean_12;
            result(j, i) = result(i, j);
        }
    }
    delete[] Means;
    return result;
}

//' combine_cov_estimates
//'
//' @param batch_mean_cov the matrix which contains estimated means and covariance for each batch of data.
//'         For n variable, and k batches, it is (n+1)*(n.k) size, with the first row being the means and
//'         rest of the rows being the covariance matrices. Covariance matrices were concatenated column-wise
//'         resulting in n.k columns.
//' @param b a vector containing number of samples in each batch of data
//'
//' @export
//'
// [[Rcpp::export]]
NumericMatrix combine_cov_estimates(NumericMatrix batch_mean_cov, NumericVector b) {
    if (b.size() == 2) {
        return combine_covs_base(batch_mean_cov(_, Range(0, batch_mean_cov.ncol() / 2 - 1)),
                                 batch_mean_cov(_, Range(batch_mean_cov.ncol() / 2, batch_mean_cov.ncol() - 1)), b[0],
                                 b[1]);
    } else {
        NumericMatrix c_left = batch_mean_cov(_, Range(0, batch_mean_cov.ncol() / 2 - 1));
        NumericMatrix c_right = batch_mean_cov(_, Range(batch_mean_cov.ncol() / 2, batch_mean_cov.ncol() - 1));
        NumericVector b_left(b.size() / 2);
        for (int i = 0; i < b.size() / 2; i++) {
            b_left[i] = b[i];
        }

        NumericVector b_right(b.size() / 2);
        for (int i = b.size() / 2; i < b.size(); i++) {
            b_right[i - b.size() / 2] = b[i];
        }

        NumericMatrix result_left = combine_cov_estimates(c_left, b_left);
        NumericMatrix result_right = combine_cov_estimates(c_right, b_right);
        int b_left_sum = 0;
        for (int i = 0; i < b_left.size(); i++) {
            b_left_sum += b_left[i];
        }
        int b_right_sum = 0;
        for (int i = 0; i < b_right.size(); i++) {
            b_right_sum += b_right[i];
        }

        return (combine_covs_base(result_left, result_right, b_left_sum, b_right_sum));
    }
}