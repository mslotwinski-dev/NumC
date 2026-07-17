#pragma once

#include "../common/vector.hpp"
#include "../common/tensor.hpp"
#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {
namespace statistics {

/// @addtogroup statistics
/// @{

/// @brief Computes the arithmetic mean of a vector.
template <typename T = double>
T mean(const vector<T>& x) {
  size_t n = x.size();
  T sum = T(0.0);
  for (size_t i = 0; i < n; ++i) sum += x[i];
  return sum / T(n);
}

/// @brief Computes the geometric mean of a vector (all elements must be positive).
template <typename T = double>
T geometric_mean(const vector<T>& x) {
  size_t n = x.size();
  T log_sum = T(0.0);
  for (size_t i = 0; i < n; ++i) log_sum += std::log(x[i]);
  return std::exp(log_sum / T(n));
}

/// @brief Computes the harmonic mean of a vector (all elements must be positive).
template <typename T = double>
T harmonic_mean(const vector<T>& x) {
  size_t n = x.size();
  T inv_sum = T(0.0);
  for (size_t i = 0; i < n; ++i) inv_sum += T(1.0) / x[i];
  return T(n) / inv_sum;
}

/// @brief Computes the median of a vector.
template <typename T = double>
T median(vector<T> x) {
  size_t n = x.size();
  std::vector<T> sorted(n);
  for (size_t i = 0; i < n; ++i) sorted[i] = x[i];
  std::sort(sorted.begin(), sorted.end());
  if (n % 2 == 0) {
    return (sorted[n / 2 - 1] + sorted[n / 2]) / T(2.0);
  }
  return sorted[n / 2];
}

/// @brief Computes the variance (population) of a vector.
template <typename T = double>
T variance(const vector<T>& x) {
  size_t n = x.size();
  T m = mean(x);
  T sum = T(0.0);
  for (size_t i = 0; i < n; ++i) sum += (x[i] - m) * (x[i] - m);
  return sum / T(n);
}

/// @brief Computes the sample variance (Bessel correction) of a vector.
template <typename T = double>
T sample_variance(const vector<T>& x) {
  size_t n = x.size();
  T m = mean(x);
  T sum = T(0.0);
  for (size_t i = 0; i < n; ++i) sum += (x[i] - m) * (x[i] - m);
  return sum / T(n - 1);
}

/// @brief Computes the standard deviation (population) of a vector.
template <typename T = double>
T std_dev(const vector<T>& x) {
  return std::sqrt(variance(x));
}

/// @brief Computes the sample standard deviation of a vector.
template <typename T = double>
T sample_std_dev(const vector<T>& x) {
  return std::sqrt(sample_variance(x));
}

/// @brief Computes the skewness of a vector.
template <typename T = double>
T skewness(const vector<T>& x) {
  size_t n = x.size();
  T m = mean(x);
  T s = std_dev(x);
  if (s < T(1e-15)) return T(0.0);
  T sum = T(0.0);
  for (size_t i = 0; i < n; ++i) {
    T d = (x[i] - m) / s;
    sum += d * d * d;
  }
  return sum / T(n);
}

/// @brief Computes the kurtosis (excess) of a vector.
template <typename T = double>
T kurtosis(const vector<T>& x) {
  size_t n = x.size();
  T m = mean(x);
  T s = std_dev(x);
  if (s < T(1e-15)) return T(0.0);
  T sum = T(0.0);
  for (size_t i = 0; i < n; ++i) {
    T d = (x[i] - m) / s;
    sum += d * d * d * d;
  }
  return sum / T(n) - T(3.0);
}

/// @brief Computes the covariance between two vectors.
template <typename T = double>
T covariance(const vector<T>& x, const vector<T>& y) {
  size_t n = x.size();
  T mx = mean(x), my = mean(y);
  T sum = T(0.0);
  for (size_t i = 0; i < n; ++i) sum += (x[i] - mx) * (y[i] - my);
  return sum / T(n);
}

/// @brief Computes the Pearson correlation coefficient between two vectors.
template <typename T = double>
T correlation(const vector<T>& x, const vector<T>& y) {
  T sx = std_dev(x), sy = std_dev(y);
  if (sx < T(1e-15) || sy < T(1e-15)) return T(0.0);
  return covariance(x, y) / (sx * sy);
}

/// @brief Computes the covariance matrix for a set of variables (columns of data).
template <typename T = double>
tensor<T> cov_matrix(const std::vector<vector<T>>& data) {
  size_t p = data.size();      // number of variables
  size_t n = data[0].size();   // number of observations
  tensor<T> C({p, p});
  std::vector<T> means(p);
  for (size_t i = 0; i < p; ++i) means[i] = mean(data[i]);

  for (size_t i = 0; i < p; ++i) {
    for (size_t j = i; j < p; ++j) {
      T sum = T(0.0);
      for (size_t k = 0; k < n; ++k) {
        sum += (data[i][k] - means[i]) * (data[j][k] - means[j]);
      }
      C(i, j) = sum / T(n - 1);
      C(j, i) = C(i, j);
    }
  }
  return C;
}

/// @brief Performs simple linear regression y = a + b*x. Returns {a, b, R^2}.
template <typename T = double>
std::tuple<T, T, T> linear_regression(const vector<T>& x, const vector<T>& y) {
  size_t n = x.size();
  T mx = mean(x), my = mean(y);
  T sxx = T(0.0), sxy = T(0.0), syy = T(0.0);
  for (size_t i = 0; i < n; ++i) {
    T dx = x[i] - mx, dy = y[i] - my;
    sxx += dx * dx;
    sxy += dx * dy;
    syy += dy * dy;
  }
  T b = sxy / sxx;
  T a = my - b * mx;
  T r2 = (sxy * sxy) / (sxx * syy);
  return {a, b, r2};
}

/// @brief Generates a vector of n linearly spaced values between a and b.
template <typename T = double>
vector<T> linspace(T a, T b, size_t n) {
  vector<T> result(n);
  if (n == 1) { result[0] = a; return result; }
  T step = (b - a) / T(n - 1);
  for (size_t i = 0; i < n; ++i) result[i] = a + T(i) * step;
  return result;
}

/// @brief Generates a vector of n logarithmically spaced values between 10^a and 10^b.
template <typename T = double>
vector<T> logspace(T a, T b, size_t n) {
  vector<T> result(n);
  if (n == 1) { result[0] = std::pow(T(10.0), a); return result; }
  T step = (b - a) / T(n - 1);
  for (size_t i = 0; i < n; ++i) result[i] = std::pow(T(10.0), a + T(i) * step);
  return result;
}

/// @brief Applies a moving average filter to a signal with window size w.
template <typename T = double>
vector<T> moving_average(const vector<T>& x, size_t w) {
  size_t n = x.size();
  vector<T> result(n);
  for (size_t i = 0; i < n; ++i) {
    size_t start = (i >= w / 2) ? i - w / 2 : 0;
    size_t end = std::min(i + w / 2 + 1, n);
    T sum = T(0.0);
    for (size_t j = start; j < end; ++j) sum += x[j];
    result[i] = sum / T(end - start);
  }
  return result;
}

/// @brief Computes the p-th percentile of a vector.
template <typename T = double>
T percentile(vector<T> x, T p) {
  size_t n = x.size();
  std::vector<T> sorted(n);
  for (size_t i = 0; i < n; ++i) sorted[i] = x[i];
  std::sort(sorted.begin(), sorted.end());
  T rank = (p / T(100.0)) * T(n - 1);
  size_t lo = static_cast<size_t>(rank);
  size_t hi = lo + 1;
  if (hi >= n) return sorted[n - 1];
  T frac = rank - T(lo);
  return sorted[lo] * (T(1.0) - frac) + sorted[hi] * frac;
}

/// @brief Returns the minimum value in a vector.
template <typename T = double>
T min(const vector<T>& x) {
  T m = x[0];
  for (size_t i = 1; i < x.size(); ++i) if (x[i] < m) m = x[i];
  return m;
}

/// @brief Returns the maximum value in a vector.
template <typename T = double>
T max(const vector<T>& x) {
  T m = x[0];
  for (size_t i = 1; i < x.size(); ++i) if (x[i] > m) m = x[i];
  return m;
}

/// @brief Returns the index of the minimum value in a vector.
template <typename T = double>
size_t argmin(const vector<T>& x) {
  size_t idx = 0;
  for (size_t i = 1; i < x.size(); ++i) if (x[i] < x[idx]) idx = i;
  return idx;
}

/// @brief Returns the index of the maximum value in a vector.
template <typename T = double>
size_t argmax(const vector<T>& x) {
  size_t idx = 0;
  for (size_t i = 1; i < x.size(); ++i) if (x[i] > x[idx]) idx = i;
  return idx;
}

/// @brief Performs a 1-sample T-test. Returns the t-statistic.
/// Null hypothesis: the sample mean is equal to the expected mean (mu).
template <typename T = double>
T t_test_1samp(const vector<T>& x, T mu = T(0.0)) {
  size_t n = x.size();
  if (n <= 1) return T(0.0);
  T m = mean(x);
  T s = std_dev(x);
  return (m - mu) / (s / std::sqrt(static_cast<T>(n)));
}

/// @brief Performs an independent 2-sample T-test (Student's t-test assuming equal variances).
/// Returns the t-statistic.
template <typename T = double>
T t_test_ind(const vector<T>& x, const vector<T>& y) {
  size_t nx = x.size(), ny = y.size();
  if (nx <= 1 || ny <= 1) return T(0.0);
  
  T mx = mean(x), my = mean(y);
  T vx = var(x) * T(nx - 1);
  T vy = var(y) * T(ny - 1);
  
  T pooled_var = (vx + vy) / T(nx + ny - 2);
  T standard_error = std::sqrt(pooled_var * (T(1.0)/T(nx) + T(1.0)/T(ny)));
  
  return (mx - my) / standard_error;
}

/// @brief Performs a Chi-Square goodness-of-fit test.
/// Returns the Chi-Square statistic.
template <typename T = double>
T chi_square_test(const vector<T>& observed, const vector<T>& expected) {
  if (observed.size() != expected.size()) return T(0.0);
  T chi_sq = T(0.0);
  for (size_t i = 0; i < observed.size(); ++i) {
    if (expected[i] > T(0.0)) {
      T diff = observed[i] - expected[i];
      chi_sq += (diff * diff) / expected[i];
    }
  }
  return chi_sq;
}

/// @brief Computes the cumulative sum of a vector.
template <typename T = double>
vector<T> cumsum(const vector<T>& x) {
  size_t n = x.size();
  vector<T> result(n);
  result[0] = x[0];
  for (size_t i = 1; i < n; ++i) result[i] = result[i - 1] + x[i];
  return result;
}

/// @brief Computes the cumulative product of a vector.
template <typename T = double>
vector<T> cumprod(const vector<T>& x) {
  size_t n = x.size();
  vector<T> result(n);
  result[0] = x[0];
  for (size_t i = 1; i < n; ++i) result[i] = result[i - 1] * x[i];
  return result;
}

/// @brief Computes the sum of all elements.
template <typename T = double>
T sum(const vector<T>& x) {
  T s = T(0.0);
  for (size_t i = 0; i < x.size(); ++i) s += x[i];
  return s;
}

/// @brief Computes the product of all elements.
template <typename T = double>
T prod(const vector<T>& x) {
  T p = T(1.0);
  for (size_t i = 0; i < x.size(); ++i) p *= x[i];
  return p;
}

/// @}

}  // namespace statistics
}  // namespace numc
