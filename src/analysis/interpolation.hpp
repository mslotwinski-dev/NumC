#pragma once

#include "../common/vector.hpp"
#include "../common/tensor.hpp"
#include "../linear_algebra/solvers.hpp"
#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {
namespace analysis {

/// @addtogroup analysis
/// @{

/// @brief Evaluates the Lagrange interpolating polynomial at x.
template <typename T = double>
T lagrange_interp(const vector<T>& xData, const vector<T>& yData, T x) {
  size_t n = xData.size();
  T sum = T(0.0);
  for (size_t i = 0; i < n; ++i) {
    T term = yData[i];
    for (size_t j = 0; j < n; ++j) {
      if (j != i) {
        term *= (x - xData[j]) / (xData[i] - xData[j]);
      }
    }
    sum += term;
  }
  return sum;
}

/// @brief Computes the divided difference coefficients for Newton's polynomial.
template <typename T = double>
vector<T> newton_coeffts(const vector<T>& xData, const vector<T>& yData) {
  size_t m = xData.size();
  vector<T> a = yData;
  for (size_t k = 1; k < m; ++k) {
    for (size_t i = m - 1; i >= k; --i) {
      a[i] = (a[i] - a[k - 1]) / (xData[i] - xData[k - 1]);
    }
  }
  return a;
}

/// @brief Evaluates Newton's interpolating polynomial at x.
template <typename T = double>
T newton_eval(const vector<T>& a, const vector<T>& xData, T x) {
  size_t n = xData.size() - 1;
  T p = a[n];
  for (size_t k = 1; k <= n; ++k) {
    p = a[n - k] + (x - xData[n - k]) * p;
  }
  return p;
}

/// @brief Evaluates polynomial interpolant at x using Neville's method.
template <typename T = double>
T neville(const vector<T>& xData, const vector<T>& yData, T x) {
  size_t m = xData.size();
  vector<T> y = yData;
  for (size_t k = 1; k < m; ++k) {
    for (size_t i = 0; i < m - k; ++i) {
      y[i] = ((x - xData[i + k]) * y[i] + (xData[i] - x) * y[i + 1]) / (xData[i] - xData[i + k]);
    }
  }
  return y[0];
}

/// @brief Computes the curvatures (second derivatives) of a natural cubic spline at the knots.
template <typename T = double>
vector<T> cubic_spline_curvatures(const vector<T>& xData, const vector<T>& yData) {
  size_t n = xData.size() - 1;
  vector<T> c(n);
  vector<T> d(n + 1);
  for (size_t i = 0; i < n + 1; ++i) d[i] = T(1.0);
  vector<T> e(n);
  vector<T> k(n + 1);

  for (size_t i = 0; i < n - 1; ++i) {
    c[i] = xData[i] - xData[i + 1];
  }
  for (size_t i = 1; i < n; ++i) {
    d[i] = T(2.0) * (xData[i - 1] - xData[i + 1]);
  }
  for (size_t i = 1; i < n; ++i) {
    e[i] = xData[i] - xData[i + 1];
  }
  for (size_t i = 1; i < n; ++i) {
    T term1 = (yData[i - 1] - yData[i]) / (xData[i - 1] - xData[i]);
    T term2 = (yData[i] - yData[i + 1]) / (xData[i] - xData[i + 1]);
    k[i] = T(6.0) * (term1 - term2);
  }

  auto [c_dec, d_dec, e_dec] = linear_algebra::lu_decomp3(c, d, e);
  k = linear_algebra::lu_solve3(c_dec, d_dec, e_dec, k);
  return k;
}

/// @brief Evaluates the cubic spline at x.
template <typename T = double>
T cubic_spline_eval(const vector<T>& xData, const vector<T>& yData, const vector<T>& k, T x) {
  size_t iLeft = 0;
  size_t iRight = xData.size() - 1;
  while (iRight - iLeft > 1) {
    size_t i = (iLeft + iRight) / 2;
    if (x < xData[i]) {
      iRight = i;
    } else {
      iLeft = i;
    }
  }
  size_t i = iLeft;
  T h = xData[i] - xData[i + 1];
  T term1 = (std::pow(x - xData[i + 1], 3) / h - (x - xData[i + 1]) * h) * k[i] / T(6.0);
  T term2 = (std::pow(x - xData[i], 3) / h - (x - xData[i]) * h) * k[i + 1] / T(6.0);
  T term3 = (yData[i] * (x - xData[i + 1]) - yData[i + 1] * (x - xData[i])) / h;
  return term1 - term2 + term3;
}

/// @brief Fits a polynomial of degree m to the data using least-squares.
template <typename T = double>
vector<T> poly_fit(const vector<T>& xData, const vector<T>& yData, int m) {
  tensor<T> a({static_cast<size_t>(m + 1), static_cast<size_t>(m + 1)});
  vector<T> b(m + 1);
  vector<T> s(2 * m + 1);

  for (size_t i = 0; i < xData.size(); ++i) {
    T temp_b = yData[i];
    for (int j = 0; j <= m; ++j) {
      b[j] += temp_b;
      temp_b *= xData[i];
    }
    T temp_s = T(1.0);
    for (int j = 0; j <= 2 * m; ++j) {
      s[j] += temp_s;
      temp_s *= xData[i];
    }
  }
  for (int i = 0; i <= m; ++i) {
    for (int j = 0; j <= m; ++j) {
      a(i, j) = s[i + j];
    }
  }
  return linear_algebra::gauss_pivot(a, b);
}

/// @brief Computes the standard deviation between the fitted polynomial and the data.
template <typename T = double>
T poly_fit_std_dev(const vector<T>& c, const vector<T>& xData, const vector<T>& yData) {
  size_t n = xData.size() - 1;
  size_t m = c.size() - 1;
  T sigma = T(0.0);
  for (size_t i = 0; i <= n; ++i) {
    T p = c[m];
    for (int j = 0; j < static_cast<int>(m); ++j) {
      p = p * xData[i] + c[m - j - 1];
    }
    sigma += std::pow(yData[i] - p, 2);
  }
  return std::sqrt(sigma / T(n - m));
}

/// @}

}  // namespace analysis
}  // namespace numc
