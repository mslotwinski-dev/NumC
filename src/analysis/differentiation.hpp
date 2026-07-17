#pragma once

#include "../common/vector.hpp"
#include "../common/function.hpp"
#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {
namespace analysis {

/// @addtogroup analysis
/// @{

/// @brief Performs Richardson extrapolation for a given function g(h).
template <typename T = double>
T richardson_extrap(numc::func<T> g, T h1, T h2, int p) {
  T ratio = h1 / h2;
  T term = std::pow(ratio, p);
  return (term * g(h2) - g(h1)) / (term - T(1.0));
}

/// @brief Computes the forward difference approximation of the first derivative.
template <typename T = double>
T derivative_forward(numc::func<T> f, T x, T h = T(1e-5)) {
  return (f(x + h) - f(x)) / h;
}

/// @brief Computes the backward difference approximation of the first derivative.
template <typename T = double>
T derivative_backward(numc::func<T> f, T x, T h = T(1e-5)) {
  return (f(x) - f(x - h)) / h;
}

/// @brief Computes the central difference approximation of the first derivative (O(h^2)).
template <typename T = double>
T derivative_central(numc::func<T> f, T x, T h = T(1e-5)) {
  return (f(x + h) - f(x - h)) / (T(2.0) * h);
}

/// @brief Computes the higher-order central difference approximation of the first derivative (O(h^4)).
template <typename T = double>
T derivative_central_O4(numc::func<T> f, T x, T h = T(1e-5)) {
  return (-f(x + T(2.0) * h) + T(8.0) * f(x + h) - T(8.0) * f(x - h) + f(x - T(2.0) * h)) / (T(12.0) * h);
}

/// @brief Computes the central difference approximation of the second derivative (O(h^2)).
template <typename T = double>
T derivative_second_central(numc::func<T> f, T x, T h = T(1e-5)) {
  return (f(x + h) - T(2.0) * f(x) + f(x - h)) / (h * h);
}

/// @brief Computes the higher-order central difference approximation of the second derivative (O(h^4)).
template <typename T = double>
T derivative_second_central_O4(numc::func<T> f, T x, T h = T(1e-5)) {
  return (-f(x + T(2.0) * h) + T(16.0) * f(x + h) - T(30.0) * f(x) + T(16.0) * f(x - h) - f(x - T(2.0) * h)) / (T(12.0) * h * h);
}

/// @brief Generates a function representing the first derivative of f.
template <typename T = double>
numc::func<T> derivative(numc::func<T> f, T h = T(1e-5), int order = 4) {
  if (order == 4) {
    return numc::func<T>([f, h](T x) { return derivative_central_O4(f, x, h); });
  }
  return numc::func<T>([f, h](T x) { return derivative_central(f, x, h); });
}

/// @brief Generates a function representing the second derivative of f.
template <typename T = double>
numc::func<T> derivative2(numc::func<T> f, T h = T(1e-5), int order = 4) {
  if (order == 4) {
    return numc::func<T>([f, h](T x) { return derivative_second_central_O4(f, x, h); });
  }
  return numc::func<T>([f, h](T x) { return derivative_second_central(f, x, h); });
}

/// @brief Computes the first derivative of a polynomial with coefficients a at x.
template <typename T = double>
T diff_poly(const vector<T>& a, T x) {
  size_t n = a.size() - 1;
  T p = a[n];
  T dp = T(0.0);
  for (size_t i = 1; i <= n; ++i) {
    dp = dp * x + p;
    p = p * x + a[n - i];
  }
  return dp;
}

/// @brief Computes the derivative (order 1 or 2) of a cubic spline at x.
template <typename T = double>
T diff_spline(const vector<T>& xData, const vector<T>& yData, const vector<T>& k, T x, int order = 1) {
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

  if (order == 1) {
    T term1 = (std::pow(x - xData[i + 1], 2) / (T(2.0) * h) - h / T(6.0)) * k[i];
    T term2 = (std::pow(x - xData[i], 2) / (T(2.0) * h) - h / T(6.0)) * k[i + 1];
    T term3 = (yData[i] - yData[i + 1]) / h;
    return term1 - term2 + term3;
  } else if (order == 2) {
    return k[i] * (x - xData[i + 1]) / h - k[i + 1] * (x - xData[i]) / h;
  }
  return T(0.0);
}

/// @}

}  // namespace analysis
}  // namespace numc
