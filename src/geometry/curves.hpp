#pragma once

#include "spatial.hpp"

namespace numc {
namespace geometry {

/// @addtogroup geometry
/// @{

/// @brief Computes the factorial of a small integer.
inline uint64_t factorial(int n) {
  uint64_t res = 1;
  for (int i = 2; i <= n; ++i) res *= i;
  return res;
}

/// @brief Computes the binomial coefficient C(n, k).
inline uint64_t binomial_coeff(int n, int k) {
  if (k < 0 || k > n) return 0;
  if (k == 0 || k == n) return 1;
  k = std::min(k, n - k); // Take advantage of symmetry
  uint64_t c = 1;
  for (int i = 0; i < k; ++i) {
    c = c * (n - i) / (i + 1);
  }
  return c;
}

/// @brief Evaluates a Bézier curve at parameter t ∈ [0, 1].
/// @param control_points A list of control points.
/// @param t The parameter, typically from 0 to 1.
template <typename T = double>
Point<T> bezier_curve(const std::vector<Point<T>>& control_points, T t) {
  if (control_points.empty()) return {};
  size_t n = control_points.size() - 1;
  
  Point<T> result(control_points[0].size());
  for (size_t i = 0; i < result.size(); ++i) result[i] = T(0.0);

  for (int i = 0; i <= static_cast<int>(n); ++i) {
    T coeff = static_cast<T>(binomial_coeff(n, i)) * std::pow(T(1.0) - t, n - i) * std::pow(t, i);
    result = result + coeff * control_points[i];
  }
  return result;
}

/// @brief Evaluates a B-Spline basis function N_{i,p}(t) using the Cox-de Boor recursion formula.
/// @param i Control point index.
/// @param p Spline degree.
/// @param t Parameter to evaluate at.
/// @param knots Knot vector.
template <typename T = double>
T b_spline_basis(int i, int p, T t, const std::vector<T>& knots) {
  if (p == 0) {
    if (knots[i] <= t && t < knots[i + 1]) return T(1.0);
    return T(0.0);
  }
  
  T left = T(0.0), right = T(0.0);
  
  T denom1 = knots[i + p] - knots[i];
  if (denom1 != T(0.0)) {
    left = ((t - knots[i]) / denom1) * b_spline_basis(i, p - 1, t, knots);
  }
  
  T denom2 = knots[i + p + 1] - knots[i + 1];
  if (denom2 != T(0.0)) {
    right = ((knots[i + p + 1] - t) / denom2) * b_spline_basis(i + 1, p - 1, t, knots);
  }
  
  return left + right;
}

/// @brief Evaluates a B-Spline curve at parameter t.
/// @param control_points A list of control points.
/// @param degree The degree of the B-Spline.
/// @param knots The knot vector.
/// @param t The parameter.
template <typename T = double>
Point<T> b_spline_curve(const std::vector<Point<T>>& control_points, int degree, const std::vector<T>& knots, T t) {
  if (control_points.empty()) return {};
  
  Point<T> result(control_points[0].size());
  for (size_t i = 0; i < result.size(); ++i) result[i] = T(0.0);
  
  for (size_t i = 0; i < control_points.size(); ++i) {
    T basis_val = b_spline_basis(static_cast<int>(i), degree, t, knots);
    result = result + basis_val * control_points[i];
  }
  return result;
}

/// @}

}  // namespace geometry
}  // namespace numc
