#pragma once

#include "../../common/function.hpp"
#include "../../inc.hpp"
#include "../../utility/log.hpp"

namespace numc {

/// @defgroup analysis Mathematical Analysis
/// @brief Mathematical analysis module, operations on continuous functions, calculus, and signal processing.

/// @defgroup functions Functions
/// @ingroup analysis
/// @brief Bunch of basic mathematical functions that can be applied to numeric calculations.

/// @defgroup elementary Elementary and Core Functions
/// @ingroup functions
/// @brief Basic mathematical operations: trigonometry, logarithms, hyperbolic functions, and rounding.
/// @{

/// @brief Computes the exponential function (e^x) for a given expression tree.
/// @tparam T Floating-point type (e.g., double).
/// @param f The input expression.
/// @return An expression representing e^f(x).
template <typename T>
numc::func<T> exp(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::exp(f(x)); });
}

/// @brief Computes the natural logarithm (base e) for a given expression tree.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing ln(f(x)).
template <typename T>
numc::func<T> ln(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::log(f(x)); });
}

/// @brief Computes the binary logarithm (base 2).
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing log2(f(x)).
template <typename T>
numc::func<T> log2(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::log2(f(x)); });
}

/// @brief Computes the common logarithm (base 10).
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing log10(f(x)).
template <typename T>
numc::func<T> log10(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::log10(f(x)); });
}

/// @brief Computes the logarithm with an arbitrary base.
/// Uses the change-of-base formula: log_a(x) = ln(x) / ln(a).
/// @tparam T Floating-point type.
/// @param base The base of the logarithm.
/// @param f The input expression.
/// @return An expression representing log_base(f(x)).
template <typename T>
numc::func<T> log(T base, const numc::func<T>& f) {
  return numc::func<T>([base, f](T x) { return std::log(f(x)) / std::log(base); });
}

/// @brief Computes the sine of a given expression (in radians).
/// @tparam T Floating-point type.
/// @param f The input expression (angle in radians).
/// @return An expression representing sin(f(x)).
template <typename T>
numc::func<T> sin(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::sin(f(x)); });
}

/// @brief Computes the cosine of a given expression (in radians).
/// @tparam T Floating-point type.
/// @param f The input expression (angle in radians).
/// @return An expression representing cos(f(x)).
template <typename T>
numc::func<T> cos(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::cos(f(x)); });
}

/// @brief Computes the tangent of a given expression (in radians).
/// @tparam T Floating-point type.
/// @param f The input expression (angle in radians).
/// @return An expression representing tan(f(x)).
template <typename T>
numc::func<T> tan(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::tan(f(x)); });
}

/// @brief Computes the cotangent of a given expression (in radians).
/// @tparam T Floating-point type.
/// @param f The input expression (angle in radians).
/// @return An expression representing cot(f(x)).
template <typename T>
numc::func<T> cot(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return T(1.0) / std::tan(f(x)); });
}

/// @brief Computes the principal value of the arc sine.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing asin(f(x)) in radians.
template <typename T>
numc::func<T> asin(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::asin(f(x)); });
}

/// @brief Computes the principal value of the arc cosine.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing acos(f(x)) in radians.
template <typename T>
numc::func<T> acos(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::acos(f(x)); });
}

/// @brief Computes the principal value of the arc tangent.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing atan(f(x)) in radians.
template <typename T>
numc::func<T> atan(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::atan(f(x)); });
}

/// @brief Computes the arc cotangent.
/// Uses a numerically safe atan2 transformation.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing acot(f(x)) in radians.
template <typename T>
numc::func<T> acot(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::atan2(T(1.0), f(x)); });
}

/// @brief Computes the hyperbolic sine.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing sinh(f(x)).
template <typename T>
numc::func<T> sinh(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::sinh(f(x)); });
}

/// @brief Computes the hyperbolic cosine.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing cosh(f(x)).
template <typename T>
numc::func<T> cosh(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::cosh(f(x)); });
}

/// @brief Computes the hyperbolic tangent.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing tanh(f(x)).
template <typename T>
numc::func<T> tanh(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::tanh(f(x)); });
}

/// @brief Computes the hyperbolic cotangent.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing coth(f(x)).
template <typename T>
numc::func<T> coth(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return T(1.0) / std::tanh(f(x)); });
}

/// @brief Computes the inverse hyperbolic sine.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing asinh(f(x)).
template <typename T>
numc::func<T> asinh(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::asinh(f(x)); });
}

/// @brief Computes the inverse hyperbolic cosine.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing acosh(f(x)).
template <typename T>
numc::func<T> acosh(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::acosh(f(x)); });
}

/// @brief Computes the inverse hyperbolic tangent.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing atanh(f(x)).
template <typename T>
numc::func<T> atanh(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::atanh(f(x)); });
}

/// @brief Computes the inverse hyperbolic cotangent.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing acoth(f(x)).
template <typename T>
numc::func<T> acoth(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::atanh(T(1.0) / f(x)); });
}

/// @brief Computes the square root of an expression.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing sqrt(f(x)).
template <typename T>
numc::func<T> sqrt(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::sqrt(f(x)); });
}

/// @brief Computes the cubic root of an expression.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing cbrt(f(x)).
template <typename T>
numc::func<T> cbrt(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::cbrt(f(x)); });
}

/// @brief Computes the arbitrary nth root of an expression.
/// @tparam T Floating-point type.
/// @param degree The degree of the root (e.g., 4.0 for the fourth root).
/// @param f The input expression.
/// @return An expression representing f(x)^(1/degree).
template <typename T>
numc::func<T> root(T degree, const numc::func<T>& f) {
  return numc::func<T>([degree, f](T x) { return std::pow(f(x), T(1.0) / degree); });
}

/// @brief Returns the absolute value (modulus) of an expression.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing |f(x)|.
template <typename T>
numc::func<T> abs(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::abs(f(x)); });
}

/// @brief Returns the sign of an expression (signum function).
/// Returns 1 for positive values, -1 for negative values, and 0 for exact zero.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing sgn(f(x)).
template <typename T>
numc::func<T> sign(const numc::func<T>& f) {
  return numc::func<T>([f](T x) {
    T val = f(x);
    if (val > T(0)) return T(1);
    if (val < T(0)) return T(-1);
    return T(0);
  });
}

/// @brief Rounds the expression down to the nearest integer.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing floor(f(x)).
template <typename T>
numc::func<T> floor(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::floor(f(x)); });
}

/// @brief Rounds the expression up to the nearest integer.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing ceil(f(x)).
template <typename T>
numc::func<T> ceil(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::ceil(f(x)); });
}

/// @brief Rounds the expression to a specified number of decimal places.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @param places The number of decimal places (default is 0, which rounds to the nearest integer).
/// @return An expression with rounded values.
template <typename T>
numc::func<T> round(const numc::func<T>& f, int places = 0) {
  return numc::func<T>([f, places](T x) {
    T factor = std::pow(T(10.0), T(places));
    return std::round(f(x) * factor) / factor;
  });
}

/// @} // End of elementary group

}  // namespace numc