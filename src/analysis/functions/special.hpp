#pragma once

#include "../../common/function.hpp"
#include "../../inc.hpp"
#include "../../utility/log.hpp"

namespace numc {

/// @defgroup special Special Functions
/// @ingroup functions
/// @brief Advanced mathematical functions including error functions, Gamma/Beta functions, elliptic integrals, and the Riemann zeta function.
/// @{

/// @brief Computes the Gauss error function, $\operatorname{erf}(x)$.
/// Widely used in probability, statistics, and partial differential equations.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing $\operatorname{erf}(f(x))$.
template <typename T>
numc::func<T> erf(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::erf(f(x)); });
}

/// @brief Computes the complementary error function, $\operatorname{erfc}(x) = 1 - \operatorname{erf}(x)$.
/// Useful for avoiding loss of precision when evaluating for large positive values of x.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing $\operatorname{erfc}(f(x))$.
template <typename T>
numc::func<T> erfc(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::erfc(f(x)); });
}

/// @brief Computes the Gamma function (true gamma), $\Gamma(x)$.
/// An extension of the factorial function to complex and real number arguments, shifted by 1: $\Gamma(n) = (n-1)!$.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing $\Gamma(f(x))$.
template <typename T>
numc::func<T> tgamma(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::tgamma(f(x)); });
}

/// @brief Computes the natural logarithm of the absolute value of the Gamma function, $\ln(|\Gamma(x)|)$.
/// Essential for avoiding floating-point overflow when working with combinatorial formulas or probability distributions.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing $\ln(|\Gamma(f(x))|)$.
template <typename T>
numc::func<T> lgamma(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::lgamma(f(x)); });
}

/// @brief Computes the Euler Beta function, $B(x, y) = \frac{\Gamma(x)\Gamma(y)}{\Gamma(x+y)}$.
/// @tparam T Floating-point type.
/// @param f1 The first input expression ($x$).
/// @param f2 The second input expression ($y$).
/// @return An expression representing the Beta function evaluated at $f_1(x)$ and $f_2(x)$.
template <typename T>
numc::func<T> beta(const numc::func<T>& f1, const numc::func<T>& f2) {
  return numc::func<T>([f1, f2](T x) { return std::beta(f1(x), f2(x)); });
}

/// @brief Computes the complete elliptic integral of the first kind, $K(k)$.
/// @tparam T Floating-point type.
/// @param k The elliptic modulus expression ($k \in [-1, 1]$).
/// @return An expression representing $K(k(x))$.
template <typename T>
numc::func<T> comp_ellint_1(const numc::func<T>& k) {
  return numc::func<T>([k](T x) { return std::comp_ellint_1(k(x)); });
}

/// @brief Computes the complete elliptic integral of the second kind, $E(k)$.
/// Represents the circumference of an ellipse.
/// @tparam T Floating-point type.
/// @param k The elliptic modulus expression ($k \in [-1, 1]$).
/// @return An expression representing $E(k(x))$.
template <typename T>
numc::func<T> comp_ellint_2(const numc::func<T>& k) {
  return numc::func<T>([k](T x) { return std::comp_ellint_2(k(x)); });
}

/// @brief Computes the complete elliptic integral of the third kind, $\Pi(\nu, k)$.
/// @tparam T Floating-point type.
/// @param nu The characteristic of the elliptic integral.
/// @param k The elliptic modulus expression ($k \in [-1, 1]$).
/// @return An expression representing $\Pi(\nu, k(x))$.
template <typename T>
numc::func<T> comp_ellint_3(T nu, const numc::func<T>& k) {
  return numc::func<T>([nu, k](T x) { return std::comp_ellint_3(nu, k(x)); });
}

/// @brief Computes the incomplete elliptic integral of the first kind, $F(k, \phi)$.
/// @tparam T Floating-point type.
/// @param k The elliptic modulus expression.
/// @param phi The amplitude expression (in radians).
/// @return An expression representing $F(k(x), \phi(x))$.
template <typename T>
numc::func<T> ellint_1(const numc::func<T>& k, const numc::func<T>& phi) {
  return numc::func<T>([k, phi](T x) { return std::ellint_1(k(x), phi(x)); });
}

/// @brief Computes the incomplete elliptic integral of the second kind, $E(k, \phi)$.
/// @tparam T Floating-point type.
/// @param k The elliptic modulus expression.
/// @param phi The amplitude expression (in radians).
/// @return An expression representing $E(k(x), \phi(x))$.
template <typename T>
numc::func<T> ellint_2(const numc::func<T>& k, const numc::func<T>& phi) {
  return numc::func<T>([k, phi](T x) { return std::ellint_2(k(x), phi(x)); });
}

/// @brief Computes the incomplete elliptic integral of the third kind, $\Pi(\nu, k, \phi)$.
/// @tparam T Floating-point type.
/// @param nu The characteristic parameter.
/// @param k The elliptic modulus expression.
/// @param phi The amplitude expression (in radians).
/// @return An expression representing $\Pi(\nu, k(x), \phi(x))$.
template <typename T>
numc::func<T> ellint_3(T nu, const numc::func<T>& k, const numc::func<T>& phi) {
  return numc::func<T>([nu, k, phi](T x) { return std::ellint_3(nu, k(x), phi(x)); });
}

/// @brief Computes the exponential integral, $\operatorname{Ei}(x)$.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing $\operatorname{Ei}(f(x))$.
template <typename T>
numc::func<T> expint(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::expint(f(x)); });
}

/// @brief Computes the Riemann zeta function, $\zeta(s)$.
/// Plays a pivotal role in analytical number theory and physics.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing $\zeta(f(x))$.
template <typename T>
numc::func<T> riemann_zeta(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::riemann_zeta(f(x)); });
}

/// @} // End of special group

}  // namespace numc