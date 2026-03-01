#pragma once

#include "../../common/function.hpp"
#include "../../inc.hpp"
#include "../../utility/log.hpp"

namespace numc {

/// @defgroup orthogonal Orthogonal Polynomials & Special Functions
/// @ingroup functions
/// @brief Cylindrical, spherical functions, and orthogonal polynomials (Bessel, Legendre, Hermite, Chebyshev, etc.).
/// @{

/// @brief Computes the cylindrical Bessel function of the first kind, $J_\nu(x)$.
/// @tparam T Floating-point type.
/// @param nu The order of the Bessel function.
/// @param f The input expression.
/// @return An expression representing $J_\nu(f(x))$.
template <typename T>
numc::func<T> bessel_j(T nu, const numc::func<T>& f) {
  return numc::func<T>([nu, f](T x) { return std::cyl_bessel_j(nu, f(x)); });
}

/// @brief Computes the cylindrical Bessel function of the second kind (Neumann function), $Y_\nu(x)$.
/// @tparam T Floating-point type.
/// @param nu The order of the Bessel function.
/// @param f The input expression.
/// @return An expression representing $Y_\nu(f(x))$.
template <typename T>
numc::func<T> bessel_y(T nu, const numc::func<T>& f) {
  // In the C++ standard, the Y (Neumann) function is named cyl_neumann
  return numc::func<T>([nu, f](T x) { return std::cyl_neumann(nu, f(x)); });
}

/// @brief Computes the regular modified cylindrical Bessel function, $I_\nu(x)$.
/// @tparam T Floating-point type.
/// @param nu The order of the Bessel function.
/// @param f The input expression.
/// @return An expression representing $I_\nu(f(x))$.
template <typename T>
numc::func<T> bessel_i(T nu, const numc::func<T>& f) {
  return numc::func<T>([nu, f](T x) { return std::cyl_bessel_i(nu, f(x)); });
}

/// @brief Computes the irregular modified cylindrical Bessel function, $K_\nu(x)$.
/// @tparam T Floating-point type.
/// @param nu The order of the Bessel function.
/// @param f The input expression.
/// @return An expression representing $K_\nu(f(x))$.
template <typename T>
numc::func<T> bessel_k(T nu, const numc::func<T>& f) {
  return numc::func<T>([nu, f](T x) { return std::cyl_bessel_k(nu, f(x)); });
}

/// @brief Computes the spherical Bessel function of the first kind, $j_n(x)$.
/// @tparam T Floating-point type.
/// @param n The order of the function (must be an integer $\ge 0$).
/// @param f The input expression.
/// @return An expression representing $j_n(f(x))$.
template <typename T>
numc::func<T> sph_bessel_j(unsigned n, const numc::func<T>& f) {
  return numc::func<T>([n, f](T x) { return std::sph_bessel(n, f(x)); });
}

/// @brief Computes the spherical Bessel function of the second kind (spherical Neumann), $y_n(x)$.
/// @tparam T Floating-point type.
/// @param n The order of the function (must be an integer $\ge 0$).
/// @param f The input expression.
/// @return An expression representing $y_n(f(x))$.
template <typename T>
numc::func<T> sph_bessel_y(unsigned n, const numc::func<T>& f) {
  return numc::func<T>([n, f](T x) { return std::sph_neumann(n, f(x)); });
}

/// @brief Computes the spherical Legendre function (spherical harmonic without the azimuthal part).
/// @tparam T Floating-point type.
/// @param l The degree of the spherical harmonic ($l \ge 0$).
/// @param m The order of the spherical harmonic ($0 \le m \le l$).
/// @param theta The polar angle expression (in radians).
/// @return An expression representing $Y_l^m(\theta, 0)$.
template <typename T>
numc::func<T> sph_legendre(unsigned l, unsigned m, const numc::func<T>& theta) {
  return numc::func<T>([l, m, theta](T x) { return std::sph_legendre(l, m, theta(x)); });
}

/// @brief Computes the Legendre polynomial of degree n, $P_n(x)$.
/// @tparam T Floating-point type.
/// @param n The degree of the polynomial.
/// @param f The input expression.
/// @return An expression representing $P_n(f(x))$.
template <typename T>
numc::func<T> legendre(unsigned n, const numc::func<T>& f) {
  return numc::func<T>([n, f](T x) { return std::legendre(n, f(x)); });
}

/// @brief Computes the associated Legendre polynomial, $P_l^m(x)$.
/// @tparam T Floating-point type.
/// @param l The degree of the polynomial.
/// @param m The order of the polynomial.
/// @param f The input expression.
/// @return An expression representing $P_l^m(f(x))$.
template <typename T>
numc::func<T> assoc_legendre(unsigned l, unsigned m, const numc::func<T>& f) {
  return numc::func<T>([l, m, f](T x) { return std::assoc_legendre(l, m, f(x)); });
}

/// @brief Computes the physicist's Hermite polynomial of degree n, $H_n(x)$.
/// @tparam T Floating-point type.
/// @param n The degree of the polynomial.
/// @param f The input expression.
/// @return An expression representing $H_n(f(x))$.
template <typename T>
numc::func<T> hermite(unsigned n, const numc::func<T>& f) {
  return numc::func<T>([n, f](T x) { return std::hermite(n, f(x)); });
}

/// @brief Computes the Laguerre polynomial of degree n, $L_n(x)$.
/// @tparam T Floating-point type.
/// @param n The degree of the polynomial.
/// @param f The input expression.
/// @return An expression representing $L_n(f(x))$.
template <typename T>
numc::func<T> laguerre(unsigned n, const numc::func<T>& f) {
  return numc::func<T>([n, f](T x) { return std::laguerre(n, f(x)); });
}

/// @brief Computes the associated Laguerre polynomial, $L_n^m(x)$.
/// @tparam T Floating-point type.
/// @param n The degree of the polynomial.
/// @param m The order of the polynomial.
/// @param f The input expression.
/// @return An expression representing $L_n^m(f(x))$.
template <typename T>
numc::func<T> assoc_laguerre(unsigned n, unsigned m, const numc::func<T>& f) {
  return numc::func<T>([n, m, f](T x) { return std::assoc_laguerre(n, m, f(x)); });
}

/// @brief Computes the Chebyshev polynomial of the first kind of degree n, $T_n(x)$.
/// Implemented using iterative recurrence relations for computational stability.
/// @tparam T Floating-point type.
/// @param n The degree of the polynomial.
/// @param f The input expression.
/// @return An expression representing $T_n(f(x))$.
template <typename T>
numc::func<T> chebyshev(unsigned n, const numc::func<T>& f) {
  return numc::func<T>([n, f](T x) {
    T val = f(x);
    if (n == 0) return T(1.0);
    if (n == 1) return val;

    T t0 = T(1.0), t1 = val, tn = T(0.0);
    for (unsigned i = 2; i <= n; ++i) {
      tn = T(2.0) * val * t1 - t0;
      t0 = t1;
      t1 = tn;
    }
    return tn;
  });
}

/// @brief Computes the Chebyshev polynomial of the second kind of degree n, $U_n(x)$.
/// Implemented using iterative recurrence relations.
/// @tparam T Floating-point type.
/// @param n The degree of the polynomial.
/// @param f The input expression.
/// @return An expression representing $U_n(f(x))$.
template <typename T>
numc::func<T> chebyshev_u(unsigned n, const numc::func<T>& f) {
  return numc::func<T>([n, f](T x) {
    T val = f(x);
    if (n == 0) return T(1.0);
    if (n == 1) return T(2.0) * val;

    T u0 = T(1.0), u1 = T(2.0) * val, un = T(0.0);
    for (unsigned i = 2; i <= n; ++i) {
      un = T(2.0) * val * u1 - u0;
      u0 = u1;
      u1 = un;
    }
    return un;
  });
}

/// @brief Computes the Gegenbauer (ultraspherical) polynomial, $C_n^{(\alpha)}(x)$.
/// Implemented using iterative recurrence relations.
/// @tparam T Floating-point type.
/// @param n The degree of the polynomial.
/// @param alpha The alpha parameter of the polynomial.
/// @param f The input expression.
/// @return An expression representing $C_n^{(\alpha)}(f(x))$.
template <typename T>
numc::func<T> gegenbauer(unsigned n, T alpha, const numc::func<T>& f) {
  return numc::func<T>([n, alpha, f](T x) {
    T val = f(x);
    if (n == 0) return T(1.0);
    if (n == 1) return T(2.0) * alpha * val;

    T c0 = T(1.0), c1 = T(2.0) * alpha * val, cn = T(0.0);
    for (unsigned i = 2; i <= n; ++i) {
      cn = (T(2.0) * val * (T(i) + alpha - T(1.0)) * c1 - (T(i) + T(2.0) * alpha - T(2.0)) * c0) / T(i);
      c0 = c1;
      c1 = cn;
    }
    return cn;
  });
}

/// @brief Computes the Jacobi polynomial, $P_n^{(\alpha, \beta)}(x)$.
/// Implemented using iterative recurrence relations for computational efficiency.
/// @tparam T Floating-point type.
/// @param n The degree of the polynomial.
/// @param alpha The alpha parameter ($\alpha > -1$).
/// @param beta The beta parameter ($\beta > -1$).
/// @param f The input expression.
/// @return An expression representing $P_n^{(\alpha, \beta)}(f(x))$.
template <typename T>
numc::func<T> jacobi(unsigned n, T alpha, T beta, const numc::func<T>& f) {
  return numc::func<T>([n, alpha, beta, f](T x) {
    T val = f(x);
    if (n == 0) return T(1.0);

    T p0 = T(1.0);
    T p1 = T(0.5) * (alpha - beta + (alpha + beta + T(2.0)) * val);
    if (n == 1) return p1;

    T pn = T(0.0);
    for (unsigned i = 2; i <= n; ++i) {
      T Ti = T(i);
      T c1 = T(2.0) * Ti * (Ti + alpha + beta) * (T(2.0) * Ti + alpha + beta - T(2.0));
      T c2 = (T(2.0) * Ti + alpha + beta - T(1.0)) *
             ((T(2.0) * Ti + alpha + beta) * (T(2.0) * Ti + alpha + beta - T(2.0)) * val + alpha * alpha - beta * beta);
      T c3 = T(2.0) * (Ti + alpha - T(1.0)) * (Ti + beta - T(1.0)) * (T(2.0) * Ti + alpha + beta);

      pn = (c2 * p1 - c3 * p0) / c1;
      p0 = p1;
      p1 = pn;
    }
    return pn;
  });
}

/// @} // End of orthogonal group

}  // namespace numc