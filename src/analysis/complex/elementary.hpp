#pragma once

#include "../../common/complex.hpp"
#include "../../common/function.hpp"
#include "../../common/vector.hpp"
#include "../../inc.hpp"
#include "../../utility/log.hpp"

namespace numc {
namespace analysis {

/// @addtogroup analysis
/// @{

/// @brief Computes the numerical derivative of a complex-valued function f(z) at z0.
template <typename T = double>
complex<T> complex_derivative(
    std::function<complex<T>(complex<T>)> f,
    complex<T> z0,
    T h = T(1e-6)) {
  complex<T> dz(h, T(0.0));
  complex<T> fz_plus = f(z0 + dz);
  complex<T> fz_minus = f(z0 - dz);
  return (fz_plus - fz_minus) / complex<T>(T(2.0) * h, T(0.0));
}

/// @brief Computes a contour integral of f(z) along a circle of radius R centered at z0.
/// Uses the trapezoidal rule on the parametric representation z = z0 + R*e^(it).
template <typename T = double>
complex<T> contour_integral_circle(
    std::function<complex<T>(complex<T>)> f,
    complex<T> z0,
    T R,
    size_t n = 1000) {
  T dt = T(2.0) * T(numc::PI) / T(n);
  complex<T> sum(T(0.0), T(0.0));
  for (size_t k = 0; k < n; ++k) {
    T t = T(k) * dt;
    complex<T> z = z0 + complex<T>(R * std::cos(t), R * std::sin(t));
    complex<T> dz(-R * std::sin(t) * dt, R * std::cos(t) * dt);
    sum = sum + f(z) * dz;
  }
  return sum;
}

/// @brief Computes the residue of f(z) at z0 using the contour integral formula.
/// Res(f, z0) = (1/2*pi*i) * integral_C f(z) dz
template <typename T = double>
complex<T> residue(
    std::function<complex<T>(complex<T>)> f,
    complex<T> z0,
    T R = T(0.01),
    size_t n = 10000) {
  complex<T> integral = contour_integral_circle(f, z0, R, n);
  T two_pi = T(2.0) * T(numc::PI);
  // Divide by 2*pi*i
  // (a + bi) / (0 + 2*pi*i) = (a + bi) * (0 - 2*pi*i) / (4*pi^2)
  // = (b*2*pi + (a * (-2*pi))i) / (4*pi^2)
  // = (b/(2*pi)) + (-a/(2*pi))i
  return complex<T>(integral.imag() / two_pi, -integral.real() / two_pi);
}

/// @brief Counts the number of zeros of f(z) inside a circle using the argument principle.
/// N = (1/2*pi*i) * integral_C f'(z)/f(z) dz
template <typename T = double>
int count_zeros_in_circle(
    std::function<complex<T>(complex<T>)> f,
    complex<T> center,
    T R,
    size_t n = 10000) {
  T h = T(1e-6);
  auto f_prime_over_f = [&](complex<T> z) -> complex<T> {
    complex<T> fz = f(z);
    complex<T> fz_plus = f(z + complex<T>(h, T(0.0)));
    complex<T> dfz = (fz_plus - fz) / complex<T>(h, T(0.0));
    return dfz / fz;
  };
  complex<T> integral = contour_integral_circle(f_prime_over_f, center, R, n);
  T two_pi = T(2.0) * T(numc::PI);
  // N = integral / (2*pi*i) => real part of integral.imag / two_pi
  return static_cast<int>(std::round(integral.imag() / two_pi));
}

/// @brief Laurent series coefficients (a_{-N} ... a_{N}) of f(z) around z0 on a circle of radius R.
/// a_n = (1/2*pi*i) * integral_C f(z) / (z - z0)^{n+1} dz
template <typename T = double>
std::vector<complex<T>> laurent_coeffs(
    std::function<complex<T>(complex<T>)> f,
    complex<T> z0,
    T R,
    int N,
    size_t n_quad = 2000) {
  std::vector<complex<T>> coeffs(2 * N + 1);
  T dt = T(2.0) * T(numc::PI) / T(n_quad);
  T two_pi = T(2.0) * T(numc::PI);

  for (int m = -N; m <= N; ++m) {
    complex<T> sum(T(0.0), T(0.0));
    for (size_t k = 0; k < n_quad; ++k) {
      T t = T(k) * dt;
      complex<T> z = z0 + complex<T>(R * std::cos(t), R * std::sin(t));
      complex<T> dz(-R * std::sin(t) * dt, R * std::cos(t) * dt);
      complex<T> w = z - z0;
      // f(z) / (z - z0)^{m+1}
      complex<T> w_pow = w;
      for (int p = 1; p < m + 1; ++p) w_pow = w_pow * w;
      if (m + 1 == 0) {
        w_pow = complex<T>(T(1.0), T(0.0));
      } else if (m + 1 < 0) {
        w_pow = complex<T>(T(1.0), T(0.0));
        for (int p = 0; p < -(m + 1); ++p) w_pow = w_pow * w;
        // invert: 1/w_pow
        T mag2 = w_pow.real() * w_pow.real() + w_pow.imag() * w_pow.imag();
        w_pow = complex<T>(w_pow.real() / mag2, -w_pow.imag() / mag2);
      }
      // Divide f(z) by w_pow
      complex<T> fz = f(z);
      T denom_mag2 = w_pow.real() * w_pow.real() + w_pow.imag() * w_pow.imag();
      if (denom_mag2 > T(1e-30)) {
        complex<T> ratio(
            (fz.real() * w_pow.real() + fz.imag() * w_pow.imag()) / denom_mag2,
            (fz.imag() * w_pow.real() - fz.real() * w_pow.imag()) / denom_mag2);
        sum = sum + ratio * dz;
      }
    }
    // Divide by 2*pi*i
    coeffs[m + N] = complex<T>(sum.imag() / two_pi, -sum.real() / two_pi);
  }
  return coeffs;
}

/// @brief Möbius transformation: w = (a*z + b) / (c*z + d).
template <typename T = double>
std::function<complex<T>(complex<T>)> mobius_transform(
    complex<T> a, complex<T> b,
    complex<T> c, complex<T> d) {
  return [a, b, c, d](complex<T> z) -> complex<T> {
    complex<T> num = a * z + b;
    complex<T> den = c * z + d;
    T mag2 = den.real() * den.real() + den.imag() * den.imag();
    return complex<T>(
        (num.real() * den.real() + num.imag() * den.imag()) / mag2,
        (num.imag() * den.real() - num.real() * den.imag()) / mag2);
  };
}

/// @}

}  // namespace analysis
}  // namespace numc
