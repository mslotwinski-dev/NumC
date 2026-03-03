#pragma once

#include "../../common/function.hpp"
#include "../../inc.hpp"
#include "../../utility/log.hpp"

namespace numc {

/// @defgroup optimize Optimization and Solvers
/// @brief Algorithms for finding roots (zeros) and local extrema of mathematical functions.
/// @{

namespace optimalization {

/// @brief Finds a root of the function f(x) using the Bisection method.
/// This method is robust and guaranteed to converge if the root is bracketed.
/// @tparam T Floating-point type.
/// @param f The function (numc::func object).
/// @param a The start of the interval.
/// @param b The end of the interval.
/// @param tol Error tolerance (default is 1e-6).
/// @param max_iter Maximum number of iterations.
/// @return The approximated root value.
template <typename T = double>
T bisection(const numc::func<T>& f, T a, T b, T tol = T(1e-6), int max_iter = 100) {
  T fa = f(a);
  T fb = f(b);

  if (fa * fb > T(0.0)) {
    Log::Error("The function must have opposite signs at the endpoints of the interval [a, b]!");
    return std::numeric_limits<T>::quiet_NaN();
  }

  T c = a;
  for (int i = 0; i < max_iter; ++i) {
    c = (a + b) / T(2.0);
    T fc = f(c);

    // Check if we hit the root exactly or are within tolerance
    if (std::abs(fc) < tol || (b - a) / T(2.0) < tol) {
      return c;
    }

    // Narrow down the interval
    if (fa * fc < T(0.0)) {
      b = c;
      fb = fc;
    } else {
      a = c;
      fa = fc;
    }
  }

  Log::Warning("Bisection method reached the maximum number of iterations.");
  return c;
}

/// @brief Finds a root of the function f(x) using the False Position method (Regula Falsi).
/// Combines the robustness of bisection with the linear interpolation of the secant method.
/// @tparam T Floating-point type.
/// @param f The function (numc::func object).
/// @param a The start of the interval.
/// @param b The end of the interval.
/// @param tol Error tolerance.
/// @param max_iter Maximum number of iterations.
/// @return The approximated root value.
template <typename T = double>
T regula_falsi(const numc::func<T>& f, T a, T b, T tol = T(1e-6), int max_iter = 100) {
  T fa = f(a);
  T fb = f(b);

  if (fa * fb > T(0.0)) {
    Log::Error("The function must have opposite signs at the endpoints of the interval [a, b]!");
    return std::numeric_limits<T>::quiet_NaN();
  }

  T c = a;
  for (int i = 0; i < max_iter; ++i) {
    if (std::abs(fa - fb) < T(1e-12)) {
      break;  // Prevent division by zero
    }

    // Secant line x-intercept
    c = (a * fb - b * fa) / (fb - fa);
    T fc = f(c);

    if (std::abs(fc) < tol) {
      return c;
    }

    // Narrow down the interval
    if (fa * fc < T(0.0)) {
      b = c;
      fb = fc;
    } else {
      a = c;
      fa = fc;
    }
  }

  Log::Warning("Regula Falsi method reached the maximum number of iterations.");
  return c;
}

/// @brief Finds the fixed point of the function, solving $x = f(x)$.
/// Uses the Fixed-point iteration method.
/// @tparam T Floating-point type.
/// @param f The function (numc::func object) defining the iteration $x_{n+1} = f(x_n)$.
/// @param x0 The initial guess.
/// @param tol Error tolerance.
/// @param max_iter Maximum number of iterations.
/// @return The approximated fixed point.
template <typename T = double>
T fixed_point(const numc::func<T>& f, T x0, T tol = T(1e-6), int max_iter = 100) {
  T x = x0;

  for (int i = 0; i < max_iter; ++i) {
    T x_next = f(x);

    if (std::abs(x_next - x) < tol) {
      return x_next;
    }
    x = x_next;
  }

  Log::Warning("Fixed-point iteration reached the maximum number of iterations.");
  return x;
}

/// @brief Finds a root of the function f(x) using the Newton-Raphson method.
/// Automatically approximates the first derivative using the central difference quotient.
/// @tparam T Floating-point type.
/// @param f The function (numc::func object).
/// @param x0 The initial guess.
/// @param tol Error tolerance.
/// @param max_iter Maximum number of iterations.
/// @return The approximated root value.
template <typename T = double>
T newton(const numc::func<T>& f, T x0, T tol = T(1e-6), int max_iter = 100) {
  T x = x0;
  T h = T(1e-5);  // Small step for numerical derivative

  for (int i = 0; i < max_iter; ++i) {
    T fx = f(x);

    if (std::abs(fx) < tol) {
      return x;
    }

    // Numerical derivative using central difference: f'(x) ≈ (f(x+h) - f(x-h)) / 2h
    T dfx = (f(x + h) - f(x - h)) / (T(2.0) * h);

    if (std::abs(dfx) < T(1e-12)) {
      Log::Error("Derivative is near zero in Newton's method. Algorithm stalled.");
      return std::numeric_limits<T>::quiet_NaN();
    }

    x = x - (fx / dfx);
  }

  Log::Warning("Newton's method reached the maximum number of iterations.");
  return x;
}

/// @brief Finds a root of the function f(x) using Halley's method.
/// A cubic convergence algorithm that uses numerical approximations for both first and second derivatives.
/// @tparam T Floating-point type.
/// @param f The function (numc::func object).
/// @param x0 The initial guess.
/// @param tol Error tolerance.
/// @param max_iter Maximum number of iterations.
/// @return The approximated root value.
template <typename T = double>
T halley(const numc::func<T>& f, T x0, T tol = T(1e-6), int max_iter = 100) {
  T x = x0;
  T h = T(1e-4);  // Slightly larger step to maintain numerical stability for second derivative

  for (int i = 0; i < max_iter; ++i) {
    T fx = f(x);

    if (std::abs(fx) < tol) {
      return x;
    }

    // First derivative approximation
    T dfx = (f(x + h) - f(x - h)) / (T(2.0) * h);

    // Second derivative approximation
    T ddfx = (f(x + h) - T(2.0) * fx + f(x - h)) / (h * h);

    T denominator = T(2.0) * dfx * dfx - fx * ddfx;
    if (std::abs(denominator) < T(1e-12)) {
      Log::Error("Denominator is near zero in Halley's method. Algorithm stalled.");
      return std::numeric_limits<T>::quiet_NaN();
    }

    T x_next = x - (T(2.0) * fx * dfx) / denominator;

    if (std::abs(x_next - x) < tol) {
      return x_next;
    }
    x = x_next;
  }

  Log::Warning("Halley's method reached the maximum number of iterations.");
  return x;
}

/// @brief Finds a root of the function f(x) using the Secant method.
/// Replaces the derivative in Newton's method with a finite difference approximation based on two prior points.
/// @tparam T Floating-point type.
/// @param f The function (numc::func object).
/// @param x0 The first initial guess.
/// @param x1 The second initial guess.
/// @param tol Error tolerance.
/// @param max_iter Maximum number of iterations.
/// @return The approximated root value.
template <typename T = double>
T secant(const numc::func<T>& f, T x0, T x1, T tol = T(1e-6), int max_iter = 100) {
  T f0 = f(x0);
  T f1 = f(x1);

  for (int i = 0; i < max_iter; ++i) {
    if (std::abs(f1 - f0) < T(1e-12)) {
      Log::Error("Difference in function values is near zero. Secant method stalled.");
      return std::numeric_limits<T>::quiet_NaN();
    }

    T x2 = x1 - f1 * (x1 - x0) / (f1 - f0);

    if (std::abs(x2 - x1) < tol || std::abs(f(x2)) < tol) {
      return x2;
    }

    // Shift points for next iteration
    x0 = x1;
    f0 = f1;
    x1 = x2;
    f1 = f(x2);
  }

  Log::Warning("Secant method reached the maximum number of iterations.");
  return x1;
}

/// @brief Finds a root of the function f(x) using Inverse Quadratic Interpolation.
/// Fits a parabola passing through three points mapping y to x, yielding rapid convergence.
/// @tparam T Floating-point type.
/// @param f The function (numc::func object).
/// @param x0 The first initial guess.
/// @param x1 The second initial guess.
/// @param x2 The third initial guess.
/// @param tol Error tolerance.
/// @param max_iter Maximum number of iterations.
/// @return The approximated root value.
template <typename T = double>
T inverse_quadratic(const numc::func<T>& f, T x0, T x1, T x2, T tol = T(1e-6), int max_iter = 100) {
  T f0 = f(x0), f1 = f(x1), f2 = f(x2);
  T x3 = x2;

  for (int i = 0; i < max_iter; ++i) {
    // Check for potential division by zero
    if (std::abs(f0 - f1) < T(1e-12) || std::abs(f1 - f2) < T(1e-12) || std::abs(f0 - f2) < T(1e-12)) {
      Log::Error("Denominator near zero in Inverse Quadratic Interpolation. Stall detected.");
      return std::numeric_limits<T>::quiet_NaN();
    }

    // Lagrange interpolation weights mapped inversely
    T L0 = (f1 * f2) / ((f0 - f1) * (f0 - f2));
    T L1 = (f0 * f2) / ((f1 - f0) * (f1 - f2));
    T L2 = (f0 * f1) / ((f2 - f0) * (f2 - f1));

    x3 = x0 * L0 + x1 * L1 + x2 * L2;
    T f3 = f(x3);

    if (std::abs(f3) < tol || std::abs(x3 - x2) < tol) {
      return x3;
    }

    // Shift points for the next iteration
    x0 = x1;
    f0 = f1;
    x1 = x2;
    f1 = f2;
    x2 = x3;
    f2 = f3;
  }

  Log::Warning("Inverse Quadratic Interpolation reached the maximum number of iterations.");
  return x3;
}

}  // namespace optimalization

/// @} // End of optimize group

}  // namespace numc