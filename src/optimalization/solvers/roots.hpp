#pragma once

#include "../../common/function.hpp"
#include "../../common/complex.hpp"
#include "../../linear_algebra/solvers.hpp"
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

  Log::Warn("Bisection method reached the maximum number of iterations.");
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

  Log::Warn("Regula Falsi method reached the maximum number of iterations.");
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

  Log::Warn("Fixed-point iteration reached the maximum number of iterations.");
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

  Log::Warn("Newton's method reached the maximum number of iterations.");
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

  Log::Warn("Halley's method reached the maximum number of iterations.");
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

  Log::Warn("Secant method reached the maximum number of iterations.");
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

  Log::Warn("Inverse Quadratic Interpolation reached the maximum number of iterations.");
  return x3;
}

/// @brief Searches for a sign change in the function f(x) over the interval [a, b] with step dx.
template <typename T = double>
std::pair<T, T> root_search(const numc::func<T>& f, T a, T b, T dx) {
  T x1 = a;
  T f1 = f(x1);
  T x2 = a + dx;
  T f2 = f(x2);
  while (f1 * f2 > T(0.0)) {
    if (x1 >= b) {
      return {std::numeric_limits<T>::quiet_NaN(), std::numeric_limits<T>::quiet_NaN()};
    }
    x1 = x2;
    f1 = f2;
    x2 = x1 + dx;
    f2 = f(x2);
  }
  return {x1, x2};
}

/// @brief Bisection method that optionally detects singularities (switch=1).
template <typename T = double>
T bisect(const numc::func<T>& f, T x1, T x2, int switch_val = 0, T epsilon = T(1e-9)) {
  T f1 = f(x1);
  if (f1 == T(0.0)) return x1;
  T f2 = f(x2);
  if (f2 == T(0.0)) return x2;
  if (f1 * f2 > T(0.0)) {
    Log::Error("Root is not bracketed in bisect.");
    throw std::runtime_error("Root is not bracketed.");
  }
  int n = static_cast<int>(std::ceil(std::log(std::abs(x2 - x1) / epsilon) / std::log(T(2.0))));
  for (int i = 0; i < n; ++i) {
    T x3 = T(0.5) * (x1 + x2);
    T f3 = f(x3);
    if (switch_val == 1 && std::abs(f3) > std::abs(f1) && std::abs(f3) > std::abs(f2)) {
      return std::numeric_limits<T>::quiet_NaN();
    }
    if (f3 == T(0.0)) return x3;
    if (f2 * f3 < T(0.0)) {
      x1 = x3;
      f1 = f3;
    } else {
      x2 = x3;
      f2 = f3;
    }
  }
  return (x1 + x2) / T(2.0);
}

/// @brief Finds a root of f(x) using Brent's method.
template <typename T = double>
T brent(const numc::func<T>& f, T a, T b, T tol = T(1e-9)) {
  T x1 = a;
  T x2 = b;
  T f1 = f(x1);
  if (f1 == T(0.0)) return x1;
  T f2 = f(x2);
  if (f2 == T(0.0)) return x2;
  if (f1 * f2 > T(0.0)) {
    Log::Error("Root is not bracketed in brent.");
    throw std::runtime_error("Root is not bracketed.");
  }
  T x3 = T(0.5) * (a + b);
  for (int i = 0; i < 30; ++i) {
    T f3 = f(x3);
    if (std::abs(f3) < tol) return x3;
    if (f1 * f3 < T(0.0)) {
      b = x3;
    } else {
      a = x3;
    }
    if (b - a < tol * std::max(std::abs(b), T(1.0))) {
      return T(0.5) * (a + b);
    }
    T denom = (f2 - f1) * (f3 - f1) * (f2 - f3);
    T numer = x3 * (f1 - f2) * (f2 - f3 + f1) + f2 * x1 * (f2 - f3) + f1 * x2 * (f3 - f1);
    T dx;
    if (std::abs(denom) < T(1e-20)) {
      dx = b - a;
    } else {
      dx = f3 * numer / denom;
    }
    T x = x3 + dx;
    if ((b - x) * (x - a) < T(0.0)) {
      dx = T(0.5) * (b - a);
      x = a + dx;
    }
    if (x < x3) {
      x2 = x3;
      f2 = f3;
    } else {
      x1 = x3;
      f1 = f3;
    }
    x3 = x;
  }
  Log::Warn("Too many iterations in brent.");
  return x3;
}

/// @brief Safe Newton-Raphson method combined with bisection.
template <typename T = double>
T newton_raphson(const numc::func<T>& f, const numc::func<T>& df, T a, T b, T tol = T(1e-9)) {
  T fa = f(a);
  if (fa == T(0.0)) return a;
  T fb = f(b);
  if (fb == T(0.0)) return b;
  if (fa * fb > T(0.0)) {
    Log::Error("Root is not bracketed in newton_raphson.");
    throw std::runtime_error("Root is not bracketed.");
  }
  T x = T(0.5) * (a + b);
  for (int i = 0; i < 30; ++i) {
    T fx = f(x);
    if (std::abs(fx) < tol) return x;
    if (fa * fx < T(0.0)) {
      b = x;
    } else {
      a = x;
      fa = fx;
    }
    T dfx = df(x);
    T dx;
    if (std::abs(dfx) < T(1e-20)) {
      dx = b - a;
    } else {
      dx = -fx / dfx;
    }
    x += dx;
    if ((b - x) * (x - a) < T(0.0)) {
      dx = T(0.5) * (b - a);
      x = a + dx;
    }
    if (std::abs(dx) < tol * std::max(std::abs(b), T(1.0))) return x;
  }
  Log::Warn("Too many iterations in newton_raphson.");
  return x;
}

/// @brief Multidimensional Newton-Raphson for systems of equations.
template <typename T = double>
vector<T> newton_raphson2(std::function<vector<T>(const vector<T>&)> f, vector<T> x, T tol = T(1e-9)) {
  size_t n = x.size();
  auto jacobian = [&](vector<T>& xv) -> std::pair<tensor<T>, vector<T>> {
    T h = T(1e-4);
    tensor<T> jac({n, n});
    vector<T> f0 = f(xv);
    for (size_t j = 0; j < n; ++j) {
      T temp = xv[j];
      xv[j] = temp + h;
      vector<T> f1 = f(xv);
      xv[j] = temp;
      for (size_t i = 0; i < n; ++i) {
        jac(i, j) = (f1[i] - f0[i]) / h;
      }
    }
    return {jac, f0};
  };

  for (int iter = 0; iter < 30; ++iter) {
    auto [jac, f0] = jacobian(x);

    T f0_norm = f0.norm();
    if (f0_norm / std::sqrt(T(n)) < tol) {
      return x;
    }

    vector<T> minus_f0 = -f0;
    vector<T> dx = linear_algebra::gauss_pivot(jac, minus_f0);
    x = x + dx;

    T dx_norm = dx.norm();
    T x_norm = x.norm();
    if (dx_norm < tol * std::max(x_norm, T(1.0))) {
      return x;
    }
  }
  Log::Warn("Too many iterations in newton_raphson2.");
  return x;
}

/// @brief Evaluates a polynomial and its first and second derivatives at x.
template <typename T = double>
std::tuple<T, T, T> eval_poly(const vector<T>& a, T x) {
  size_t n = a.size() - 1;
  T p = a[n];
  T dp = T(0.0);
  T ddp = T(0.0);
  for (size_t i = 1; i <= n; ++i) {
    ddp = ddp * x + T(2.0) * dp;
    dp = dp * x + p;
    p = p * x + a[n - i];
  }
  return {p, dp, ddp};
}

/// @brief Evaluates a complex polynomial and its first and second derivatives at complex x.
template <typename T = double>
std::tuple<complex<T>, complex<T>, complex<T>> eval_poly(const std::vector<complex<T>>& a, complex<T> x) {
  size_t n = a.size() - 1;
  complex<T> p = a[n];
  complex<T> dp(0.0, 0.0);
  complex<T> ddp(0.0, 0.0);
  for (size_t i = 1; i <= n; ++i) {
    ddp = ddp * x + T(2.0) * dp;
    dp = dp * x + p;
    p = p * x + a[n - i];
  }
  return {p, dp, ddp};
}

/// @brief Helper for Laguerre's root finding algorithm on complex polynomials.
template <typename T = double>
complex<T> laguerre(const std::vector<complex<T>>& a, complex<T> x, T tol) {
  size_t n = a.size() - 1;
  T n_val = T(n);
  for (int iter = 0; iter < 30; ++iter) {
    auto [p, dp, ddp] = eval_poly(a, x);
    if (p.abs() < tol) return x;
    complex<T> g = dp / p;
    complex<T> h = g * g - ddp / p;

    complex<T> inside = (n_val - T(1.0)) * (n_val * h - g * g);
    T r_ins = inside.abs();
    T theta_ins = inside.arg();
    T root_r = std::sqrt(r_ins);
    T root_theta = theta_ins / T(2.0);
    complex<T> f(root_r * std::cos(root_theta), root_r * std::sin(root_theta));

    complex<T> denom1 = g + f;
    complex<T> denom2 = g - f;
    complex<T> dx;
    if (denom1.abs() > denom2.abs()) {
      dx = n_val / denom1;
    } else {
      dx = n_val / denom2;
    }
    x = x - dx;
    if (dx.abs() < tol * std::max(x.abs(), T(1.0))) return x;
  }
  Log::Warn("Too many iterations in Laguerre.");
  return x;
}

/// @brief Deflates polynomial a by removing a linear factor (x - root).
template <typename T = double>
std::vector<complex<T>> defl_poly(const std::vector<complex<T>>& a, complex<T> root) {
  size_t n = a.size() - 1;
  std::vector<complex<T>> b(n);
  b[n - 1] = a[n];
  for (int i = static_cast<int>(n) - 2; i >= 0; --i) {
    b[i] = a[i + 1] + root * b[i + 1];
  }
  return b;
}

/// @brief Computes all roots of a real polynomial using Laguerre's method.
template <typename T = double>
std::vector<complex<T>> poly_roots(const vector<T>& a, T tol = T(1e-12)) {
  size_t n = a.size() - 1;
  std::vector<complex<T>> a_comp(n + 1);
  for (size_t i = 0; i <= n; ++i) {
    a_comp[i] = complex<T>(a[i], T(0.0));
  }

  std::vector<complex<T>> roots_list(n);
  std::vector<complex<T>> current_a = a_comp;

  std::srand(12345);

  for (size_t i = 0; i < n; ++i) {
    T r_val = T(std::rand()) / T(RAND_MAX);
    complex<T> x(r_val, T(0.0));
    x = laguerre(current_a, x, tol);

    if (std::abs(x.imag()) < tol) {
      x = complex<T>(x.real(), T(0.0));
    }
    roots_list[i] = x;
    current_a = defl_poly(current_a, x);
  }
  return roots_list;
}

/// @brief Computes all roots of a complex polynomial using Laguerre's method.
template <typename T = double>
std::vector<complex<T>> poly_roots(const std::vector<complex<T>>& a, T tol = T(1e-12)) {
  size_t n = a.size() - 1;
  std::vector<complex<T>> roots_list(n);
  std::vector<complex<T>> current_a = a;

  std::srand(12345);

  for (size_t i = 0; i < n; ++i) {
    T r_val = T(std::rand()) / T(RAND_MAX);
    complex<T> x(r_val, T(0.0));
    x = laguerre(current_a, x, tol);

    if (std::abs(x.imag()) < tol) {
      x = complex<T>(x.real(), T(0.0));
    }
    roots_list[i] = x;
    current_a = defl_poly(current_a, x);
  }
  return roots_list;
}

}  // namespace optimalization

/// @} // End of optimize group

}  // namespace numc