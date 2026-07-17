#pragma once

#include "../common/vector.hpp"
#include "../common/tensor.hpp"
#include "../linear_algebra/solvers.hpp"
#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {
namespace differential_equations {

/// @addtogroup differential_equations
/// @{

/// @brief Solves a two-point BVP using the finite difference method.
/// Solves y'' = f(x, y, y') with boundary conditions y(a) = alpha, y(b) = beta.
/// Uses the linearized form: y'' = p(x)*y' + q(x)*y + r(x).
template <typename T = double>
std::pair<vector<T>, vector<T>> finite_diff_bvp(
    std::function<T(T)> p,
    std::function<T(T)> q,
    std::function<T(T)> r,
    T a, T b,
    T alpha, T beta,
    size_t n) {
  T h = (b - a) / T(n);

  // Setup tridiagonal system
  vector<T> sub(n - 1);
  vector<T> diag(n - 1);
  vector<T> sup(n - 1);
  vector<T> rhs(n - 1);

  for (size_t i = 0; i < n - 1; ++i) {
    T x = a + T(i + 1) * h;
    T pi = p(x);
    T qi = q(x);
    T ri = r(x);

    sub[i] = T(1.0) - h * pi / T(2.0);
    diag[i] = -T(2.0) + h * h * qi;
    sup[i] = T(1.0) + h * pi / T(2.0);
    rhs[i] = h * h * ri;
  }

  // Apply boundary conditions
  rhs[0] -= sub[0] * alpha;
  rhs[n - 2] -= sup[n - 2] * beta;

  // Solve tridiagonal system (Thomas algorithm)
  // Forward sweep
  for (size_t i = 1; i < n - 1; ++i) {
    T m = sub[i] / diag[i - 1];
    diag[i] -= m * sup[i - 1];
    rhs[i] -= m * rhs[i - 1];
  }
  // Back substitution
  vector<T> y_inner(n - 1);
  y_inner[n - 2] = rhs[n - 2] / diag[n - 2];
  for (int i = static_cast<int>(n) - 3; i >= 0; --i) {
    y_inner[i] = (rhs[i] - sup[i] * y_inner[i + 1]) / diag[i];
  }

  // Build full solution
  vector<T> x_full(n + 1);
  vector<T> y_full(n + 1);
  for (size_t i = 0; i <= n; ++i) x_full[i] = a + T(i) * h;
  y_full[0] = alpha;
  for (size_t i = 0; i < n - 1; ++i) y_full[i + 1] = y_inner[i];
  y_full[n] = beta;

  return {x_full, y_full};
}

/// @brief Solves a two-point BVP using the shooting method combined with RK4.
template <typename T = double>
std::pair<vector<T>, vector<T>> shooting_method(
    std::function<vector<T>(T, const vector<T>&)> F,
    T a, T b,
    T alpha, T beta,
    T slope_guess1 = T(0.0),
    T slope_guess2 = T(1.0),
    T tol = T(1e-6),
    size_t n_steps = 100) {
  T h = (b - a) / T(n_steps);

  auto integrate_rk4 = [&](T s0) -> std::pair<std::vector<T>, std::vector<vector<T>>> {
    std::vector<T> X;
    std::vector<vector<T>> Y;
    T x = a;
    vector<T> y = {alpha, s0};
    X.push_back(x);
    Y.push_back(y);
    for (size_t i = 0; i < n_steps; ++i) {
      vector<T> K0 = h * F(x, y);
      vector<T> K1 = h * F(x + h / T(2.0), y + K0 / T(2.0));
      vector<T> K2 = h * F(x + h / T(2.0), y + K1 / T(2.0));
      vector<T> K3 = h * F(x + h, y + K2);
      y = y + (K0 + T(2.0) * K1 + T(2.0) * K2 + K3) / T(6.0);
      x += h;
      X.push_back(x);
      Y.push_back(y);
    }
    return {X, Y};
  };

  // Secant method on the residual r(s) = y(b; s) - beta
  T s1 = slope_guess1, s2 = slope_guess2;
  auto [X1, Y1] = integrate_rk4(s1);
  T r1 = Y1.back()[0] - beta;
  auto [X2, Y2] = integrate_rk4(s2);
  T r2 = Y2.back()[0] - beta;

  for (int iter = 0; iter < 50; ++iter) {
    if (std::abs(r2) < tol) {
      // Extract just x and y[0]
      vector<T> x_out(X2.size());
      vector<T> y_out(X2.size());
      for (size_t i = 0; i < X2.size(); ++i) {
        x_out[i] = X2[i];
        y_out[i] = Y2[i][0];
      }
      return {x_out, y_out};
    }
    if (std::abs(r2 - r1) < T(1e-15)) break;
    T s_new = s2 - r2 * (s2 - s1) / (r2 - r1);
    s1 = s2; r1 = r2;
    s2 = s_new;
    auto [Xn, Yn] = integrate_rk4(s2);
    X2 = Xn; Y2 = Yn;
    r2 = Y2.back()[0] - beta;
  }
  Log::Warn("Shooting method did not converge.");
  vector<T> x_out(X2.size());
  vector<T> y_out(X2.size());
  for (size_t i = 0; i < X2.size(); ++i) {
    x_out[i] = X2[i];
    y_out[i] = Y2[i][0];
  }
  return {x_out, y_out};
}

/// @}

}  // namespace differential_equations
}  // namespace numc
