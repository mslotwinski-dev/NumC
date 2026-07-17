#pragma once

#include "../common/vector.hpp"
#include "../common/tensor.hpp"
#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {
namespace differential_equations {

/// @addtogroup differential_equations
/// @{

/// @brief Solves the 1D heat equation u_t = alpha * u_xx using explicit (FTCS) method.
/// @return Grid of u values at each time step: result[t][x].
template <typename T = double>
std::vector<vector<T>> heat_explicit_1d(
    const vector<T>& u0,
    T alpha,
    T dx,
    T dt,
    size_t n_steps,
    T bc_left = T(0.0),
    T bc_right = T(0.0)) {
  size_t nx = u0.size();
  T r = alpha * dt / (dx * dx);
  if (r > T(0.5)) {
    Log::Warn("FTCS stability condition violated (r > 0.5). Results may be unstable.");
  }

  std::vector<vector<T>> result;
  vector<T> u = u0;
  result.push_back(u);

  for (size_t t = 0; t < n_steps; ++t) {
    vector<T> u_new(nx);
    u_new[0] = bc_left;
    u_new[nx - 1] = bc_right;
    for (size_t i = 1; i < nx - 1; ++i) {
      u_new[i] = u[i] + r * (u[i + 1] - T(2.0) * u[i] + u[i - 1]);
    }
    u = u_new;
    result.push_back(u);
  }
  return result;
}

/// @brief Solves the 1D heat equation using implicit (BTCS) method with Thomas algorithm.
template <typename T = double>
std::vector<vector<T>> heat_implicit_1d(
    const vector<T>& u0,
    T alpha,
    T dx,
    T dt,
    size_t n_steps,
    T bc_left = T(0.0),
    T bc_right = T(0.0)) {
  size_t nx = u0.size();
  T r = alpha * dt / (dx * dx);
  size_t n = nx - 2;  // interior points

  std::vector<vector<T>> result;
  vector<T> u = u0;
  result.push_back(u);

  for (size_t t = 0; t < n_steps; ++t) {
    // Setup tridiagonal system
    vector<T> a(n), b(n), c(n), d(n);
    for (size_t i = 0; i < n; ++i) {
      a[i] = -r;
      b[i] = T(1.0) + T(2.0) * r;
      c[i] = -r;
      d[i] = u[i + 1];
    }
    d[0] += r * bc_left;
    d[n - 1] += r * bc_right;

    // Thomas algorithm
    for (size_t i = 1; i < n; ++i) {
      T m = a[i] / b[i - 1];
      b[i] -= m * c[i - 1];
      d[i] -= m * d[i - 1];
    }
    vector<T> x_sol(n);
    x_sol[n - 1] = d[n - 1] / b[n - 1];
    for (int i = static_cast<int>(n) - 2; i >= 0; --i) {
      x_sol[i] = (d[i] - c[i] * x_sol[i + 1]) / b[i];
    }

    vector<T> u_new(nx);
    u_new[0] = bc_left;
    for (size_t i = 0; i < n; ++i) u_new[i + 1] = x_sol[i];
    u_new[nx - 1] = bc_right;
    u = u_new;
    result.push_back(u);
  }
  return result;
}

/// @brief Solves the 1D wave equation u_tt = c^2 * u_xx using explicit method.
template <typename T = double>
std::vector<vector<T>> wave_explicit_1d(
    const vector<T>& u0,
    const vector<T>& v0,
    T c,
    T dx,
    T dt,
    size_t n_steps,
    T bc_left = T(0.0),
    T bc_right = T(0.0)) {
  size_t nx = u0.size();
  T r = c * dt / dx;
  T r2 = r * r;

  std::vector<vector<T>> result;
  result.push_back(u0);

  // First step using initial velocity
  vector<T> u_prev = u0;
  vector<T> u_curr(nx);
  u_curr[0] = bc_left;
  u_curr[nx - 1] = bc_right;
  for (size_t i = 1; i < nx - 1; ++i) {
    u_curr[i] = u0[i] + dt * v0[i] +
                T(0.5) * r2 * (u0[i + 1] - T(2.0) * u0[i] + u0[i - 1]);
  }
  result.push_back(u_curr);

  for (size_t t = 1; t < n_steps; ++t) {
    vector<T> u_next(nx);
    u_next[0] = bc_left;
    u_next[nx - 1] = bc_right;
    for (size_t i = 1; i < nx - 1; ++i) {
      u_next[i] = T(2.0) * u_curr[i] - u_prev[i] +
                  r2 * (u_curr[i + 1] - T(2.0) * u_curr[i] + u_curr[i - 1]);
    }
    u_prev = u_curr;
    u_curr = u_next;
    result.push_back(u_curr);
  }
  return result;
}

/// @brief Solves the 2D Laplace equation on a rectangular grid using Gauss-Seidel iteration.
/// u_xx + u_yy = 0 with Dirichlet boundary conditions.
template <typename T = double>
tensor<T> laplace_2d(
    size_t nx, size_t ny,
    std::function<T(size_t, size_t)> bc,
    T tol = T(1e-6),
    int max_iter = 10000) {
  tensor<T> u({nx, ny});

  // Apply boundary conditions
  for (size_t i = 0; i < nx; ++i) {
    u(i, 0) = bc(i, 0);
    u(i, ny - 1) = bc(i, ny - 1);
  }
  for (size_t j = 0; j < ny; ++j) {
    u(0, j) = bc(0, j);
    u(nx - 1, j) = bc(nx - 1, j);
  }

  for (int iter = 0; iter < max_iter; ++iter) {
    T max_diff = T(0.0);
    for (size_t i = 1; i < nx - 1; ++i) {
      for (size_t j = 1; j < ny - 1; ++j) {
        T u_old = u(i, j);
        u(i, j) = T(0.25) * (u(i + 1, j) + u(i - 1, j) + u(i, j + 1) + u(i, j - 1));
        max_diff = std::max(max_diff, std::abs(u(i, j) - u_old));
      }
    }
    if (max_diff < tol) break;
  }
  return u;
}

/// @brief Solves the 2D Poisson equation u_xx + u_yy = f(x, y) using SOR iteration.
template <typename T = double>
tensor<T> poisson_2d(
    size_t nx, size_t ny,
    T dx, T dy,
    std::function<T(T, T)> source,
    std::function<T(size_t, size_t)> bc,
    T omega = T(1.5),
    T tol = T(1e-6),
    int max_iter = 10000) {
  tensor<T> u({nx, ny});

  // Apply boundary conditions
  for (size_t i = 0; i < nx; ++i) {
    u(i, 0) = bc(i, 0);
    u(i, ny - 1) = bc(i, ny - 1);
  }
  for (size_t j = 0; j < ny; ++j) {
    u(0, j) = bc(0, j);
    u(nx - 1, j) = bc(nx - 1, j);
  }

  T dx2 = dx * dx, dy2 = dy * dy;
  T denom = T(2.0) * (dx2 + dy2);

  for (int iter = 0; iter < max_iter; ++iter) {
    T max_diff = T(0.0);
    for (size_t i = 1; i < nx - 1; ++i) {
      for (size_t j = 1; j < ny - 1; ++j) {
        T u_old = u(i, j);
        T x = T(i) * dx, y = T(j) * dy;
        T u_gs = (dy2 * (u(i + 1, j) + u(i - 1, j)) +
                  dx2 * (u(i, j + 1) + u(i, j - 1)) -
                  dx2 * dy2 * source(x, y)) / denom;
        u(i, j) = (T(1.0) - omega) * u_old + omega * u_gs;
        max_diff = std::max(max_diff, std::abs(u(i, j) - u_old));
      }
    }
    if (max_diff < tol) break;
  }
  return u;
}

/// @}

}  // namespace differential_equations
}  // namespace numc
