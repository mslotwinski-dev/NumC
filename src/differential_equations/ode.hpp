#pragma once

#include "../common/vector.hpp"
#include "../common/tensor.hpp"
#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {
namespace differential_equations {

/// @addtogroup differential_equations
/// @{

/// @brief Solves initial value problem using 4th-order Taylor series method.
template <typename T = double>
std::pair<std::vector<T>, std::vector<vector<T>>> taylor(
    std::function<tensor<T>(T, const vector<T>&)> deriv,
    T x,
    vector<T> y,
    T x_stop,
    T h) {
  std::vector<T> X;
  std::vector<vector<T>> Y;
  X.push_back(x);
  Y.push_back(y);
  size_t n = y.size();

  while (x < x_stop) {
    h = std::min(h, x_stop - x);
    tensor<T> D = deriv(x, y);
    T H = T(1.0);
    for (size_t j = 0; j < 4; ++j) {
      H = H * h / T(j + 1);
      for (size_t i = 0; i < n; ++i) {
        y[i] += D(j, i) * H;
      }
    }
    x += h;
    X.push_back(x);
    Y.push_back(y);
  }
  return {X, Y};
}

/// @brief Solves initial value problem using 4th-order Runge-Kutta method.
template <typename T = double>
std::pair<std::vector<T>, std::vector<vector<T>>> run_kut4(
    std::function<vector<T>(T, const vector<T>&)> F,
    T x,
    vector<T> y,
    T x_stop,
    T h) {
  auto step = [&](T xv, const vector<T>& yv, T hv) -> vector<T> {
    vector<T> K0 = hv * F(xv, yv);
    vector<T> K1 = hv * F(xv + hv / T(2.0), yv + K0 / T(2.0));
    vector<T> K2 = hv * F(xv + hv / T(2.0), yv + K1 / T(2.0));
    vector<T> K3 = hv * F(xv + hv, yv + K2);
    return (K0 + T(2.0) * K1 + T(2.0) * K2 + K3) / T(6.0);
  };

  std::vector<T> X;
  std::vector<vector<T>> Y;
  X.push_back(x);
  Y.push_back(y);

  while (x < x_stop) {
    h = std::min(h, x_stop - x);
    y = y + step(x, y, h);
    x += h;
    X.push_back(x);
    Y.push_back(y);
  }
  return {X, Y};
}

/// @brief Solves initial value problem using adaptive Runge-Kutta method (Cash-Karp coefficients).
template <typename T = double>
std::pair<std::vector<T>, std::vector<vector<T>>> run_kut5(
    std::function<vector<T>(T, const vector<T>&)> F,
    T x,
    vector<T> y,
    T x_stop,
    T h,
    T tol = T(1e-6)) {
  size_t n = y.size();
  T C[6] = {T(37.0) / T(378.0), T(0.0), T(250.0) / T(621.0), T(125.0) / T(594.0), T(0.0), T(512.0) / T(1771.0)};
  T D[6] = {T(2825.0) / T(27648.0), T(0.0), T(18575.0) / T(48384.0), T(13525.0) / T(55296.0), T(277.0) / T(14336.0), T(0.25)};

  auto step = [&](T xv, const vector<T>& yv, T hv, vector<T>& dy, T& e) {
    std::vector<vector<T>> K(6, vector<T>(n));
    K[0] = hv * F(xv, yv);
    K[1] = hv * F(xv + T(0.2) * hv, yv + T(0.2) * K[0]);
    K[2] = hv * F(xv + T(0.3) * hv, yv + T(0.075) * K[0] + T(0.225) * K[1]);
    K[3] = hv * F(xv + T(0.6) * hv, yv + T(0.3) * K[0] - T(0.9) * K[1] + T(1.2) * K[2]);
    K[4] = hv * F(xv + hv, yv - (T(11.0) / T(54.0)) * K[0] + T(2.5) * K[1] - (T(70.0) / T(27.0)) * K[2] + (T(35.0) / T(27.0)) * K[3]);
    K[5] = hv * F(xv + T(0.875) * hv, yv + (T(1631.0) / T(55296.0)) * K[0] + (T(175.0) / T(512.0)) * K[1] + (T(575.0) / T(13824.0)) * K[2] +
                              (T(44275.0) / T(110592.0)) * K[3] + (T(253.0) / T(4096.0)) * K[4]);

    dy = vector<T>(n);
    vector<T> E_vec(n);
    for (int i = 0; i < 6; ++i) {
      dy = dy + C[i] * K[i];
      E_vec = E_vec + (C[i] - D[i]) * K[i];
    }
    T sum_sq = T(0.0);
    for (size_t i = 0; i < n; ++i) {
      sum_sq += E_vec[i] * E_vec[i];
    }
    e = std::sqrt(sum_sq / T(n));
  };

  std::vector<T> X;
  std::vector<vector<T>> Y;
  X.push_back(x);
  Y.push_back(y);

  int stopper = 0;
  for (int iter = 0; iter < 10000; ++iter) {
    vector<T> dy(n);
    T e = T(0.0);
    step(x, y, h, dy, e);

    if (e <= tol) {
      y = y + dy;
      x += h;
      X.push_back(x);
      Y.push_back(y);
    }
    if (stopper == 1) break;
    T hNext;
    if (e != T(0.0)) {
      hNext = T(0.9) * h * std::pow(tol / e, T(0.2));
    } else {
      hNext = h;
    }
    if (x + hNext >= x_stop) {
      hNext = x_stop - x;
      stopper = 1;
    }
    h = hNext;
  }
  return {X, Y};
}

/// @brief Advances y from x to x_stop in n_steps using the midpoint method.
template <typename T = double>
vector<T> midpoint_step(
    std::function<vector<T>(T, const vector<T>&)> F,
    T x,
    const vector<T>& y,
    T x_stop,
    int n_steps) {
  T h = (x_stop - x) / T(n_steps);
  vector<T> y0 = y;
  vector<T> y1 = y0 + h * F(x, y0);
  vector<T> y2 = y1;
  for (int i = 0; i < n_steps - 1; ++i) {
    x += h;
    y2 = y0 + T(2.0) * h * F(x, y1);
    y0 = y1;
    y1 = y2;
  }
  return T(0.5) * (y1 + y0 + h * F(x, y2));
}

/// @brief Inner integration loop of the Bulirsch-Stoer method.
template <typename T = double>
vector<T> midpoint_integrate(
    std::function<vector<T>(T, const vector<T>&)> F,
    T x,
    const vector<T>& y,
    T x_stop,
    T tol = T(1e-6)) {
  size_t n = y.size();
  std::vector<vector<T>> r(51, vector<T>(n));
  int n_steps = 2;
  r[1] = midpoint_step(F, x, y, x_stop, n_steps);
  vector<T> r_old = r[1];

  for (int k = 2; k < 51; ++k) {
    n_steps *= 2;
    r[k] = midpoint_step(F, x, y, x_stop, n_steps);
    // Richardson extrapolation
    for (int j = k - 1; j > 0; --j) {
      T const_factor = std::pow(T(4.0), k - j);
      r[j] = (const_factor * r[j + 1] - r[j]) / (const_factor - T(1.0));
    }
    T sum_sq = T(0.0);
    for (size_t i = 0; i < n; ++i) {
      T diff = r[1][i] - r_old[i];
      sum_sq += diff * diff;
    }
    T e = std::sqrt(sum_sq / T(n));
    if (e < tol) return r[1];
    r_old = r[1];
  }
  Log::Warn("Midpoint integration did not converge.");
  return r[1];
}

/// @brief Solves initial value problem using the Bulirsch-Stoer method.
template <typename T = double>
std::pair<std::vector<T>, std::vector<vector<T>>> bul_stoer(
    std::function<vector<T>(T, const vector<T>&)> F,
    T x,
    vector<T> y,
    T x_stop,
    T H,
    T tol = T(1e-6)) {
  std::vector<T> X;
  std::vector<vector<T>> Y;
  X.push_back(x);
  Y.push_back(y);
  while (x < x_stop) {
    H = std::min(H, x_stop - x);
    y = midpoint_integrate(F, x, y, x + H, tol);
    x += H;
    X.push_back(x);
    Y.push_back(y);
  }
  return {X, Y};
}

/// @brief Finds the root of f(x) = 0 using linear interpolation.
template <typename T = double>
T lin_interp(std::function<T(T)> f, T x1, T x2) {
  T f1 = f(x1);
  T f2 = f(x2);
  return x2 - f2 * (x2 - x1) / (f2 - f1);
}

/// @}

}  // namespace differential_equations
}  // namespace numc
