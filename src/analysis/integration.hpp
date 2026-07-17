#pragma once

#include "../common/vector.hpp"
#include "../common/function.hpp"
#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {
namespace analysis {

/// @addtogroup analysis
/// @{

/// @brief Recursive trapezoidal rule evaluation with 2^(k-1) panels.
template <typename T = double>
T trapezoid(numc::func<T> f, T a, T b, T old_integral, int k) {
  if (k == 1) {
    return (f(a) + f(b)) * (b - a) / T(2.0);
  } else {
    size_t n = 1ULL << (k - 2);
    T h = (b - a) / T(n);
    T x = a + h / T(2.0);
    T sum = T(0.0);
    for (size_t i = 0; i < n; ++i) {
      sum += f(x);
      x += h;
    }
    return (old_integral + h * sum) / T(2.0);
  }
}

/// @brief Performs Romberg integration of f(x) from a to b.
template <typename T = double>
std::pair<T, int> romberg(numc::func<T> f, T a, T b, T tol = T(1e-6)) {
  std::vector<T> r(21, T(0.0));
  r[1] = trapezoid(f, a, b, T(0.0), 1);
  T r_old = r[1];
  for (int k = 2; k < 21; ++k) {
    r[k] = trapezoid(f, a, b, r[k - 1], k);
    // Richardson extrapolation
    for (int j = k - 1; j > 0; --j) {
      T const_factor = std::pow(T(4.0), k - j);
      r[j] = (const_factor * r[j + 1] - r[j]) / (const_factor - T(1.0));
    }
    if (std::abs(r[1] - r_old) < tol * std::max(std::abs(r[1]), T(1.0))) {
      return {r[1], 1 << (k - 1)};
    }
    r_old = r[1];
  }
  Log::Warn("Romberg quadrature did not converge.");
  return {r[1], 1 << 19};
}

/// @brief Computes the integral of f(x) from a to b using Simpson's 1/3 rule.
template <typename T = double>
T simpson_1_3(numc::func<T> f, T a, T b, size_t n = 100) {
  if (n % 2 != 0) ++n; // n must be even
  T h = (b - a) / T(n);
  T sum = f(a) + f(b);
  for (size_t i = 1; i < n; i += 2) sum += T(4.0) * f(a + T(i) * h);
  for (size_t i = 2; i < n - 1; i += 2) sum += T(2.0) * f(a + T(i) * h);
  return sum * h / T(3.0);
}

/// @brief Computes the integral of f(x) from a to b using Simpson's 3/8 rule.
template <typename T = double>
T simpson_3_8(numc::func<T> f, T a, T b, size_t n = 99) {
  if (n % 3 != 0) n += 3 - (n % 3); // n must be multiple of 3
  T h = (b - a) / T(n);
  T sum = f(a) + f(b);
  for (size_t i = 1; i < n; ++i) {
    if (i % 3 == 0) sum += T(2.0) * f(a + T(i) * h);
    else sum += T(3.0) * f(a + T(i) * h);
  }
  return sum * T(3.0) * h / T(8.0);
}

/// @brief Computes the integral of f(x) from a to b using Boole's rule.
template <typename T = double>
T boole(numc::func<T> f, T a, T b, size_t n = 100) {
  if (n % 4 != 0) n += 4 - (n % 4); // n must be multiple of 4
  T h = (b - a) / T(n);
  T sum = T(7.0) * (f(a) + f(b));
  for (size_t i = 1; i < n; ++i) {
    if (i % 2 != 0) sum += T(32.0) * f(a + T(i) * h);
    else if (i % 4 == 2) sum += T(12.0) * f(a + T(i) * h);
    else if (i % 4 == 0) sum += T(14.0) * f(a + T(i) * h);
  }
  return sum * T(2.0) * h / T(45.0);
}

/// @brief Helper for adaptive simpson
template <typename T = double>
T adaptive_simpson_helper(numc::func<T> f, T a, T b, T eps, T S, T fa, T fc, T fb, int bottom) {
  T c = (a + b) / T(2.0);
  T h = b - a;
  T d = (a + c) / T(2.0);
  T e = (c + b) / T(2.0);
  T fd = f(d);
  T fe = f(e);
  T Sleft = (h / T(12.0)) * (fa + T(4.0) * fd + fc);
  T Sright = (h / T(12.0)) * (fc + T(4.0) * fe + fb);
  T S2 = Sleft + Sright;
  if (bottom <= 0 || std::abs(S2 - S) <= T(15.0) * eps) {
    return S2 + (S2 - S) / T(15.0);
  }
  return adaptive_simpson_helper(f, a, c, eps / T(2.0), Sleft, fa, fd, fc, bottom - 1) +
         adaptive_simpson_helper(f, c, b, eps / T(2.0), Sright, fc, fe, fb, bottom - 1);
}

/// @brief Computes the integral of f(x) from a to b using adaptive Simpson's quadrature.
template <typename T = double>
T adaptive_simpson(numc::func<T> f, T a, T b, T tol = T(1e-6), int max_depth = 15) {
  T c = (a + b) / T(2.0);
  T h = b - a;
  T fa = f(a);
  T fb = f(b);
  T fc = f(c);
  T S = (h / T(6.0)) * (fa + T(4.0) * fc + fb);
  return adaptive_simpson_helper(f, a, b, tol, S, fa, fc, fb, max_depth);
}

/// @brief Computes the 2D integral of f(x,y) over the rectangle [ax, bx] x [ay, by] using nested adaptive Simpson's quadrature.
template <typename T = double>
T adaptive_simpson_2d(std::function<T(T, T)> f, T ax, T bx, T ay, T by, T tol = T(1e-4), int max_depth = 10) {
  auto inner_integral = [&](T x) {
    auto fx = [&](T y) { return f(x, y); };
    return adaptive_simpson<T>(fx, ay, by, tol, max_depth);
  };
  return adaptive_simpson<T>(inner_integral, ax, bx, tol, max_depth);
}

/// @brief Generates an integral function F(x) = \int_{start_point}^{x} f(t) dt.
template <typename T = double>
numc::func<T> integral(numc::func<T> f, T start_point = T(0.0), T tol = T(1e-6)) {
  return numc::func<T>([f, start_point, tol](T x) {
    if (x == start_point) return T(0.0);
    if (x < start_point) {
      return -adaptive_simpson(f, x, start_point, tol);
    }
    return adaptive_simpson(f, start_point, x, tol);
  });
}

/// @brief Computes the integral of f(x) from a to b using the Monte Carlo method.
/// @return A pair containing the integral value and the estimated error.
template <typename T = double>
std::pair<T, T> monte_carlo_integrate(numc::func<T> f, T a, T b, size_t num_samples) {
  std::mt19937_64 rng(std::random_device{}());
  std::uniform_real_distribution<T> dist(a, b);
  T sum = T(0.0);
  T sum_sq = T(0.0);
  for (size_t i = 0; i < num_samples; ++i) {
    T val = f(dist(rng));
    sum += val;
    sum_sq += val * val;
  }
  T mean = sum / T(num_samples);
  T variance = (sum_sq / T(num_samples)) - (mean * mean);
  T integral_val = (b - a) * mean;
  T error = (b - a) * std::sqrt(variance / T(num_samples));
  return {integral_val, error};
}

/// @brief Computes the nodal abscissas and weights for Gauss-Legendre quadrature.
template <typename T = double>
std::pair<vector<T>, vector<T>> gauss_nodes(int m, T tol = T(1e-9)) {
  auto legendre = [](T t, int m_val) -> std::pair<T, T> {
    T p0 = T(1.0);
    T p1 = t;
    for (int k = 1; k < m_val; ++k) {
      T p = ((T(2.0) * k + T(1.0)) * t * p1 - T(k) * p0) / T(1.0 + k);
      p0 = p1;
      p1 = p;
    }
    T dp = T(m_val) * (p0 - t * p1) / (T(1.0) - t * t);
    return {p1, dp};
  };

  vector<T> x(m);
  vector<T> A(m);
  int nRoots = (m + 1) / 2;
  T pi_val = T(3.14159265358979323846L);
  for (int i = 0; i < nRoots; ++i) {
    T t = std::cos(pi_val * (T(i) + T(0.75)) / (T(m) + T(0.5)));
    for (int j = 0; j < 30; ++j) {
      auto [p, dp] = legendre(t, m);
      T dt = -p / dp;
      t += dt;
      if (std::abs(dt) < tol) {
        x[i] = t;
        x[m - i - 1] = -t;
        A[i] = T(2.0) / ((T(1.0) - t * t) * dp * dp);
        A[m - i - 1] = A[i];
        break;
      }
    }
  }
  return {x, A};
}

/// @brief Computes the integral of f(x) from a to b using Gauss-Legendre quadrature.
template <typename T = double>
T gauss_quad(numc::func<T> f, T a, T b, int m) {
  T c1 = (b + a) / T(2.0);
  T c2 = (b - a) / T(2.0);
  auto [x, A] = gauss_nodes(m);
  T sum = T(0.0);
  for (size_t i = 0; i < x.size(); ++i) {
    sum += A[i] * f(c1 + c2 * x[i]);
  }
  return c2 * sum;
}

/// @brief Integrates f(x, y) over a quadrilateral element.
template <typename T = double>
T gauss_quad2(std::function<T(T, T)> f, const vector<T>& x, const vector<T>& y, int m) {
  auto jac = [](const vector<T>& xv, const vector<T>& yv, T s, T t) -> T {
    T J00 = -(T(1.0) - t) * xv[0] + (T(1.0) - t) * xv[1] + (T(1.0) + t) * xv[2] - (T(1.0) + t) * xv[3];
    T J01 = -(T(1.0) - t) * yv[0] + (T(1.0) - t) * yv[1] + (T(1.0) + t) * yv[2] - (T(1.0) + t) * yv[3];
    T J10 = -(T(1.0) - s) * xv[0] - (T(1.0) + s) * xv[1] + (T(1.0) + s) * xv[2] + (T(1.0) - s) * xv[3];
    T J11 = -(T(1.0) - s) * yv[0] - (T(1.0) + s) * yv[1] + (T(1.0) + s) * yv[2] + (T(1.0) - s) * yv[3];
    return (J00 * J11 - J01 * J10) / T(16.0);
  };
  auto map_coords = [](const vector<T>& xv, const vector<T>& yv, T s, T t) -> std::pair<T, T> {
    T N0 = (T(1.0) - s) * (T(1.0) - t) / T(4.0);
    T N1 = (T(1.0) + s) * (T(1.0) - t) / T(4.0);
    T N2 = (T(1.0) + s) * (T(1.0) + t) / T(4.0);
    T N3 = (T(1.0) - s) * (T(1.0) + t) / T(4.0);
    T xCoord = N0 * xv[0] + N1 * xv[1] + N2 * xv[2] + N3 * xv[3];
    T yCoord = N0 * yv[0] + N1 * yv[1] + N2 * yv[2] + N3 * yv[3];
    return {xCoord, yCoord};
  };

  auto [s_nodes, A_w] = gauss_nodes(m);
  T sum = T(0.0);
  for (int i = 0; i < m; ++i) {
    for (int j = 0; j < m; ++j) {
      auto [xCoord, yCoord] = map_coords(x, y, s_nodes[i], s_nodes[j]);
      sum += A_w[i] * A_w[j] * jac(x, y, s_nodes[i], s_nodes[j]) * f(xCoord, yCoord);
    }
  }
  return sum;
}

/// @brief Integrates f(x, y) over a triangular element using the cubic formula.
template <typename T = double>
T triangle_quad(std::function<T(T, T)> f, const vector<T>& xc, const vector<T>& yc) {
  T alpha[4][3] = {
      {T(1.0) / T(3.0), T(1.0) / T(3.0), T(1.0) / T(3.0)},
      {T(0.2), T(0.2), T(0.6)},
      {T(0.6), T(0.2), T(0.2)},
      {T(0.2), T(0.6), T(0.2)}};
  T W[4] = {
      T(-27.0) / T(48.0),
      T(25.0) / T(48.0),
      T(25.0) / T(48.0),
      T(25.0) / T(48.0)};
  vector<T> x(4);
  vector<T> y(4);
  for (int i = 0; i < 4; ++i) {
    x[i] = alpha[i][0] * xc[0] + alpha[i][1] * xc[1] + alpha[i][2] * xc[2];
    y[i] = alpha[i][0] * yc[0] + alpha[i][1] * yc[1] + alpha[i][2] * yc[2];
  }
  T A = (xc[1] * yc[2] - xc[2] * yc[1] -
         xc[0] * yc[2] + xc[2] * yc[0] +
         xc[0] * yc[1] - xc[1] * yc[0]) /
        T(2.0);
  T sum = T(0.0);
  for (int i = 0; i < 4; ++i) {
    sum += W[i] * f(x[i], y[i]);
  }
  return A * sum;
}

/// @}

}  // namespace analysis

template <typename T>
T func<T>::integral(T a, T b, T tol) const {
  return analysis::adaptive_simpson<T>(*this, a, b, tol);
}

template <typename T>
T func<T>::integral_2d(T ax, T bx, T ay, T by, T tol) const {
  auto wrapper = [this](T x, T y) { return (*this)(x, y); };
  return analysis::adaptive_simpson_2d<T>(wrapper, ax, bx, ay, by, tol);
}

template <typename T>
std::pair<T, T> func<T>::integral_mc(const vector<T>& a, const vector<T>& b, size_t n_samples) const {
  if (a.size() != b.size() || a.size() == 0) {
    throw std::invalid_argument("integral_mc: bounds must have same non-zero size");
  }
  size_t dim = a.size();
  T v_box = T(1.0);
  for (size_t i = 0; i < dim; ++i) {
    v_box *= (b[i] - a[i]);
  }

  T sum = T(0.0);
  T sum_sq = T(0.0);

  // use simple random generator
  std::mt19937 gen(42); // fixed seed for reproducibility, or random_device
  std::vector<std::uniform_real_distribution<T>> dists;
  for (size_t i = 0; i < dim; ++i) {
    dists.emplace_back(a[i], b[i]);
  }

  vector<T> pt(dim);
  for (size_t i = 0; i < n_samples; ++i) {
    for (size_t j = 0; j < dim; ++j) {
      pt[j] = dists[j](gen);
    }
    T val = (*this)(pt);
    sum += val;
    sum_sq += val * val;
  }

  T mean = sum / T(n_samples);
  T variance = (sum_sq / T(n_samples)) - (mean * mean);
  T error = v_box * std::sqrt(variance / T(n_samples));
  T integral = v_box * mean;

  return {integral, error};
}

}  // namespace numc
