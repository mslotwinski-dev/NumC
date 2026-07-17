#pragma once

#include "../../common/vector.hpp"
#include "../../common/function.hpp"
#include "../../inc.hpp"
#include "../../utility/log.hpp"

namespace numc {
namespace optimalization {

/// @addtogroup optimize
/// @{

/// @brief Brackets a minimum of the function f(x) starting from x_start with step h.
template <typename T = double>
std::pair<T, T> gold_bracket(numc::func<T> f, T x_start, T h) {
  T c = T(1.618033988749895);
  T x1 = x_start;
  T f1 = f(x1);
  T x2 = x1 + h;
  T f2 = f(x2);
  if (f2 > f1) {
    h = -h;
    x2 = x1 + h;
    f2 = f(x2);
    if (f2 > f1) {
      return {x2, x1 - h};
    }
  }
  for (int i = 0; i < 100; ++i) {
    h = c * h;
    T x3 = x2 + h;
    T f3 = f(x3);
    if (f3 > f2) {
      return {x1, x3};
    }
    x1 = x2;
    x2 = x3;
    f1 = f2;
    f2 = f3;
  }
  Log::Warn("Bracket did not find a minimum.");
  return {x1, x2};
}

/// @brief Finds the minimum of f(x) in the bracket [a, b] using Golden Section search.
template <typename T = double>
std::pair<T, T> gold_search(numc::func<T> f, T a, T b, T tol = T(1e-6)) {
  T R = T(0.618033988749895);
  T C = T(1.0) - R;
  T diff = std::abs(b - a);
  if (diff < tol) {
    T mid = T(0.5) * (a + b);
    return {mid, f(mid)};
  }
  T nIter_val = T(-2.078087) * std::log(tol / diff);
  int nIter = static_cast<int>(std::ceil(nIter_val));
  if (nIter < 1) nIter = 1;

  T x1 = R * a + C * b;
  T x2 = C * a + R * b;
  T f1 = f(x1);
  T f2 = f(x2);

  for (int i = 0; i < nIter; ++i) {
    if (f1 > f2) {
      a = x1;
      x1 = x2;
      f1 = f2;
      x2 = C * a + R * b;
      f2 = f(x2);
    } else {
      b = x2;
      x2 = x1;
      f2 = f1;
      x1 = R * a + C * b;
      f1 = f(x1);
    }
  }
  if (f1 < f2) return {x1, f1};
  return {x2, f2};
}

/// @brief Minimizes a function F(x) of n variables using Powell's direction set method.
template <typename T = double>
std::pair<vector<T>, int> powell(std::function<T(const vector<T>&)> F, vector<T> x, T h = T(0.1), T tol = T(1e-6)) {
  size_t n = x.size();
  vector<T> df(n);
  std::vector<vector<T>> u(n);
  for (size_t i = 0; i < n; ++i) {
    u[i] = vector<T>(n);
    u[i][i] = T(1.0);
  }

  for (int j = 0; j < 30; ++j) {
    vector<T> xOld = x;
    T fOld = F(xOld);

    for (size_t i = 0; i < n; ++i) {
      vector<T> v = u[i];
      numc::func<T> f_line([&F, &x, &v](T s) -> T { return F(x + s * v); });
      auto [a, b] = gold_bracket<T>(f_line, T(0.0), h);
      auto [s, fMin] = gold_search<T>(f_line, a, b);
      df[i] = fOld - fMin;
      fOld = fMin;
      x = x + s * v;
    }

    vector<T> v = x - xOld;
    numc::func<T> f_line([&F, &x, &v](T s) -> T { return F(x + s * v); });
    auto [a, b] = gold_bracket<T>(f_line, T(0.0), h);
    auto [s, fLast] = gold_search<T>(f_line, a, b);
    x = x + s * v;

    T diff_norm = (x - xOld).norm();
    if (diff_norm / std::sqrt(T(n)) < tol) {
      return {x, j + 1};
    }

    size_t iMax = 0;
    T max_decrease = df[0];
    for (size_t i = 1; i < n; ++i) {
      if (df[i] > max_decrease) {
        max_decrease = df[i];
        iMax = i;
      }
    }

    for (size_t i = iMax; i < n - 1; ++i) {
      u[i] = u[i + 1];
    }
    u[n - 1] = v;
  }
  Log::Warn("Powell did not converge.");
  return {x, 30};
}

/// @brief Minimizes a function F(x) of n variables using the Fletcher-Reeves conjugate gradient method.
template <typename T = double>
std::pair<vector<T>, int> fletcher_reeves(
    std::function<T(const vector<T>&)> F,
    std::function<vector<T>(const vector<T>&)> gradF,
    vector<T> x,
    T h = T(0.01),
    T tol = T(1e-6)) {
  size_t n = x.size();
  vector<T> g0 = -gradF(x);
  vector<T> v = g0;
  T F0 = F(x);

  for (int i = 0; i < 200; ++i) {
    numc::func<T> f_line([&F, &x, &v](T s) -> T { return F(x + s * v); });
    auto [a, b] = gold_bracket<T>(f_line, T(0.0), h);
    auto [s, fMin] = gold_search<T>(f_line, a, b);
    x = x + s * v;
    T F1 = F(x);
    vector<T> g1 = -gradF(x);

    if (g1.norm() <= tol || std::abs(F0 - F1) < tol) {
      return {x, i + 1};
    }

    T g1_minus_g0_dot_g1 = (g1 - g0) * g1;
    T g0_dot_g0 = g0 * g0;
    T gamma = T(0.0);
    if (std::abs(g0_dot_g0) > T(1e-20)) {
      gamma = g1_minus_g0_dot_g1 / g0_dot_g0;
    }

    v = g1 + gamma * v;
    g0 = g1;
    F0 = F1;
  }
  Log::Warn("Fletcher-Reeves did not converge.");
  return {x, 200};
}

/// @brief Nelder-Mead simplex method for multidimensional minimization (derivative-free).
template <typename T = double>
std::pair<vector<T>, int> nelder_mead(
    std::function<T(const vector<T>&)> F,
    vector<T> x0,
    T step = T(0.5),
    T tol = T(1e-6),
    int max_iter = 1000) {
  size_t n = x0.size();
  size_t np = n + 1;

  // Build initial simplex
  std::vector<vector<T>> simplex(np, x0);
  for (size_t i = 0; i < n; ++i) {
    simplex[i + 1][i] += step;
  }

  std::vector<T> fvals(np);
  for (size_t i = 0; i < np; ++i) {
    fvals[i] = F(simplex[i]);
  }

  T alpha = T(1.0), gamma_c = T(2.0), rho = T(0.5), sigma = T(0.5);

  for (int iter = 0; iter < max_iter; ++iter) {
    // Sort vertices
    std::vector<size_t> order(np);
    std::iota(order.begin(), order.end(), 0);
    std::sort(order.begin(), order.end(), [&](size_t a, size_t b) { return fvals[a] < fvals[b]; });
    std::vector<vector<T>> sorted_simplex(np);
    std::vector<T> sorted_fvals(np);
    for (size_t i = 0; i < np; ++i) {
      sorted_simplex[i] = simplex[order[i]];
      sorted_fvals[i] = fvals[order[i]];
    }
    simplex = sorted_simplex;
    fvals = sorted_fvals;

    // Convergence check
    T frange = std::abs(fvals[np - 1] - fvals[0]);
    if (frange < tol) return {simplex[0], iter + 1};

    // Centroid (excluding worst)
    vector<T> centroid(n);
    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        centroid[j] += simplex[i][j];
      }
    }
    for (size_t j = 0; j < n; ++j) centroid[j] /= T(n);

    // Reflection
    vector<T> xr = centroid + alpha * (centroid - simplex[np - 1]);
    T fr = F(xr);

    if (fr < fvals[n - 1] && fr >= fvals[0]) {
      simplex[np - 1] = xr;
      fvals[np - 1] = fr;
      continue;
    }

    // Expansion
    if (fr < fvals[0]) {
      vector<T> xe = centroid + gamma_c * (xr - centroid);
      T fe = F(xe);
      if (fe < fr) {
        simplex[np - 1] = xe;
        fvals[np - 1] = fe;
      } else {
        simplex[np - 1] = xr;
        fvals[np - 1] = fr;
      }
      continue;
    }

    // Contraction
    vector<T> xc = centroid + rho * (simplex[np - 1] - centroid);
    T fc = F(xc);
    if (fc < fvals[np - 1]) {
      simplex[np - 1] = xc;
      fvals[np - 1] = fc;
      continue;
    }

    // Shrink
    for (size_t i = 1; i < np; ++i) {
      simplex[i] = simplex[0] + sigma * (simplex[i] - simplex[0]);
      fvals[i] = F(simplex[i]);
    }
  }
  Log::Warn("Nelder-Mead did not converge.");
  return {simplex[0], max_iter};
}

/// @brief Gradient descent with optional momentum for multidimensional minimization.
template <typename T = double>
std::pair<vector<T>, int> gradient_descent(
    std::function<T(const vector<T>&)> F,
    std::function<vector<T>(const vector<T>&)> gradF,
    vector<T> x,
    T lr = T(0.01),
    T momentum = T(0.0),
    T tol = T(1e-6),
    int max_iter = 10000) {
  size_t n = x.size();
  vector<T> v(n);  // velocity for momentum

  for (int i = 0; i < max_iter; ++i) {
    vector<T> g = gradF(x);
    if (g.norm() < tol) return {x, i + 1};

    v = momentum * v - lr * g;
    x = x + v;
  }
  Log::Warn("Gradient descent did not converge.");
  return {x, max_iter};
}

/// @brief BFGS quasi-Newton method for multidimensional minimization.
template <typename T = double>
std::pair<vector<T>, int> bfgs(
    std::function<T(const vector<T>&)> F,
    std::function<vector<T>(const vector<T>&)> gradF,
    vector<T> x,
    T tol = T(1e-6),
    int max_iter = 200) {
  size_t n = x.size();
  // Initialize H as identity (inverse Hessian approximation)
  tensor<T> H = tensor<T>::eye(n);

  vector<T> g = gradF(x);

  for (int iter = 0; iter < max_iter; ++iter) {
    if (g.norm() < tol) return {x, iter + 1};

    // Search direction: p = -H * g
    vector<T> p(n);
    for (size_t i = 0; i < n; ++i) {
      T sum = T(0.0);
      for (size_t j = 0; j < n; ++j) {
        sum += H(i, j) * g[j];
      }
      p[i] = -sum;
    }

    // Line search using golden section
    numc::func<T> f_line([&F, &x, &p](T s) -> T { return F(x + s * p); });
    auto [a, b] = gold_bracket<T>(f_line, T(0.0), T(0.1));
    auto [alpha_val, fMin] = gold_search<T>(f_line, a, b);

    vector<T> s_vec = alpha_val * p;
    vector<T> x_new = x + s_vec;
    vector<T> g_new = gradF(x_new);
    vector<T> y_vec = g_new - g;

    T sy = s_vec * y_vec;  // dot product
    if (std::abs(sy) < T(1e-20)) {
      x = x_new;
      g = g_new;
      continue;
    }

    // BFGS update of H
    // H = (I - rho*s*y^T) * H * (I - rho*y*s^T) + rho*s*s^T
    T rho = T(1.0) / sy;

    // Compute Hy
    vector<T> Hy(n);
    for (size_t i = 0; i < n; ++i) {
      T sum = T(0.0);
      for (size_t j = 0; j < n; ++j) {
        sum += H(i, j) * y_vec[j];
      }
      Hy[i] = sum;
    }

    T yHy = y_vec * Hy;

    for (size_t i = 0; i < n; ++i) {
      for (size_t j = 0; j < n; ++j) {
        H(i, j) += rho * ((T(1.0) + rho * yHy) * s_vec[i] * s_vec[j]
                          - Hy[i] * s_vec[j] - s_vec[i] * Hy[j]);
      }
    }

    x = x_new;
    g = g_new;
  }
  Log::Warn("BFGS did not converge.");
  return {x, max_iter};
}

/// @brief Simulated annealing for global optimization (stochastic, derivative-free).
template <typename T = double>
std::pair<vector<T>, T> simulated_annealing(
    std::function<T(const vector<T>&)> F,
    vector<T> x0,
    T temp_init = T(100.0),
    T temp_min = T(1e-8),
    T cooling = T(0.995),
    int iter_per_temp = 100,
    unsigned seed = 42) {
  size_t n = x0.size();
  std::mt19937 rng(seed);
  std::uniform_real_distribution<T> uni(T(0.0), T(1.0));
  std::normal_distribution<T> norm(T(0.0), T(1.0));

  vector<T> x_best = x0;
  T f_best = F(x0);
  vector<T> x_curr = x0;
  T f_curr = f_best;
  T temp = temp_init;

  while (temp > temp_min) {
    for (int i = 0; i < iter_per_temp; ++i) {
      // Perturb
      vector<T> x_new = x_curr;
      for (size_t j = 0; j < n; ++j) {
        x_new[j] += temp * norm(rng);
      }
      T f_new = F(x_new);
      T delta = f_new - f_curr;

      if (delta < T(0.0) || uni(rng) < std::exp(-delta / temp)) {
        x_curr = x_new;
        f_curr = f_new;
        if (f_curr < f_best) {
          x_best = x_curr;
          f_best = f_curr;
        }
      }
    }
    temp *= cooling;
  }
  return {x_best, f_best};
}

/// @}

}  // namespace optimalization
}  // namespace numc
