#pragma once

#include "../common/vector.hpp"
#include "../common/function.hpp"
#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {
namespace analysis {

/// @addtogroup analysis
/// @{

/// @brief Computes the numerical limit of f(x) as x -> x0.
template <typename T = double>
T limit(numc::func<T> f, T x0, T h = T(0.01), T tol = T(1e-10), int max_steps = 30) {
  T prev = f(x0 + h);
  for (int i = 0; i < max_steps; ++i) {
    h /= T(2.0);
    T curr = f(x0 + h);
    if (std::abs(curr - prev) < tol) return curr;
    prev = curr;
  }
  return prev;
}

/// @brief Computes the numerical limit from the left: lim_{x->x0^-} f(x).
template <typename T = double>
T limit_left(numc::func<T> f, T x0, T h = T(0.01), T tol = T(1e-10)) {
  return limit(f, x0, -std::abs(h), tol);
}

/// @brief Computes the numerical limit from the right: lim_{x->x0^+} f(x).
template <typename T = double>
T limit_right(numc::func<T> f, T x0, T h = T(0.01), T tol = T(1e-10)) {
  return limit(f, x0, std::abs(h), tol);
}

/// @brief Computes the numerical partial sum of a series: sum_{n=start}^{N} f(n).
template <typename T = double>
T series_sum(numc::func<T> f, int start, int N) {
  T sum = T(0.0);
  for (int n = start; n <= N; ++n) sum += f(T(n));
  return sum;
}

/// @brief Computes the sum of a convergent infinite series using Aitken's delta-squared acceleration.
template <typename T = double>
T series_sum_convergent(numc::func<T> f, int start = 0, T tol = T(1e-10), int max_terms = 100000) {
  T S0 = T(0.0), S1 = T(0.0), S2 = T(0.0);
  S0 = f(T(start));
  S1 = S0 + f(T(start + 1));
  for (int n = start + 2; n < start + max_terms; ++n) {
    S2 = S1 + f(T(n));
    // Aitken acceleration
    T denom = S2 - T(2.0) * S1 + S0;
    if (std::abs(denom) > T(1e-20)) {
      T accelerated = S0 - (S1 - S0) * (S1 - S0) / denom;
      if (std::abs(S2 - S1) < tol) return accelerated;
    } else {
      if (std::abs(S2 - S1) < tol) return S2;
    }
    S0 = S1;
    S1 = S2;
  }
  Log::Warn("Series did not converge.");
  return S2;
}

/// @brief Computes the numerical product: prod_{n=start}^{N} f(n).
template <typename T = double>
T series_product(numc::func<T> f, int start, int N) {
  T prod = T(1.0);
  for (int n = start; n <= N; ++n) prod *= f(T(n));
  return prod;
}

/// @brief Solves equation f(x) = target by generating a table of values and interpolating.
template <typename T = double>
T find_value(numc::func<T> f, T target, T a, T b, T tol = T(1e-9)) {
  // Use bisection internally
  T fa = f(a) - target;
  T fb = f(b) - target;
  if (fa * fb > T(0.0)) {
    Log::Warn("find_value: target not bracketed in [a, b].");
    return std::numeric_limits<T>::quiet_NaN();
  }
  for (int i = 0; i < 100; ++i) {
    T c = (a + b) / T(2.0);
    T fc = f(c) - target;
    if (std::abs(fc) < tol || (b - a) / T(2.0) < tol) return c;
    if (fa * fc < T(0.0)) { b = c; fb = fc; }
    else { a = c; fa = fc; }
  }
  return (a + b) / T(2.0);
}

/// @brief Computes the n-th order Taylor polynomial coefficients of f around x0.
/// Returns vector of coefficients [f(x0), f'(x0)/1!, f''(x0)/2!, ...].
template <typename T = double>
vector<T> taylor_coeffs(numc::func<T> f, T x0, int n, T h = T(1e-3)) {
  vector<T> coeffs(n + 1);
  coeffs[0] = f(x0);

  // Use finite differences for higher derivatives
  // n-th derivative using central differences
  for (int k = 1; k <= n; ++k) {
    // k-th derivative via recursion of central differences
    T deriv = T(0.0);
    T sign = T(1.0);
    T binom = T(1.0);
    for (int j = 0; j <= k; ++j) {
      deriv += sign * binom * f(x0 + (T(k) / T(2.0) - T(j)) * h);
      sign = -sign;
      binom *= T(k - j) / T(j + 1);
    }
    deriv /= std::pow(h, k);
    // factorial
    T fact = T(1.0);
    for (int j = 1; j <= k; ++j) fact *= T(j);
    coeffs[k] = deriv / fact;
  }
  return coeffs;
}

/// @brief Evaluates a Taylor polynomial given its coefficients at point x around x0.
template <typename T = double>
T taylor_eval(const vector<T>& coeffs, T x0, T x) {
  T dx = x - x0;
  T result = T(0.0);
  T dx_pow = T(1.0);
  for (size_t i = 0; i < coeffs.size(); ++i) {
    result += coeffs[i] * dx_pow;
    dx_pow *= dx;
  }
  return result;
}

/// @brief Computes the arc length of f(x) from a to b using numerical integration.
template <typename T = double>
T arc_length(numc::func<T> f, T a, T b, size_t n = 1000) {
  T h = (b - a) / T(n);
  T length = T(0.0);
  T x_prev = a;
  T y_prev = f(a);
  for (size_t i = 1; i <= n; ++i) {
    T x_curr = a + T(i) * h;
    T y_curr = f(x_curr);
    T dx = x_curr - x_prev;
    T dy = y_curr - y_prev;
    length += std::sqrt(dx * dx + dy * dy);
    x_prev = x_curr;
    y_prev = y_curr;
  }
  return length;
}

/// @brief Computes the curvature of f(x) at point x.
template <typename T = double>
T curvature(numc::func<T> f, T x, T h = T(1e-5)) {
  T fp = (f(x + h) - f(x - h)) / (T(2.0) * h);
  T fpp = (f(x + h) - T(2.0) * f(x) + f(x - h)) / (h * h);
  return std::abs(fpp) / std::pow(T(1.0) + fp * fp, T(1.5));
}

/// @brief Generates a table of function values. Returns {x_values, y_values}.
template <typename T = double>
std::pair<vector<T>, vector<T>> tabulate(numc::func<T> f, T a, T b, size_t n) {
  vector<T> xs(n);
  vector<T> ys(n);
  T h = (b - a) / T(n - 1);
  for (size_t i = 0; i < n; ++i) {
    xs[i] = a + T(i) * h;
    ys[i] = f(xs[i]);
  }
  return {xs, ys};
}

/// @}

}  // namespace analysis
}  // namespace numc
