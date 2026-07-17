#pragma once

#include "../common/tensor.hpp"
#include "../common/vector.hpp"
#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {
namespace linear_algebra {

/// @addtogroup linear_algebra
/// @{

/// @brief Solves the linear system Ax = b using Gauss elimination without pivoting.
template <typename T = double>
vector<T> gauss_elimin(tensor<T> a, vector<T> b) {
  size_t n = b.size();
  // Elimination phase
  for (size_t k = 0; k < n - 1; ++k) {
    for (size_t i = k + 1; i < n; ++i) {
      if (a(i, k) != T(0.0)) {
        T lam = a(i, k) / a(k, k);
        for (size_t j = k + 1; j < n; ++j) {
          a(i, j) -= lam * a(k, j);
        }
        b[i] -= lam * b[k];
      }
    }
  }
  // Back substitution
  for (int k = static_cast<int>(n) - 1; k >= 0; --k) {
    T sum = T(0.0);
    for (size_t j = k + 1; j < n; ++j) {
      sum += a(k, j) * b[j];
    }
    b[k] = (b[k] - sum) / a(k, k);
  }
  return b;
}

/// @brief Performs LU decomposition (Doolittle method) in-place on matrix A.
template <typename T = double>
tensor<T> lu_decomp(tensor<T> a) {
  size_t n = a.shape()[0];
  for (size_t k = 0; k < n - 1; ++k) {
    for (size_t i = k + 1; i < n; ++i) {
      if (a(i, k) != T(0.0)) {
        T lam = a(i, k) / a(k, k);
        for (size_t j = k + 1; j < n; ++j) {
          a(i, j) -= lam * a(k, j);
        }
        a(i, k) = lam;
      }
    }
  }
  return a;
}

/// @brief Solves the linear system LUx = b where LU is the decomposed matrix.
template <typename T = double>
vector<T> lu_solve(const tensor<T>& lu, vector<T> b) {
  size_t n = lu.shape()[0];
  for (size_t k = 1; k < n; ++k) {
    T sum = T(0.0);
    for (size_t j = 0; j < k; ++j) {
      sum += lu(k, j) * b[j];
    }
    b[k] -= sum;
  }
  for (int k = static_cast<int>(n) - 1; k >= 0; --k) {
    T sum = T(0.0);
    for (size_t j = k + 1; j < n; ++j) {
      sum += lu(k, j) * b[j];
    }
    b[k] = (b[k] - sum) / lu(k, k);
  }
  return b;
}

/// @brief Performs Cholesky decomposition A = LL^T. Returns the lower triangular matrix L.
template <typename T = double>
tensor<T> choleski(tensor<T> a) {
  size_t n = a.shape()[0];
  for (size_t k = 0; k < n; ++k) {
    T sum_diag = T(0.0);
    for (size_t j = 0; j < k; ++j) {
      sum_diag += a(k, j) * a(k, j);
    }
    T val = a(k, k) - sum_diag;
    if (val <= T(0.0)) {
      Log::Error("Cholesky decomposition failed: Matrix is not positive definite.");
      throw std::domain_error("Matrix is not positive definite.");
    }
    a(k, k) = std::sqrt(val);
    for (size_t i = k + 1; i < n; ++i) {
      T sum_row = T(0.0);
      for (size_t j = 0; j < k; ++j) {
        sum_row += a(i, j) * a(k, j);
      }
      a(i, k) = (a(i, k) - sum_row) / a(k, k);
    }
  }
  // Zero out the upper triangle
  for (size_t k = 1; k < n; ++k) {
    for (size_t i = 0; i < k; ++i) {
      a(i, k) = T(0.0);
    }
  }
  return a;
}

/// @brief Performs LU decomposition of a tridiagonal matrix defined by diagonals c, d, e.
template <typename T = double>
std::tuple<vector<T>, vector<T>, vector<T>> lu_decomp3(vector<T> c, vector<T> d, vector<T> e) {
  size_t n = d.size();
  for (size_t k = 1; k < n; ++k) {
    T lam = c[k - 1] / d[k - 1];
    d[k] -= lam * e[k - 1];
    c[k - 1] = lam;
  }
  return {c, d, e};
}

/// @brief Solves a tridiagonal system using the decomposed diagonals c, d, e.
template <typename T = double>
vector<T> lu_solve3(const vector<T>& c, const vector<T>& d, const vector<T>& e, vector<T> b) {
  size_t n = d.size();
  for (size_t k = 1; k < n; ++k) {
    b[k] -= c[k - 1] * b[k - 1];
  }
  b[n - 1] /= d[n - 1];
  for (int k = static_cast<int>(n) - 2; k >= 0; --k) {
    b[k] = (b[k] - e[k] * b[k + 1]) / d[k];
  }
  return b;
}

/// @brief Performs LU decomposition of a symmetric pentadiagonal matrix defined by diagonals d, e, f.
template <typename T = double>
std::tuple<vector<T>, vector<T>, vector<T>> lu_decomp5(vector<T> d, vector<T> e, vector<T> f) {
  size_t n = d.size();
  for (size_t k = 0; k < n - 2; ++k) {
    T lam = e[k] / d[k];
    d[k + 1] -= lam * e[k];
    e[k + 1] -= lam * f[k];
    e[k] = lam;

    lam = f[k] / d[k];
    d[k + 2] -= lam * f[k];
    f[k] = lam;
  }
  T lam = e[n - 2] / d[n - 2];
  d[n - 1] -= lam * e[n - 2];
  e[n - 2] = lam;
  return {d, e, f};
}

/// @brief Solves a symmetric pentadiagonal system using the decomposed diagonals d, e, f.
template <typename T = double>
vector<T> lu_solve5(const vector<T>& d, const vector<T>& e, const vector<T>& f, vector<T> b) {
  size_t n = d.size();
  b[1] -= e[0] * b[0];
  for (size_t k = 2; k < n; ++k) {
    b[k] -= e[k - 1] * b[k - 1] + f[k - 2] * b[k - 2];
  }
  b[n - 1] /= d[n - 1];
  b[n - 2] = b[n - 2] / d[n - 2] - e[n - 2] * b[n - 1];
  for (int k = static_cast<int>(n) - 3; k >= 0; --k) {
    b[k] = b[k] / d[k] - e[k] * b[k + 1] - f[k] * b[k + 2];
  }
  return b;
}

/// @brief Swaps two rows of a 2D tensor.
template <typename T = double>
void swap_rows(tensor<T>& a, size_t i, size_t j) {
  size_t cols = a.shape()[1];
  for (size_t k = 0; k < cols; ++k) {
    std::swap(a(i, k), a(j, k));
  }
}

/// @brief Swaps two elements of a vector.
template <typename T = double>
void swap_rows(vector<T>& v, size_t i, size_t j) {
  std::swap(v[i], v[j]);
}

/// @brief Swaps two columns of a 2D tensor.
template <typename T = double>
void swap_cols(tensor<T>& a, size_t i, size_t j) {
  size_t rows = a.shape()[0];
  for (size_t k = 0; k < rows; ++k) {
    std::swap(a(k, i), a(k, j));
  }
}

/// @brief Solves the linear system Ax = b using Gauss elimination with scaled row pivoting.
template <typename T = double>
vector<T> gauss_pivot(tensor<T> a, vector<T> b, T tol = T(1e-9)) {
  size_t n = b.size();
  vector<T> s(n);
  for (size_t i = 0; i < n; ++i) {
    T max_val = T(0.0);
    for (size_t j = 0; j < n; ++j) {
      max_val = std::max(max_val, std::abs(a(i, j)));
    }
    s[i] = max_val;
  }
  for (size_t k = 0; k < n - 1; ++k) {
    // Row interchange
    size_t p = k;
    T max_rel = T(0.0);
    for (size_t i = k; i < n; ++i) {
      T rel = std::abs(a(i, k)) / s[i];
      if (rel > max_rel) {
        max_rel = rel;
        p = i;
      }
    }
    if (std::abs(a(p, k)) < tol) {
      Log::Error("Matrix is singular in gauss_pivot.");
      throw std::runtime_error("Matrix is singular.");
    }
    if (p != k) {
      swap_rows(b, k, p);
      swap_rows(s, k, p);
      swap_rows(a, k, p);
    }
    // Elimination
    for (size_t i = k + 1; i < n; ++i) {
      if (a(i, k) != T(0.0)) {
        T lam = a(i, k) / a(k, k);
        for (size_t j = k + 1; j < n; ++j) {
          a(i, j) -= lam * a(k, j);
        }
        b[i] -= lam * b[k];
      }
    }
  }
  if (std::abs(a(n - 1, n - 1)) < tol) {
    Log::Error("Matrix is singular in gauss_pivot.");
    throw std::runtime_error("Matrix is singular.");
  }
  // Back substitution
  for (int k = static_cast<int>(n) - 1; k >= 0; --k) {
    T sum = T(0.0);
    for (size_t j = k + 1; j < n; ++j) {
      sum += a(k, j) * b[j];
    }
    b[k] = (b[k] - sum) / a(k, k);
  }
  return b;
}

/// @brief Performs LU decomposition with scaled row pivoting.
/// Returns the decomposed matrix and the row permutation sequence.
template <typename T = double>
std::pair<tensor<T>, std::vector<size_t>> lu_decomp_pivot(tensor<T> a, T tol = T(1e-9)) {
  size_t n = a.shape()[0];
  std::vector<size_t> seq(n);
  for (size_t i = 0; i < n; ++i) seq[i] = i;

  vector<T> s(n);
  for (size_t i = 0; i < n; ++i) {
    T max_val = T(0.0);
    for (size_t j = 0; j < n; ++j) {
      max_val = std::max(max_val, std::abs(a(i, j)));
    }
    s[i] = max_val;
  }

  for (size_t k = 0; k < n - 1; ++k) {
    size_t p = k;
    T max_rel = T(0.0);
    for (size_t i = k; i < n; ++i) {
      T rel = std::abs(a(i, k)) / s[i];
      if (rel > max_rel) {
        max_rel = rel;
        p = i;
      }
    }
    if (std::abs(a(p, k)) < tol) {
      Log::Error("Matrix is singular in lu_decomp_pivot.");
      throw std::runtime_error("Matrix is singular.");
    }
    if (p != k) {
      swap_rows(s, k, p);
      swap_rows(a, k, p);
      std::swap(seq[k], seq[p]);
    }
    // Elimination
    for (size_t i = k + 1; i < n; ++i) {
      if (a(i, k) != T(0.0)) {
        T lam = a(i, k) / a(k, k);
        for (size_t j = k + 1; j < n; ++j) {
          a(i, j) -= lam * a(k, j);
        }
        a(i, k) = lam;
      }
    }
  }
  return {a, seq};
}

/// @brief Solves the linear system LUx = b using the decomposed matrix with pivoting.
template <typename T = double>
vector<T> lu_solve_pivot(const tensor<T>& lu, const vector<T>& b, const std::vector<size_t>& seq) {
  size_t n = lu.shape()[0];
  vector<T> x(n);
  for (size_t i = 0; i < n; ++i) {
    x[i] = b[seq[i]];
  }
  for (size_t k = 1; k < n; ++k) {
    T sum = T(0.0);
    for (size_t j = 0; j < k; ++j) {
      sum += lu(k, j) * x[j];
    }
    x[k] -= sum;
  }
  for (int k = static_cast<int>(n) - 1; k >= 0; --k) {
    T sum = T(0.0);
    for (size_t j = k + 1; j < n; ++j) {
      sum += lu(k, j) * x[j];
    }
    x[k] = (x[k] - sum) / lu(k, k);
  }
  return x;
}

/// @brief Computes the inverse of a matrix using LU decomposition with pivoting.
template <typename T = double>
tensor<T> mat_inv(tensor<T> a) {
  size_t n = a.shape()[0];
  tensor<T> aInv = tensor<T>::eye(n);
  auto [lu, seq] = lu_decomp_pivot(a);
  for (size_t i = 0; i < n; ++i) {
    vector<T> col(n);
    for (size_t j = 0; j < n; ++j) {
      col[j] = aInv(j, i);
    }
    vector<T> sol = lu_solve_pivot(lu, col, seq);
    for (size_t j = 0; j < n; ++j) {
      aInv(j, i) = sol[j];
    }
  }
  return aInv;
}

/// @brief Solves a sparse linear system using Gauss-Seidel method with relaxation.
template <typename T = double>
std::tuple<vector<T>, int, T> gauss_seidel(
    std::function<vector<T>(const vector<T>&, T)> iter_eqs,
    vector<T> x,
    T tol = T(1e-9)) {
  T omega = T(1.0);
  int k = 10;
  int p = 1;
  T dx1 = T(0.0);
  T dx2 = T(0.0);
  int num_iter = 500;
  int final_iter = 0;

  for (int i = 1; i <= num_iter; ++i) {
    vector<T> x_old = x;
    x = iter_eqs(x, omega);

    T sum_diff = T(0.0);
    for (size_t j = 0; j < x.size(); ++j) {
      T diff = x[j] - x_old[j];
      sum_diff += diff * diff;
    }
    T dx = std::sqrt(sum_diff);

    if (dx < tol) {
      final_iter = i;
      return {x, final_iter, omega};
    }
    if (i == k) dx1 = dx;
    if (i == k + p) {
      dx2 = dx;
      T ratio = dx2 / dx1;
      T power = T(1.0) / T(p);
      omega = T(2.0) / (T(1.0) + std::sqrt(T(1.0) - std::pow(ratio, power)));
    }
  }
  Log::Warn("Gauss-Seidel failed to converge.");
  return {x, num_iter, omega};
}

/// @brief Solves a linear system Ax = b using the Conjugate Gradient method.
template <typename T = double>
std::pair<vector<T>, int> conj_grad(
    std::function<vector<T>(const vector<T>&)> Av,
    vector<T> x,
    const vector<T>& b,
    T tol = T(1e-9)) {
  size_t n = b.size();
  vector<T> Ax = Av(x);
  vector<T> r = b - Ax;
  vector<T> s = r;
  int it = 0;

  for (size_t i = 0; i < n; ++i) {
    vector<T> u = Av(s);
    T s_dot_r = s * r;
    T s_dot_u = s * u;
    if (std::abs(s_dot_u) < T(1e-12)) break;
    T alpha = s_dot_r / s_dot_u;
    x = x + alpha * s;

    vector<T> Av_x = Av(x);
    r = b - Av_x;

    T r_norm = r.norm();
    if (r_norm < tol) {
      it = i + 1;
      break;
    }

    T r_dot_u = r * u;
    T beta = -r_dot_u / s_dot_u;
    s = r + beta * s;
    it = i + 1;
  }
  return {x, it};
}

/// @}

}  // namespace linear_algebra
}  // namespace numc
