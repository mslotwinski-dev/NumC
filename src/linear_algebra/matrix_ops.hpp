#pragma once

#include "../common/tensor.hpp"
#include "../common/vector.hpp"
#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {
namespace linear_algebra {

/// @addtogroup linear_algebra
/// @{

/// @brief Computes the determinant of a square matrix using LU decomposition.
template <typename T = double>
T det(tensor<T> a) {
  size_t n = a.shape()[0];
  T d = T(1.0);
  for (size_t k = 0; k < n - 1; ++k) {
    // Partial pivoting
    size_t iMax = k;
    T maxVal = std::abs(a(k, k));
    for (size_t i = k + 1; i < n; ++i) {
      if (std::abs(a(i, k)) > maxVal) {
        maxVal = std::abs(a(i, k));
        iMax = i;
      }
    }
    if (iMax != k) {
      for (size_t j = 0; j < n; ++j) {
        T tmp = a(k, j);
        a(k, j) = a(iMax, j);
        a(iMax, j) = tmp;
      }
      d = -d;
    }
    if (std::abs(a(k, k)) < T(1e-15)) return T(0.0);
    d *= a(k, k);
    for (size_t i = k + 1; i < n; ++i) {
      T lam = a(i, k) / a(k, k);
      for (size_t j = k + 1; j < n; ++j) {
        a(i, j) -= lam * a(k, j);
      }
    }
  }
  d *= a(n - 1, n - 1);
  return d;
}

/// @brief Computes the trace of a square matrix (sum of diagonal elements).
template <typename T = double>
T trace(const tensor<T>& a) {
  size_t n = a.shape()[0];
  T sum = T(0.0);
  for (size_t i = 0; i < n; ++i) sum += a(i, i);
  return sum;
}

/// @brief Computes the transpose of a matrix.
template <typename T = double>
tensor<T> transpose(const tensor<T>& a) {
  size_t m = a.shape()[0], n = a.shape()[1];
  tensor<T> result({n, m});
  for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < n; ++j)
      result(j, i) = a(i, j);
  return result;
}

/// @brief Computes the matrix product C = A * B.
template <typename T = double>
tensor<T> matmul(const tensor<T>& a, const tensor<T>& b) {
  size_t m = a.shape()[0], k = a.shape()[1], n = b.shape()[1];
  tensor<T> c({m, n});
  for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < n; ++j) {
      T sum = T(0.0);
      for (size_t p = 0; p < k; ++p)
        sum += a(i, p) * b(p, j);
      c(i, j) = sum;
    }
  return c;
}

/// @brief Computes the matrix-vector product y = A * x.
template <typename T = double>
vector<T> matvec(const tensor<T>& a, const vector<T>& x) {
  size_t m = a.shape()[0], n = a.shape()[1];
  vector<T> y(m);
  for (size_t i = 0; i < m; ++i) {
    T sum = T(0.0);
    for (size_t j = 0; j < n; ++j)
      sum += a(i, j) * x[j];
    y[i] = sum;
  }
  return y;
}

/// @brief Computes the outer product of two vectors: M = u * v^T.
template <typename T = double>
tensor<T> outer(const vector<T>& u, const vector<T>& v) {
  size_t m = u.size(), n = v.size();
  tensor<T> result({m, n});
  for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < n; ++j)
      result(i, j) = u[i] * v[j];
  return result;
}

/// @brief Computes the Frobenius norm of a matrix.
template <typename T = double>
T frobenius_norm(const tensor<T>& a) {
  size_t m = a.shape()[0], n = a.shape()[1];
  T sum = T(0.0);
  for (size_t i = 0; i < m; ++i)
    for (size_t j = 0; j < n; ++j)
      sum += a(i, j) * a(i, j);
  return std::sqrt(sum);
}

/// @brief Computes the infinity norm of a matrix (max row sum).
template <typename T = double>
T inf_norm(const tensor<T>& a) {
  size_t m = a.shape()[0], n = a.shape()[1];
  T maxVal = T(0.0);
  for (size_t i = 0; i < m; ++i) {
    T rowSum = T(0.0);
    for (size_t j = 0; j < n; ++j)
      rowSum += std::abs(a(i, j));
    maxVal = std::max(maxVal, rowSum);
  }
  return maxVal;
}

/// @brief Computes the condition number of a matrix using the 1-norm estimation.
template <typename T = double>
T cond(tensor<T> a) {
  T normA = inf_norm(a);
  size_t n = a.shape()[0];
  // Compute inverse via Gauss-Jordan
  tensor<T> augmented({n, 2 * n});
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) augmented(i, j) = a(i, j);
    augmented(i, n + i) = T(1.0);
  }
  for (size_t k = 0; k < n; ++k) {
    // Partial pivot
    size_t iMax = k;
    for (size_t i = k + 1; i < n; ++i)
      if (std::abs(augmented(i, k)) > std::abs(augmented(iMax, k))) iMax = i;
    if (iMax != k)
      for (size_t j = 0; j < 2 * n; ++j) {
        T tmp = augmented(k, j);
        augmented(k, j) = augmented(iMax, j);
        augmented(iMax, j) = tmp;
      }
    T pivot = augmented(k, k);
    if (std::abs(pivot) < T(1e-15)) return std::numeric_limits<T>::infinity();
    for (size_t j = 0; j < 2 * n; ++j) augmented(k, j) /= pivot;
    for (size_t i = 0; i < n; ++i) {
      if (i == k) continue;
      T factor = augmented(i, k);
      for (size_t j = 0; j < 2 * n; ++j) augmented(i, j) -= factor * augmented(k, j);
    }
  }
  tensor<T> inv({n, n});
  for (size_t i = 0; i < n; ++i)
    for (size_t j = 0; j < n; ++j)
      inv(i, j) = augmented(i, n + j);
  T normInv = inf_norm(inv);
  return normA * normInv;
}

/// @brief Computes the rank of a matrix via row echelon form.
template <typename T = double>
size_t rank(tensor<T> a, T tol = T(1e-10)) {
  size_t m = a.shape()[0], n = a.shape()[1];
  size_t r = 0;
  for (size_t col = 0; col < n && r < m; ++col) {
    // Find pivot
    size_t pivot = m;
    for (size_t row = r; row < m; ++row) {
      if (std::abs(a(row, col)) > tol) { pivot = row; break; }
    }
    if (pivot == m) continue;
    // Swap
    if (pivot != r)
      for (size_t j = 0; j < n; ++j) {
        T tmp = a(r, j);
        a(r, j) = a(pivot, j);
        a(pivot, j) = tmp;
      }
    // Eliminate
    for (size_t i = r + 1; i < m; ++i) {
      T factor = a(i, col) / a(r, col);
      for (size_t j = col; j < n; ++j) a(i, j) -= factor * a(r, j);
    }
    ++r;
  }
  return r;
}

/// @brief QR decomposition using Gram-Schmidt orthogonalization (A = Q * R).
template <typename T = double>
std::pair<tensor<T>, tensor<T>> qr_decomp(const tensor<T>& A) {
  size_t m = A.shape()[0], n = A.shape()[1];
  tensor<T> Q({m, n});
  tensor<T> R({n, n});

  for (size_t j = 0; j < n; ++j) {
    // Copy column j of A
    vector<T> v(m);
    for (size_t i = 0; i < m; ++i) v[i] = A(i, j);

    for (size_t k = 0; k < j; ++k) {
      T dot = T(0.0);
      for (size_t i = 0; i < m; ++i) dot += Q(i, k) * v[i];
      R(k, j) = dot;
      for (size_t i = 0; i < m; ++i) v[i] -= dot * Q(i, k);
    }

    T norm_v = T(0.0);
    for (size_t i = 0; i < m; ++i) norm_v += v[i] * v[i];
    norm_v = std::sqrt(norm_v);
    R(j, j) = norm_v;

    if (norm_v > T(1e-15)) {
      for (size_t i = 0; i < m; ++i) Q(i, j) = v[i] / norm_v;
    }
  }
  return {Q, R};
}

/// @brief Singular Value Decomposition using one-sided Jacobi (A = U * S * V^T).
/// Returns {U, singular_values, V}.
template <typename T = double>
std::tuple<tensor<T>, vector<T>, tensor<T>> svd(const tensor<T>& A, T tol = T(1e-10), int max_iter = 100) {
  size_t m = A.shape()[0], n = A.shape()[1];
  // Work on A^T * A for the right singular vectors
  tensor<T> ATA = matmul(transpose(A), A);
  size_t p = n;

  // Initialize V as identity
  tensor<T> V = tensor<T>::eye(p);

  // Jacobi iterations on ATA
  for (int iter = 0; iter < max_iter; ++iter) {
    T off_diag = T(0.0);
    for (size_t i = 0; i < p; ++i)
      for (size_t j = i + 1; j < p; ++j)
        off_diag += ATA(i, j) * ATA(i, j);
    if (std::sqrt(off_diag) < tol) break;

    for (size_t i = 0; i < p - 1; ++i) {
      for (size_t j = i + 1; j < p; ++j) {
        if (std::abs(ATA(i, j)) < tol) continue;
        T theta = T(0.5) * std::atan2(T(2.0) * ATA(i, j), ATA(i, i) - ATA(j, j));
        T c = std::cos(theta), s = std::sin(theta);

        // Apply Givens rotation to ATA
        for (size_t k = 0; k < p; ++k) {
          T aik = ATA(k, i), ajk = ATA(k, j);
          ATA(k, i) = c * aik + s * ajk;
          ATA(k, j) = -s * aik + c * ajk;
        }
        for (size_t k = 0; k < p; ++k) {
          T aki = ATA(i, k), akj = ATA(j, k);
          ATA(i, k) = c * aki + s * akj;
          ATA(j, k) = -s * aki + c * akj;
        }

        // Accumulate V
        for (size_t k = 0; k < p; ++k) {
          T vki = V(k, i), vkj = V(k, j);
          V(k, i) = c * vki + s * vkj;
          V(k, j) = -s * vki + c * vkj;
        }
      }
    }
  }

  // Singular values are sqrt of eigenvalues of ATA
  vector<T> sigma(p);
  for (size_t i = 0; i < p; ++i) {
    sigma[i] = std::sqrt(std::max(ATA(i, i), T(0.0)));
  }

  // U = A * V * Sigma^(-1)
  tensor<T> U({m, p});
  tensor<T> AV = matmul(A, V);
  for (size_t j = 0; j < p; ++j) {
    if (sigma[j] > tol) {
      for (size_t i = 0; i < m; ++i) {
        U(i, j) = AV(i, j) / sigma[j];
      }
    }
  }

  return {U, sigma, V};
}

/// @brief Computes the pseudo-inverse (Moore-Penrose) of a matrix using SVD.
template <typename T = double>
tensor<T> pinv(const tensor<T>& A, T tol = T(1e-10)) {
  auto [U, sigma, V] = svd(A, tol);
  size_t m = U.shape()[0], n = V.shape()[0], k = sigma.size();

  // Sigma_plus = diag(1/sigma_i) for non-zero sigma_i
  tensor<T> result({n, m});
  for (size_t i = 0; i < k; ++i) {
    if (sigma[i] > tol) {
      T inv_s = T(1.0) / sigma[i];
      for (size_t row = 0; row < n; ++row)
        for (size_t col = 0; col < m; ++col)
          result(row, col) += V(row, i) * inv_s * U(col, i);
    }
  }
  return result;
}

/// @brief Solves the least-squares problem min ||Ax - b||^2 using SVD.
template <typename T = double>
vector<T> lstsq(const tensor<T>& A, const vector<T>& b, T tol = T(1e-10)) {
  tensor<T> Ap = pinv(A, tol);
  return matvec(Ap, b);
}

/// @}

}  // namespace linear_algebra
}  // namespace numc
