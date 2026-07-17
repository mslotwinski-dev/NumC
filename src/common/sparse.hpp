#pragma once

#include "../common/vector.hpp"
#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {

/// @addtogroup structures
/// @{

/// @brief Sparse matrix in Compressed Sparse Row (CSR) format.
template <typename T = double>
class sparse_matrix {
  static_assert(std::is_floating_point<T>::value, "numc::sparse_matrix supports floating-point types only!");

 private:
  size_t _rows, _cols;
  std::vector<T> _values;
  std::vector<size_t> _col_indices;
  std::vector<size_t> _row_ptr;

 public:
  /// @brief Constructs an empty sparse matrix.
  sparse_matrix(size_t rows, size_t cols)
      : _rows(rows), _cols(cols), _row_ptr(rows + 1, 0) {}

  /// @brief Constructs a sparse matrix from an initializer list of triplets.
  sparse_matrix(size_t rows, size_t cols, std::initializer_list<std::tuple<size_t, size_t, T>> triplets)
      : sparse_matrix(from_triplets(rows, cols, std::vector<std::tuple<size_t, size_t, T>>(triplets))) {}

  /// @brief Returns the number of rows.
  size_t rows() const { return _rows; }

  /// @brief Returns the number of columns.
  size_t cols() const { return _cols; }

  /// @brief Returns the number of non-zero elements.
  size_t nnz() const { return _values.size(); }

  /// @brief Constructs a sparse matrix from triplets (row, col, value).
  static sparse_matrix from_triplets(
      size_t rows, size_t cols,
      const std::vector<std::tuple<size_t, size_t, T>>& triplets) {
    sparse_matrix m(rows, cols);

    // Sort triplets by row then column
    auto sorted = triplets;
    std::sort(sorted.begin(), sorted.end(),
        [](const auto& a, const auto& b) {
          if (std::get<0>(a) != std::get<0>(b)) return std::get<0>(a) < std::get<0>(b);
          return std::get<1>(a) < std::get<1>(b);
        });

    m._values.reserve(sorted.size());
    m._col_indices.reserve(sorted.size());

    size_t current_row = 0;
    m._row_ptr[0] = 0;
    for (const auto& [r, c, v] : sorted) {
      while (current_row < r) {
        ++current_row;
        m._row_ptr[current_row] = m._values.size();
      }
      m._values.push_back(v);
      m._col_indices.push_back(c);
    }
    for (size_t i = current_row + 1; i <= rows; ++i) {
      m._row_ptr[i] = m._values.size();
    }
    return m;
  }

  /// @brief Constructs a sparse identity matrix.
  static sparse_matrix eye(size_t n) {
    std::vector<std::tuple<size_t, size_t, T>> triplets;
    for (size_t i = 0; i < n; ++i) triplets.push_back({i, i, T(1.0)});
    return from_triplets(n, n, triplets);
  }

  /// @brief Sparse matrix-vector multiplication: y = A * x.
  vector<T> operator*(const vector<T>& x) const {
    vector<T> y(_rows);
    for (size_t i = 0; i < _rows; ++i) {
      T sum = T(0.0);
      for (size_t k = _row_ptr[i]; k < _row_ptr[i + 1]; ++k) {
        sum += _values[k] * x[_col_indices[k]];
      }
      y[i] = sum;
    }
    return y;
  }

  /// @brief Solves Ax = b using the Conjugate Gradient method (for SPD matrices).
  vector<T> solve_cg(const vector<T>& b, T tol = T(1e-8), int max_iter = 10000) const {
    vector<T> x(_cols);  // initial guess = 0
    vector<T> r = b;     // r = b - A*x = b
    vector<T> p = r;
    T rs_old = r * r;

    for (int i = 0; i < max_iter; ++i) {
      vector<T> Ap = (*this) * p;
      T pAp = p * Ap;
      if (std::abs(pAp) < T(1e-20)) break;
      T alpha = rs_old / pAp;
      x = x + alpha * p;
      r = r - alpha * Ap;
      T rs_new = r * r;
      if (std::sqrt(rs_new) < tol) return x;
      p = r + (rs_new / rs_old) * p;
      rs_old = rs_new;
    }
    Log::Warn("Sparse CG did not converge.");
    return x;
  }

  /// @brief Solves Ax = b using Jacobi iterative method.
  vector<T> solve_jacobi(const vector<T>& b, T tol = T(1e-8), int max_iter = 10000) const {
    vector<T> x(_cols);
    vector<T> x_new(_cols);

    for (int iter = 0; iter < max_iter; ++iter) {
      for (size_t i = 0; i < _rows; ++i) {
        T sum = T(0.0);
        T diag = T(0.0);
        for (size_t k = _row_ptr[i]; k < _row_ptr[i + 1]; ++k) {
          if (_col_indices[k] == i) diag = _values[k];
          else sum += _values[k] * x[_col_indices[k]];
        }
        if (std::abs(diag) > T(1e-20)) {
          x_new[i] = (b[i] - sum) / diag;
        }
      }
      T norm_diff = T(0.0);
      for (size_t i = 0; i < _cols; ++i) {
        T d = x_new[i] - x[i];
        norm_diff += d * d;
      }
      x = x_new;
      if (std::sqrt(norm_diff) < tol) return x;
    }
    Log::Warn("Sparse Jacobi did not converge.");
    return x;
  }

  /// @brief Solves Ax = b using Gauss-Seidel iterative method.
  vector<T> solve_gauss_seidel(const vector<T>& b, T tol = T(1e-8), int max_iter = 10000) const {
    vector<T> x(_cols);

    for (int iter = 0; iter < max_iter; ++iter) {
      T max_diff = T(0.0);
      for (size_t i = 0; i < _rows; ++i) {
        T sum = T(0.0);
        T diag = T(0.0);
        for (size_t k = _row_ptr[i]; k < _row_ptr[i + 1]; ++k) {
          if (_col_indices[k] == i) diag = _values[k];
          else sum += _values[k] * x[_col_indices[k]];
        }
        if (std::abs(diag) > T(1e-20)) {
          T x_old = x[i];
          x[i] = (b[i] - sum) / diag;
          max_diff = std::max(max_diff, std::abs(x[i] - x_old));
        }
      }
      if (max_diff < tol) return x;
    }
    Log::Warn("Sparse Gauss-Seidel did not converge.");
    return x;
  }

  /// @brief Prints the sparse matrix in a readable format.
  friend std::ostream& operator<<(std::ostream& os, const sparse_matrix& m) {
    os << "Sparse matrix (" << m._rows << " x " << m._cols << "), nnz = " << m.nnz() << "\n";
    for (size_t i = 0; i < m._rows; ++i) {
      for (size_t k = m._row_ptr[i]; k < m._row_ptr[i + 1]; ++k) {
        os << "  (" << i << ", " << m._col_indices[k] << ") = " << m._values[k] << "\n";
      }
    }
    return os;
  }
};

/// @}

}  // namespace numc
