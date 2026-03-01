#pragma once

#include <cmath>
#include <iostream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {

/// @addtogroup structures
/// @{

/// @brief Class representing an N-dimensional multi-dimensional array (tensor).
/// Provides support for linear algebra operations like matrix multiplication, inversion, and element-wise arithmetic.
/// @tparam T The underlying floating-point type (e.g., float, double).
template <typename T = double>
class tensor {
  /// @brief Static assertion to ensure the template parameter T is a floating-point type.
  static_assert(std::is_floating_point<T>::value, "numc::tensor mozna uzywac TYLKO z typami zmiennoprzecinkowymi!");

 private:
  /// @brief The dimensions of the tensor (e.g., {rows, cols} for a 2D matrix).
  std::vector<size_t> _shape;

  /// @brief The memory strides for each dimension, used to compute flat indices.
  std::vector<size_t> _strides;

  /// @brief The flat 1D array storing the actual tensor data.
  std::vector<T> _data;

  /// @brief Computes the memory strides based on the current shape.
  void compute_strides() {
    _strides.resize(_shape.size());
    if (_shape.empty()) return;
    size_t current_stride = 1;
    for (int i = static_cast<int>(_shape.size()) - 1; i >= 0; --i) {
      _strides[i] = current_stride;
      current_stride *= _shape[i];
    }
  }

  /// @brief Converts multi-dimensional indices into a 1D flat index for memory access.
  /// @param indices The multi-dimensional coordinates.
  /// @return The computed flat index.
  /// @throws std::invalid_argument If the number of indices does not match the tensor dimensions.
  /// @throws std::out_of_range If any index exceeds its corresponding dimension bounds.
  size_t get_flat_index(const std::vector<size_t>& indices) const {
    if (indices.size() != _shape.size()) {
      Log::Error("Wymiary nie pasuja!");
      throw std::invalid_argument("Incorrect number of dimensions accessed.");
    }
    size_t flat_idx = 0;
    for (size_t i = 0; i < indices.size(); ++i) {
      if (indices[i] >= _shape[i]) throw std::out_of_range("Tensor index out of bounds.");
      flat_idx += indices[i] * _strides[i];
    }
    return flat_idx;
  }

 public:
  /// @brief Default constructor. Creates an empty tensor.
  tensor() = default;

  /// @brief Constructs a tensor with a specific shape, initialized to zeros.
  /// @param shape A vector specifying the size of each dimension.
  explicit tensor(std::vector<size_t> shape) : _shape(shape) {
    size_t total_size = 1;
    for (size_t dim : _shape) total_size *= dim;
    _data.resize(total_size, T(0.0));
    compute_strides();
  }

  /// @brief Constructs a tensor with a specific shape and fills it with provided flat data.
  /// @param shape A vector specifying the size of each dimension.
  /// @param flat_data A flat 1D vector containing the data elements.
  /// @throws std::invalid_argument If the size of flat_data does not match the total size defined by the shape.
  tensor(std::vector<size_t> shape, std::vector<T> flat_data) : _shape(shape), _data(std::move(flat_data)) {
    size_t total_size = 1;
    for (size_t dim : _shape) total_size *= dim;
    if (_data.size() != total_size) throw std::invalid_argument("Data size does not match shape total size.");
    compute_strides();
  }

  /// @brief Gets the shape (dimensions) of the tensor.
  /// @return A constant reference to the shape vector.
  inline const std::vector<size_t>& shape() const { return _shape; }

  /// @brief Gets the total number of elements in the tensor.
  /// @return The total element count.
  inline size_t total_size() const { return _data.size(); }

  /// @brief Gets the number of dimensions (rank) of the tensor.
  /// @return The dimension count (e.g., 2 for a matrix).
  inline size_t dimensions() const { return _shape.size(); }

  /// @brief Creates an identity matrix (2D tensor) of size n x n.
  /// @param n The number of rows and columns.
  /// @return An identity matrix tensor.
  static tensor eye(size_t n) {
    tensor res({n, n});
    for (size_t i = 0; i < n; ++i) res(i, i) = T(1.0);
    return res;
  }

  /// @brief Accesses an element in the tensor using multi-dimensional coordinates.
  /// @tparam Args Integer types representing the coordinates.
  /// @param args The coordinates of the desired element.
  /// @return A reference to the element.
  template <typename... Args>
  T& operator()(Args... args) {
    std::vector<size_t> indices = {static_cast<size_t>(args)...};
    return _data[get_flat_index(indices)];
  }

  /// @brief Accesses an element in the tensor using multi-dimensional coordinates (const version).
  /// @tparam Args Integer types representing the coordinates.
  /// @param args The coordinates of the desired element.
  /// @return A constant reference to the element.
  template <typename... Args>
  const T& operator()(Args... args) const {
    std::vector<size_t> indices = {static_cast<size_t>(args)...};
    return _data[get_flat_index(indices)];
  }

  /// @brief Compound element-wise addition operator.
  /// @param other The tensor to add. Must have the exact same shape.
  /// @return A reference to the current tensor.
  /// @throws std::invalid_argument If shapes do not match.
  tensor& operator+=(const tensor& other) {
    if (_shape != other._shape) throw std::invalid_argument("Shape mismatch for addition.");
    for (size_t i = 0; i < _data.size(); ++i) _data[i] += other._data[i];
    return *this;
  }

  /// @brief Compound element-wise subtraction operator.
  /// @param other The tensor to subtract. Must have the exact same shape.
  /// @return A reference to the current tensor.
  /// @throws std::invalid_argument If shapes do not match.
  tensor& operator-=(const tensor& other) {
    if (_shape != other._shape) throw std::invalid_argument("Shape mismatch for subtraction.");
    for (size_t i = 0; i < _data.size(); ++i) _data[i] -= other._data[i];
    return *this;
  }

  /// @brief Computes the Hadamard (element-wise) product of two tensors.
  /// @param other The tensor to multiply with. Must have the exact same shape.
  /// @return A new tensor representing the element-wise product.
  /// @throws std::invalid_argument If shapes do not match.
  tensor hadamard(const tensor& other) const {
    if (_shape != other._shape) throw std::invalid_argument("Shape mismatch for Hadamard product.");
    tensor result(_shape);
    for (size_t i = 0; i < _data.size(); ++i) result._data[i] = _data[i] * other._data[i];
    return result;
  }

  /// @brief Performs classic 2D matrix multiplication (dot product).
  /// @param other The right-hand side 2D tensor.
  /// @return A new 2D tensor representing the matrix product.
  /// @throws std::invalid_argument If tensors are not 2D or inner dimensions do not match.
  tensor matmul(const tensor& other) const {
    if (this->dimensions() != 2 || other.dimensions() != 2) {
      Log::Error("Obecnie mnozenie (*) obsluguje tylko macierze 2D!");
      throw std::invalid_argument("Matmul requires 2D tensors.");
    }
    if (_shape[1] != other._shape[0]) {
      Log::Error("Wymiary wewnetrzne macierzy nie pasuja do mnozenia!");
      throw std::invalid_argument("Inner dimensions must match for matrix multiplication.");
    }

    size_t rows = _shape[0];
    size_t common_dim = _shape[1];
    size_t cols = other._shape[1];

    tensor result({rows, cols});

    for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < cols; ++j) {
        T sum = T(0.0);
        for (size_t k = 0; k < common_dim; ++k) {
          sum += (*this)(i, k) * other(k, j);
        }
        result(i, j) = sum;
      }
    }
    return result;
  }

  /// @brief Computes the mathematical inverse of a square 2D matrix using Gauss-Jordan elimination.
  /// @return A new tensor representing the inverted matrix.
  /// @throws std::invalid_argument If the tensor is not 2D or not square.
  /// @throws std::runtime_error If the matrix is singular (determinant = 0) and cannot be inverted.
  tensor inverse() const {
    if (dimensions() != 2 || _shape[0] != _shape[1]) {
      Log::Error("Tylko kwadratowe macierze 2D moga byc odwrocone!");
      throw std::invalid_argument("Inverse requires a square 2D matrix.");
    }

    size_t n = _shape[0];
    tensor result = tensor::eye(n);
    tensor tmp = *this;

    for (size_t i = 0; i < n; ++i) {
      size_t pivot = i;
      for (size_t j = i + 1; j < n; ++j) {
        if (std::abs(tmp(j, i)) > std::abs(tmp(pivot, i))) pivot = j;
      }

      if (tmp(pivot, i) == T(0.0)) {
        Log::Error("Macierz jest osobliwa (wyznacznik = 0). Nie mozna odwrocic!");
        throw std::runtime_error("Singular matrix cannot be inverted.");
      }

      if (pivot != i) {
        for (size_t k = 0; k < n; ++k) {
          std::swap(tmp(i, k), tmp(pivot, k));
          std::swap(result(i, k), result(pivot, k));
        }
      }

      T pivot_val = tmp(i, i);
      for (size_t k = 0; k < n; ++k) {
        tmp(i, k) /= pivot_val;
        result(i, k) /= pivot_val;
      }

      for (size_t j = 0; j < n; ++j) {
        if (j != i) {
          T factor = tmp(j, i);
          for (size_t k = 0; k < n; ++k) {
            tmp(j, k) -= factor * tmp(i, k);
            result(j, k) -= factor * result(i, k);
          }
        }
      }
    }

    return result;
  }

  /// @brief Compound scalar multiplication operator.
  /// @param scalar The scalar value to multiply each element by.
  /// @return A reference to the current tensor.
  tensor& operator*=(T scalar) {
    for (T& val : _data) val *= scalar;
    return *this;
  }

  /// @brief Compound scalar division operator.
  /// @param scalar The scalar value to divide each element by.
  /// @return A reference to the current tensor.
  /// @throws std::domain_error If the scalar is exactly zero.
  tensor& operator/=(T scalar) {
    if (scalar == T(0.0)) throw std::domain_error("Division by zero scalar.");
    for (T& val : _data) val /= scalar;
    return *this;
  }

  /// @brief Element-wise addition of two tensors.
  friend tensor operator+(tensor lhs, const tensor& rhs) { return lhs += rhs; }

  /// @brief Element-wise subtraction of two tensors.
  friend tensor operator-(tensor lhs, const tensor& rhs) { return lhs -= rhs; }

  /// @brief Matrix multiplication of two 2D tensors. Note: This is NOT element-wise (Hadamard) multiplication.
  friend tensor operator*(const tensor& lhs, const tensor& rhs) { return lhs.matmul(rhs); }

  /// @brief Matrix "division" (lhs multiplied by the inverse of rhs). Note: rhs must be a square invertible matrix.
  friend tensor operator/(const tensor& lhs, const tensor& rhs) { return lhs.matmul(rhs.inverse()); }

  /// @brief Multiplication of a tensor by a scalar.
  friend tensor operator*(tensor t, T scalar) { return t *= scalar; }

  /// @brief Multiplication of a scalar by a tensor.
  friend tensor operator*(T scalar, tensor t) { return t *= scalar; }

  /// @brief Division of a tensor by a scalar.
  friend tensor operator/(tensor t, T scalar) { return t /= scalar; }

  /// @brief Division of a scalar by a tensor (returns scalar * inverse of t).
  friend tensor operator/(T scalar, const tensor& t) {
    tensor inv = t.inverse();
    inv *= scalar;
    return inv;
  }
};

/// @brief Output stream operator for tensors. Neatly formats 1D and 2D tensors, and provides summary info for higher dimensions.
/// @tparam T The underlying floating-point type.
/// @param os The output stream.
/// @param t The tensor to be printed.
/// @return The output stream.
template <typename T>
std::ostream& operator<<(std::ostream& os, const tensor<T>& t) {
  if (t.dimensions() == 1) {
    os << "[";
    for (size_t i = 0; i < t.shape()[0]; ++i) {
      os << t(i) << (i + 1 < t.shape()[0] ? ", " : "");
    }
    os << "]";
  } else if (t.dimensions() == 2) {
    os << "[\n";
    for (size_t i = 0; i < t.shape()[0]; ++i) {
      os << "  [";
      for (size_t j = 0; j < t.shape()[1]; ++j) {
        os << t(i, j) << (j + 1 < t.shape()[1] ? ", " : "");
      }
      os << "]" << (i + 1 < t.shape()[0] ? ",\n" : "\n");
    }
    os << "]";
  } else {
    os << "[Tensor " << t.dimensions() << "D: total elements = " << t.total_size() << "]";
  }
  return os;
}

/// @} // Koniec grupy common

}  // namespace numc