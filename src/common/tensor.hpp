#pragma once

#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {

template <typename T = double>
class tensor {
  static_assert(std::is_floating_point<T>::value, "numc::tensor mozna uzywac TYLKO z typami zmiennoprzecinkowymi!");

 private:
  std::vector<size_t> _shape;
  std::vector<size_t> _strides;
  std::vector<T> _data;

  void compute_strides() {
    _strides.resize(_shape.size());
    if (_shape.empty()) return;
    size_t current_stride = 1;
    for (int i = static_cast<int>(_shape.size()) - 1; i >= 0; --i) {
      _strides[i] = current_stride;
      current_stride *= _shape[i];
    }
  }

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
  tensor() = default;

  explicit tensor(std::vector<size_t> shape) : _shape(shape) {
    size_t total_size = 1;
    for (size_t dim : _shape) total_size *= dim;
    _data.resize(total_size, T(0.0));
    compute_strides();
  }

  tensor(std::vector<size_t> shape, std::vector<T> flat_data) : _shape(shape), _data(std::move(flat_data)) {
    size_t total_size = 1;
    for (size_t dim : _shape) total_size *= dim;
    if (_data.size() != total_size) throw std::invalid_argument("Data size does not match shape total size.");
    compute_strides();
  }

  inline const std::vector<size_t>& shape() const { return _shape; }
  inline size_t total_size() const { return _data.size(); }
  inline size_t dimensions() const { return _shape.size(); }

  static tensor eye(size_t n) {
    tensor res({n, n});
    for (size_t i = 0; i < n; ++i) res(i, i) = T(1.0);
    return res;
  }

  template <typename... Args>
  T& operator()(Args... args) {
    std::vector<size_t> indices = {static_cast<size_t>(args)...};
    return _data[get_flat_index(indices)];
  }

  template <typename... Args>
  const T& operator()(Args... args) const {
    std::vector<size_t> indices = {static_cast<size_t>(args)...};
    return _data[get_flat_index(indices)];
  }

  tensor& operator+=(const tensor& other) {
    if (_shape != other._shape) throw std::invalid_argument("Shape mismatch for addition.");
    for (size_t i = 0; i < _data.size(); ++i) _data[i] += other._data[i];
    return *this;
  }

  tensor& operator-=(const tensor& other) {
    if (_shape != other._shape) throw std::invalid_argument("Shape mismatch for subtraction.");
    for (size_t i = 0; i < _data.size(); ++i) _data[i] -= other._data[i];
    return *this;
  }

  tensor hadamard(const tensor& other) const {
    if (_shape != other._shape) throw std::invalid_argument("Shape mismatch for Hadamard product.");
    tensor result(_shape);
    for (size_t i = 0; i < _data.size(); ++i) result._data[i] = _data[i] * other._data[i];
    return result;
  }

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
          // (*this)(i, k) razy other(k, j)
          sum += (*this)(i, k) * other(k, j);
        }
        result(i, j) = sum;
      }
    }
    return result;
  }

  tensor inverse() const {
    if (dimensions() != 2 || _shape[0] != _shape[1]) {
      Log::Error("Tylko kwadratowe macierze 2D moga byc odwrocone!");
      throw std::invalid_argument("Inverse requires a square 2D matrix.");
    }

    size_t n = _shape[0];
    tensor result = tensor::eye(n);  // Macierz jednostkowa
    tensor tmp = *this;              // Kopia robocza (augmented matrix)

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

  tensor& operator*=(T scalar) {
    for (T& val : _data) val *= scalar;
    return *this;
  }

  tensor& operator/=(T scalar) {
    if (scalar == T(0.0)) throw std::domain_error("Division by zero scalar.");
    for (T& val : _data) val /= scalar;
    return *this;
  }

  friend tensor operator+(tensor lhs, const tensor& rhs) { return lhs += rhs; }
  friend tensor operator-(tensor lhs, const tensor& rhs) { return lhs -= rhs; }

  friend tensor operator*(const tensor& lhs, const tensor& rhs) { return lhs.matmul(rhs); }

  friend tensor operator/(const tensor& lhs, const tensor& rhs) { return lhs.matmul(rhs.inverse()); }

  friend tensor operator*(tensor t, T scalar) { return t *= scalar; }
  friend tensor operator*(T scalar, tensor t) { return t *= scalar; }
  friend tensor operator/(tensor t, T scalar) { return t /= scalar; }

  friend tensor operator/(T scalar, const tensor& t) {
    tensor inv = t.inverse();
    inv *= scalar;
    return inv;
  }
};

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

}  // namespace numc