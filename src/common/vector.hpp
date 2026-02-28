#pragma once

#include "../inc.hpp"
#include "../utility/log.hpp"  // Twój logger

namespace numc {

template <typename T = double>
class vector {
  static_assert(std::is_floating_point<T>::value, "numc::vector mozna uzywac TYLKO z typami zmiennoprzecinkowymi (float, double, long double)!");

  std::vector<T> data;

 public:
  vector() = default;
  explicit vector(size_t size) : data(size, T(0.0)) {}
  vector(std::initializer_list<T> values) : data(values) {}

  ~vector() = default;

  inline size_t size() const { return data.size(); }

  T& operator[](size_t index) {
    if (index >= data.size()) {
      Log::Error("Index out of bounds in vector::operator[]");
      throw std::out_of_range("Index out of bounds");
    }
    return data[index];
  }

  const T& operator[](size_t index) const {
    if (index >= data.size()) {
      Log::Error("Index out of bounds in vector::operator[]");
      throw std::out_of_range("Index out of bounds");
    }
    return data[index];
  }

  T& operator()(size_t index) { return (*this)[index]; }

  const T& operator()(size_t index) const { return (*this)[index]; }

  T dot(const vector& other) const {
    if (this->size() != other.size()) {
      Log::Error("Dot product failed: vectors have different sizes (" + std::to_string(this->size()) + " vs " + std::to_string(other.size()) + ").");
      throw std::invalid_argument("Vectors must be of the same size for dot product.");
    }
    T result = T(0.0);
    for (size_t i = 0; i < this->size(); ++i) {
      result += data[i] * other.data[i];
    }
    return result;
  }

  vector cross(const vector& other) const {
    if (this->size() != 3 || other.size() != 3) {
      Log::Error("Cross product failed: both vectors must be 3D. (Got " + std::to_string(this->size()) + "D and " + std::to_string(other.size()) +
                 "D).");
      throw std::invalid_argument("Cross product requires two vectors of size 3.");
    }
    return {data[1] * other.data[2] - data[2] * other.data[1],
            data[2] * other.data[0] - data[0] * other.data[2],
            data[0] * other.data[1] - data[1] * other.data[0]};
  }

  T norm() const {
    T sum = T(0.0);
    for (T val : data) sum += val * val;
    return std::sqrt(sum);
  }

  vector normalize() const {
    T n = this->norm();
    if (n == T(0.0)) {
      Log::Error("Normalization failed: cannot normalize a zero vector.");
      throw std::domain_error("Cannot normalize a zero vector.");
    }
    return *this / n;
  }

  vector& operator+=(const vector& other) {
    if (this->size() != other.size()) {
      Log::Error("Addition failed: vectors have different sizes.");
      throw std::invalid_argument("Vector sizes must match for addition.");
    }
    for (size_t i = 0; i < this->size(); ++i) data[i] += other.data[i];
    return *this;
  }

  vector& operator-=(const vector& other) {
    if (this->size() != other.size()) {
      Log::Error("Subtraction failed: vectors have different sizes.");
      throw std::invalid_argument("Vector sizes must match for subtraction.");
    }
    for (size_t i = 0; i < this->size(); ++i) data[i] -= other.data[i];
    return *this;
  }

  vector& operator^=(const vector& other) {
    *this = this->cross(other);
    return *this;
  }

  vector& operator*=(T scalar) {
    for (T& val : data) val *= scalar;
    return *this;
  }

  vector& operator/=(T scalar) {
    if (scalar == T(0.0)) {
      Log::Error("Division by zero failed: scalar is 0.0.");
      throw std::domain_error("Division by zero scalar.");
    }
    for (T& val : data) val /= scalar;
    return *this;
  }

  friend vector operator+(vector lhs, const vector& rhs) { return lhs += rhs; }
  friend vector operator-(vector lhs, const vector& rhs) { return lhs -= rhs; }

  friend vector operator*(vector v, T scalar) { return v *= scalar; }
  friend vector operator*(T scalar, vector v) { return v *= scalar; }
  friend vector operator/(vector v, T scalar) { return v /= scalar; }

  friend T operator*(const vector& lhs, const vector& rhs) { return lhs.dot(rhs); }
  friend vector operator^(const vector& lhs, const vector& rhs) { return lhs.cross(rhs); }

  template <typename U>
  friend std::ostream& operator<<(std::ostream& os, const vector<U>& v);
};

template <typename T>
std::ostream& operator<<(std::ostream& os, const vector<T>& v) {
  os << "(";
  for (size_t i = 0; i < v.size(); ++i) {
    os << v[i];
    if (i + 1 < v.size()) os << ", ";
  }
  os << ")";
  return os;
}

}  // namespace numc