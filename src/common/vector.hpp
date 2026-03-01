#pragma once

#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {

/// @addtogroup structures
/// @{

/// @brief Class representing a mathematical vector (1D array) for linear algebra and physics calculations.
/// Supports vector arithmetic, dot products, and 3D cross products.
/// @tparam T The underlying floating-point type (e.g., float, double, long double).
template <typename T = double>
class vector {
  /// @brief Static assertion to ensure the template parameter T is a floating-point type.
  static_assert(std::is_floating_point<T>::value, "numc::vector mozna uzywac TYLKO z typami zmiennoprzecinkowymi (float, double, long double)!");

 private:
  /// @brief The internal storage for the vector's components.
  std::vector<T> data;

 public:
  /// @brief Default constructor. Creates an empty vector.
  vector() = default;

  /// @brief Constructs a zero-initialized vector of a specific size.
  /// @param size The number of elements in the vector.
  explicit vector(size_t size) : data(size, T(0.0)) {}

  /// @brief Constructs a vector from an initializer list (e.g., vector v = {1.0, 2.0, 3.0}).
  /// @param values The initializer list of values.
  vector(std::initializer_list<T> values) : data(values) {}

  /// @brief Default destructor.
  ~vector() = default;

  /// @brief Gets the number of elements (dimensionality) of the vector.
  /// @return The size of the vector.
  inline size_t size() const { return data.size(); }

  /// @brief Array subscript operator for element access.
  /// @param index The zero-based index of the element.
  /// @return A reference to the element at the specified index.
  /// @throws std::out_of_range If the index is out of bounds.
  T& operator[](size_t index) {
    if (index >= data.size()) {
      Log::Error("Index out of bounds in vector::operator[]");
      throw std::out_of_range("Index out of bounds");
    }
    return data[index];
  }

  /// @brief Array subscript operator for element access (const version).
  /// @param index The zero-based index of the element.
  /// @return A constant reference to the element at the specified index.
  /// @throws std::out_of_range If the index is out of bounds.
  const T& operator[](size_t index) const {
    if (index >= data.size()) {
      Log::Error("Index out of bounds in vector::operator[]");
      throw std::out_of_range("Index out of bounds");
    }
    return data[index];
  }

  /// @brief Function call operator for element access (alias for operator[]).
  /// @param index The zero-based index of the element.
  /// @return A reference to the element.
  T& operator()(size_t index) { return (*this)[index]; }

  /// @brief Function call operator for element access (const version, alias for operator[]).
  /// @param index The zero-based index of the element.
  /// @return A constant reference to the element.
  const T& operator()(size_t index) const { return (*this)[index]; }

  /// @brief Computes the scalar dot product of this vector with another vector.
  /// @param other The second vector.
  /// @return The scalar result of the dot product.
  /// @throws std::invalid_argument If the vectors have different sizes.
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

  /// @brief Computes the cross product of this vector with another vector. Strictly for 3D vectors.
  /// @param other The second 3D vector.
  /// @return A new 3D vector orthogonal to both input vectors.
  /// @throws std::invalid_argument If either vector is not strictly 3-dimensional.
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

  /// @brief Computes the Euclidean norm (magnitude/length) of the vector.
  /// @return The length of the vector.
  T norm() const {
    T sum = T(0.0);
    for (T val : data) sum += val * val;
    return std::sqrt(sum);
  }

  /// @brief Returns a normalized copy of the vector (unit vector with length 1.0).
  /// @return A new normalized vector.
  /// @throws std::domain_error If the vector is a zero vector (magnitude is 0).
  vector normalize() const {
    T n = this->norm();
    if (n == T(0.0)) {
      Log::Error("Normalization failed: cannot normalize a zero vector.");
      throw std::domain_error("Cannot normalize a zero vector.");
    }
    return *this / n;
  }

  /// @brief Compound element-wise addition operator.
  /// @param other The vector to add.
  /// @return A reference to the current vector.
  /// @throws std::invalid_argument If vector sizes mismatch.
  vector& operator+=(const vector& other) {
    if (this->size() != other.size()) {
      Log::Error("Addition failed: vectors have different sizes.");
      throw std::invalid_argument("Vector sizes must match for addition.");
    }
    for (size_t i = 0; i < this->size(); ++i) data[i] += other.data[i];
    return *this;
  }

  /// @brief Compound element-wise subtraction operator.
  /// @param other The vector to subtract.
  /// @return A reference to the current vector.
  /// @throws std::invalid_argument If vector sizes mismatch.
  vector& operator-=(const vector& other) {
    if (this->size() != other.size()) {
      Log::Error("Subtraction failed: vectors have different sizes.");
      throw std::invalid_argument("Vector sizes must match for subtraction.");
    }
    for (size_t i = 0; i < this->size(); ++i) data[i] -= other.data[i];
    return *this;
  }

  /// @brief Compound cross product operator. Replaces current vector with the cross product result.
  /// @param other The vector to cross with.
  /// @return A reference to the current vector.
  vector& operator^=(const vector& other) {
    *this = this->cross(other);
    return *this;
  }

  /// @brief Compound scalar multiplication operator.
  /// @param scalar The scalar to multiply by.
  /// @return A reference to the current vector.
  vector& operator*=(T scalar) {
    for (T& val : data) val *= scalar;
    return *this;
  }

  /// @brief Compound scalar division operator.
  /// @param scalar The scalar to divide by.
  /// @return A reference to the current vector.
  /// @throws std::domain_error If scalar is zero.
  vector& operator/=(T scalar) {
    if (scalar == T(0.0)) {
      Log::Error("Division by zero failed: scalar is 0.0.");
      throw std::domain_error("Division by zero scalar.");
    }
    for (T& val : data) val /= scalar;
    return *this;
  }

  /// @brief Element-wise vector addition.
  friend vector operator+(vector lhs, const vector& rhs) { return lhs += rhs; }

  /// @brief Element-wise vector subtraction.
  friend vector operator-(vector lhs, const vector& rhs) { return lhs -= rhs; }

  /// @brief Scalar multiplication (vector * scalar).
  friend vector operator*(vector v, T scalar) { return v *= scalar; }

  /// @brief Scalar multiplication (scalar * vector).
  friend vector operator*(T scalar, vector v) { return v *= scalar; }

  /// @brief Scalar division (vector / scalar).
  friend vector operator/(vector v, T scalar) { return v /= scalar; }

  /// @brief Vector dot product. Returns a scalar.
  friend T operator*(const vector& lhs, const vector& rhs) { return lhs.dot(rhs); }

  /// @brief Vector cross product. Returns an orthogonal vector.
  friend vector operator^(const vector& lhs, const vector& rhs) { return lhs.cross(rhs); }

  /// @brief Output stream operator. Formats the vector as a tuple, e.g., (1, 2, 3).
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

/// @} // Koniec grupy common

}  // namespace numc