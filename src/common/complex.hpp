#pragma once

#include "../inc.hpp"
#include "../utility/log.hpp"
#include "./constants.hpp"

namespace numc {

/// @brief Class representing complex numbers with a template parameter for the underlying floating-point type.
/// @tparam T The underlying floating-point type (e.g., float, double, long double).
template <typename T = double>
class complex {
  /// @brief Static assertion to ensure that the template parameter T is a floating-point type.
  static_assert(std::is_floating_point<T>::value, "numc::complex mozna uzywac TYLKO z typami zmiennoprzecinkowymi (float, double, long double)!");

  /// @brief The real part of the complex number.
  T x;

  /// @brief The imaginary part of the complex number.
  T y;

 public:
  /// @brief Constructor for complex numbers.
  /// @param real The real part of the complex number. Default is 0.0.
  /// @param imag The imaginary part of the complex number. Default is 0.0.
  complex(T real = T(0.0), T imag = T(0.0)) : x(real), y(imag) {}

  /// @brief Copy constructor for complex numbers.
  /// @param other The complex number to copy.
  complex(const complex& other) : x(other.x), y(other.y) {}

  /// @brief Default destructor.
  ~complex() = default;

  /// @brief Static method to create a complex number representing the imaginary unit (0 + 1i).
  /// @return An imaginary unit complex number.
  static inline complex i() { return complex(T(0.0), T(1.0)); }

  /// @brief Get the real part of the complex number.
  /// @return The real part.
  inline T real() const { return x; }

  /// @brief Get the imaginary part of the complex number.
  /// @return The imaginary part.
  inline T imag() const { return y; }

  /// @brief Calculate the magnitude (absolute value) of the complex number.
  /// @return The magnitude of the complex number.
  T abs() const;

  /// @brief Calculate the argument (angle, phase) of the complex number.
  /// @return The argument of the complex number in radians.
  T arg() const;

  /// @brief Calculate the complex conjugate of the complex number.
  /// @return The complex conjugate.
  complex conjugate() const;

  /// @brief Raise the complex number to the power of another complex number.
  /// @param other The exponent complex number.
  /// @return The result of exponentiation.
  complex pow(const complex& other) const;

  /// @brief Calculate the nth roots of the complex number.
  /// @param n The degree of the root to calculate. Must be a positive integer.
  /// @return A vector containing the nth roots of the complex number.
  std::vector<complex> roots(int n) const;

  /// @brief Round the real and imaginary parts of the complex number to a specified number of decimal places.
  /// @param places The number of decimal places to round to. Default is 0 (round to the nearest integer).
  /// @return A new complex number with rounded components.
  complex round(int places = 0) const;

  /// @brief Compound addition assignment operator.
  /// @param other The complex number to add.
  /// @return A reference to the current complex number.
  complex& operator+=(const complex& other);

  /// @brief Compound subtraction assignment operator.
  /// @param other The complex number to subtract.
  /// @return A reference to the current complex number.
  complex& operator-=(const complex& other);

  /// @brief Compound multiplication assignment operator.
  /// @param other The complex number to multiply by.
  /// @return A reference to the current complex number.
  complex& operator*=(const complex& other);

  /// @brief Compound division assignment operator.
  /// @param other The complex number to divide by.
  /// @return A reference to the current complex number.
  complex& operator/=(const complex& other);

  /// @brief Compound exponentiation assignment operator.
  /// @param other The exponent complex number.
  /// @return A reference to the current complex number.
  complex& operator^=(const complex& other);

  /// @brief Equality operator. Compares exact values of real and imaginary parts.
  /// @param other The complex number to compare with.
  /// @return True if strictly equal, false otherwise.
  bool operator==(const complex& other) const;

  /// @brief Inequality operator.
  /// @param other The complex number to compare with.
  /// @return True if not equal, false otherwise.
  bool operator!=(const complex& other) const;

  /// @brief Less-than operator. Compares the magnitudes (abs) of the complex numbers. Use with caution.
  /// @param other The complex number to compare with.
  /// @return True if the magnitude of this is less than the magnitude of other.
  bool operator<(const complex& other) const;

  /// @brief Greater-than operator. Compares the magnitudes (abs) of the complex numbers. Use with caution.
  /// @param other The complex number to compare with.
  /// @return True if the magnitude of this is greater than the magnitude of other.
  bool operator>(const complex& other) const;

  /// @brief Less-than-or-equal operator. Compares the magnitudes of the complex numbers. Use with caution.
  /// @param other The complex number to compare with.
  /// @return True if the magnitude of this is less than or equal to the magnitude of other.
  bool operator<=(const complex& other) const;

  /// @brief Greater-than-or-equal operator. Compares the magnitudes of the complex numbers. Use with caution.
  /// @param other The complex number to compare with.
  /// @return True if the magnitude of this is greater than or equal to the magnitude of other.
  bool operator>=(const complex& other) const;

  /// @brief Unary bitwise NOT operator, overloaded to return the complex conjugate.
  /// @return The complex conjugate.
  complex operator~() const;

  /// @brief Unary negation operator. Negates both real and imaginary parts.
  /// @return The negated complex number.
  complex operator-() const;

  /// @brief Binary addition operator.
  friend complex operator+(complex lhs, const complex& rhs) { return lhs += rhs; }

  /// @brief Binary subtraction operator.
  friend complex operator-(complex lhs, const complex& rhs) { return lhs -= rhs; }

  /// @brief Binary multiplication operator.
  friend complex operator*(complex lhs, const complex& rhs) { return lhs *= rhs; }

  /// @brief Binary division operator.
  friend complex operator/(complex lhs, const complex& rhs) { return lhs /= rhs; }

  /// @brief Binary exponentiation operator.
  friend complex operator^(complex lhs, const complex& rhs) { return lhs ^= rhs; }

  /// @brief Output stream operator. Formats the complex number in a human-readable form (e.g., 3 + 4i, -2i, 5).
  /// @tparam U The type of the complex number.
  /// @param os The output stream.
  /// @param c The complex number to output.
  /// @return The output stream.
  template <typename U>
  friend std::ostream& operator<<(std::ostream& os, const complex<U>& c);
};

template <typename T>
inline T Re(const complex<T>& c) {
  return c.real();
}
template <typename T>
inline T Im(const complex<T>& c) {
  return c.imag();
}
template <typename T>
inline complex<T> pow(const complex<T>& base, const complex<T>& exp) {
  return base.pow(exp);
}

template <typename T>
T complex<T>::abs() const {
  return std::hypot(x, y);
}

template <typename T>
T complex<T>::arg() const {
  return std::atan2(y, x);
}

template <typename T>
complex<T> complex<T>::round(int places) const {
  T factor = std::pow(T(10.0), T(places));
  return complex<T>(std::round(x * factor) / factor, std::round(y * factor) / factor);
}

template <typename T>
complex<T>& complex<T>::operator+=(const complex<T>& other) {
  x += other.x;
  y += other.y;
  return *this;
}

template <typename T>
complex<T>& complex<T>::operator-=(const complex<T>& other) {
  x -= other.x;
  y -= other.y;
  return *this;
}

template <typename T>
complex<T>& complex<T>::operator*=(const complex<T>& other) {
  T new_x = x * other.x - y * other.y;
  T new_y = x * other.y + y * other.x;
  x = new_x;
  y = new_y;
  return *this;
}

template <typename T>
complex<T>& complex<T>::operator/=(const complex<T>& other) {
  T denom = other.x * other.x + other.y * other.y;
  T new_x = (x * other.x + y * other.y) / denom;
  T new_y = (y * other.x - x * other.y) / denom;
  x = new_x;
  y = new_y;
  return *this;
}

template <typename T>
complex<T> complex<T>::operator-() const {
  return complex<T>(-x, -y);
}

template <typename T>
complex<T> complex<T>::operator~() const {
  return complex<T>(x, -y);
}

template <typename T>
complex<T> complex<T>::conjugate() const {
  return complex<T>(x, -y);
}

template <typename T>
complex<T> complex<T>::pow(const complex<T>& other) const {
  T r = this->abs();
  T theta = this->arg();

  if (r == T(0.0)) return complex<T>(T(0.0), T(0.0));

  T new_r = std::pow(r, other.x) * std::exp(-other.y * theta);
  T new_theta = other.x * theta + other.y * std::log(r);

  return complex<T>(new_r * std::cos(new_theta), new_r * std::sin(new_theta));
}

template <typename T>
complex<T>& complex<T>::operator^=(const complex<T>& other) {
  *this = this->pow(other);
  return *this;
}

template <typename T>
std::vector<complex<T>> complex<T>::roots(int n) const {
  if (n <= 0) {
    Log::Error("Number of roots must be positive.");
    return {};
  }

  std::vector<complex<T>> result(n);

  T r = std::pow(this->abs(), T(1.0) / T(n));
  T theta = this->arg();

  for (int k = 0; k < n; ++k) {
    T root_theta = (theta + T(2.0) * T(numc::PI) * T(k)) / T(n);
    result[k] = complex<T>(r * std::cos(root_theta), r * std::sin(root_theta));
  }

  return result;
}

template <typename T>
bool complex<T>::operator==(const complex<T>& other) const {
  return x == other.x && y == other.y;
}

template <typename T>
bool complex<T>::operator!=(const complex<T>& other) const {
  return !(*this == other);
}

template <typename T>
bool complex<T>::operator<(const complex<T>& other) const {
  Log::Warn(
      "Comparison of complex numbers is not well-defined. Comparing by magnitude. If you want to compare by real or imaginary part, use Re() or Im() "
      "methods.");
  return this->abs() < other.abs();
}

template <typename T>
bool complex<T>::operator>(const complex<T>& other) const {
  Log::Warn(
      "Comparison of complex numbers is not well-defined. Comparing by magnitude. If you want to compare by real or imaginary part, use real() or "
      "imag() methods or Re() or Im() functions.");
  return this->abs() > other.abs();
}

template <typename T>
bool complex<T>::operator<=(const complex<T>& other) const {
  if (*this == other) return true;
  return *this < other;
}

template <typename T>
bool complex<T>::operator>=(const complex<T>& other) const {
  if (*this == other) return true;
  return *this > other;
}

template <typename T>
std::ostream& operator<<(std::ostream& os, const complex<T>& c) {
  if (c.real() == T(0.0) && c.imag() == T(0.0)) {
    return os << "0";
  }

  if (c.real() != T(0.0)) {
    os << c.real();
    if (c.imag() > T(0.0))
      os << " + ";
    else if (c.imag() < T(0.0))
      os << " - ";
  } else {
    if (c.imag() < T(0.0)) os << "-";
  }

  if (c.imag() != T(0.0)) {
    T abs_y = std::abs(c.imag());
    if (abs_y != T(1.0)) {
      os << abs_y;
    }
    os << "i";
  }

  return os;
}

}  // namespace numc