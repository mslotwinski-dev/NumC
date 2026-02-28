#pragma once

#include "../inc.hpp"
#include "../utility/log.hpp"

// // Zabezpieczenie, gdyby kompilator nie definiował M_PI (standard C++ formalnie tego nie wymaga)
// #ifndef M_PI
// #define M_PI 3.14159265358979323846
// #endif

namespace numc {

template <typename T = double>
class complex {
  static_assert(std::is_floating_point<T>::value, "numc::complex mozna uzywac TYLKO z typami zmiennoprzecinkowymi (float, double, long double)!");

  T x;
  T y;

 public:
  complex(T real = T(0.0), T imag = T(0.0)) : x(real), y(imag) {}
  ~complex() = default;

  static inline complex i() { return complex(T(0.0), T(1.0)); }

  inline T real() const { return x; }
  inline T imag() const { return y; }

  T abs() const;
  T arg() const;

  complex conjugate() const;
  complex pow(const complex& other) const;
  std::vector<complex> roots(int n) const;
  complex round(int places = 0) const;

  complex& operator+=(const complex& other);
  complex& operator-=(const complex& other);
  complex& operator*=(const complex& other);
  complex& operator/=(const complex& other);
  complex& operator^=(const complex& other);

  bool operator==(const complex& other) const;
  bool operator!=(const complex& other) const;
  bool operator<(const complex& other) const;
  bool operator>(const complex& other) const;
  bool operator<=(const complex& other) const;
  bool operator>=(const complex& other) const;

  complex operator~() const;
  complex operator-() const;

  friend complex operator+(complex lhs, const complex& rhs) { return lhs += rhs; }
  friend complex operator-(complex lhs, const complex& rhs) { return lhs -= rhs; }
  friend complex operator*(complex lhs, const complex& rhs) { return lhs *= rhs; }
  friend complex operator/(complex lhs, const complex& rhs) { return lhs /= rhs; }
  friend complex operator^(complex lhs, const complex& rhs) { return lhs ^= rhs; }

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
    T root_theta = (theta + T(2.0) * T(M_PI) * T(k)) / T(n);
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

} 