#pragma once

#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {

class complex {
  double x;
  double y;

 public:
  complex(double real = 0.0, double imag = 0.0);
  ~complex();

  complex operator+(const complex& other) const;
  complex operator-(const complex& other) const;
  complex operator*(const complex& other) const;
  complex operator/(const complex& other) const;

  complex& operator+=(const complex& other);
  complex& operator-=(const complex& other);
  complex& operator*=(const complex& other);
  complex& operator/=(const complex& other);

  bool operator==(const complex& other) const;
  bool operator!=(const complex& other) const;
  bool operator<(const complex& other) const;
  bool operator>(const complex& other) const;
  bool operator<=(const complex& other) const;
  bool operator>=(const complex& other) const;

  complex operator~() const;
  complex operator^(const complex& other) const;

  inline double real() const { return x; }
  inline double imag() const { return y; }

  double abs() const;
  double arg() const;

  complex conjugate() const;
  complex pow(const complex& other) const;
  std::vector<complex> roots(int n) const;

  complex round(int places = 0);

  friend std::ostream& operator<<(std::ostream& os, const complex& c);

  friend double Re(const complex& c);
  friend double Im(const complex& c);
};

double Re(const complex& c);
double Im(const complex& c);

std::ostream& operator<<(std::ostream& os, const complex& c);

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
inline complex operator+(T lhs, const complex& rhs) {
  return complex(lhs) + rhs;
}

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
inline complex operator-(T lhs, const complex& rhs) {
  return complex(lhs) - rhs;
}

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
inline complex operator*(T lhs, const complex& rhs) {
  return complex(lhs) * rhs;
}

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
inline complex operator/(T lhs, const complex& rhs) {
  return complex(lhs) / rhs;
}

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
inline complex operator+=(T lhs, const complex& rhs) {
  return complex(lhs) += rhs;
}

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
inline complex operator-=(T lhs, const complex& rhs) {
  return complex(lhs) -= rhs;
}

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
inline complex operator*=(T lhs, const complex& rhs) {
  return complex(lhs) *= rhs;
}

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
inline complex operator/=(T lhs, const complex& rhs) {
  return complex(lhs) /= rhs;
}

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
inline complex operator^(T lhs, const complex& rhs) {
  return complex(lhs) ^ rhs;
}

template <typename T, typename = std::enable_if_t<std::is_arithmetic<T>::value>>
inline complex pow(T lhs, const complex& rhs) {
  return complex(lhs) ^ rhs;
}

const complex i(0.0, 1.0);

}  // namespace numc
