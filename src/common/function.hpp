#pragma once

#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {

template <typename T = double>
class func {
  static_assert(std::is_floating_point<T>::value, "numc::func mozna uzywac TYLKO z typami zmiennoprzecinkowymi!");

 private:
  std::function<T(T)> f;

 public:
  template <typename F>
  func(F func) : f(func) {}

  static func x() {
    return func([](T val) { return val; });
  }

  [[nodiscard]] T operator()(T x) const { return f(x); }

  func operator-() const {
    auto func = this->f;
    return func([func](T val) { return -func(val); });
  }

  friend func operator+(const func& lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return lhs(val) + rhs(val); });
  }
  friend func operator-(const func& lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return lhs(val) - rhs(val); });
  }
  friend func operator*(const func& lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return lhs(val) * rhs(val); });
  }
  friend func operator/(const func& lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return lhs(val) / rhs(val); });
  }
  friend func operator^(const func& lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return std::pow(lhs(val), rhs(val)); });
  }

  friend func operator+(const func& lhs, T rhs) {
    return func([lhs, rhs](T val) { return lhs(val) + rhs; });
  }
  friend func operator-(const func& lhs, T rhs) {
    return func([lhs, rhs](T val) { return lhs(val) - rhs; });
  }
  friend func operator*(const func& lhs, T rhs) {
    return func([lhs, rhs](T val) { return lhs(val) * rhs; });
  }
  friend func operator/(const func& lhs, T rhs) {
    return func([lhs, rhs](T val) { return lhs(val) / rhs; });
  }
  friend func operator^(const func& lhs, T rhs) {
    return func([lhs, rhs](T val) { return std::pow(lhs(val), rhs); });
  }

  friend func operator+(T lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return lhs + rhs(val); });
  }
  friend func operator-(T lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return lhs - rhs(val); });
  }
  friend func operator*(T lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return lhs * rhs(val); });
  }
  friend func operator/(T lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return lhs / rhs(val); });
  }
  friend func operator^(T lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return std::pow(lhs, rhs(val)); });
  }
};

template <typename T>
inline func<T> pow(const func<T>& lhs, const func<T>& rhs) {
  return lhs ^ rhs;
}

template <typename T>
inline func<T> pow(const func<T>& lhs, T rhs) {
  return lhs ^ rhs;
}

template <typename T>
inline func<T> pow(T lhs, const func<T>& rhs) {
  return lhs ^ rhs;
}

// Przydatne funkcje trygonometryczne!
// template <typename T>
// inline func<T> sin(const func<T>& f) {
//   return func<T>([f](T val) { return std::sin(f(val)); });
// }

// template <typename T>
// inline func<T> cos(const func<T>& f) {
//   return func<T>([f](T val) { return std::cos(f(val)); });
// }

// template <typename T>
// inline func<T> exp(const func<T>& f) {
//   return func<T>([f](T val) { return std::exp(f(val)); });
// }

}  // namespace numc