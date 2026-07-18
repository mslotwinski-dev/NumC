#pragma once

#include "../inc.hpp"
#include "../utility/log.hpp"
#include "vector.hpp"

namespace numc {

/// @addtogroup structures
/// @{

/// @brief Class representing a mathematical function of one or multiple variables.
/// Enables lazy evaluation and the construction of expression trees.
/// @tparam T The underlying floating-point type (e.g., float, double, long double).
template <typename T = double>
class func {
  static_assert(std::is_floating_point<T>::value, "numc::func mozna uzywac TYLKO z typami zmiennoprzecinkowymi!");

 protected:
  std::function<T(const vector<T>&)> f;

 public:
  /// @brief Constructs a function object from a generic callable (e.g., lambda, std::function).
  template <typename F>
  func(F fn) {
    if constexpr (std::is_invocable_r_v<T, F, const vector<T>&>) {
      f = fn;
    } else if constexpr (std::is_invocable_r_v<T, F, T>) {
      f = [fn](const vector<T>& v) -> T { return fn(v.size() > 0 ? v[0] : T(0.0)); };
    } else if constexpr (std::is_invocable_r_v<T, F, T, T>) {
      f = [fn](const vector<T>& v) -> T { return fn(v.size() > 0 ? v[0] : T(0.0), v.size() > 1 ? v[1] : T(0.0)); };
    } else {
      // Fallback compilation error helper
      f = [](const vector<T>&) { return T(0.0); };
    }
  }

  /// @brief Default constructor
  func() : f([](const vector<T>&) { return T(0.0); }) {}

  /// @brief Implicit conversion to std::function
  operator std::function<T(const vector<T>&)>() const { return f; }

  /// @brief Creates an identity function f(x) = x representing the 1st independent variable (x).
  static func x() {
    return func([](const vector<T>& v) { return v.size() > 0 ? v[0] : T(0.0); });
  }

  /// @brief Creates an identity function f(y) = y representing the 2nd independent variable (y).
  static func y() {
    return func([](const vector<T>& v) { return v.size() > 1 ? v[1] : T(0.0); });
  }

  /// @brief Creates an identity function f(z) = z representing the 3rd independent variable (z).
  static func z() {
    return func([](const vector<T>& v) { return v.size() > 2 ? v[2] : T(0.0); });
  }

  /// @brief Evaluates the function at a specific multi-dimensional point.
  [[nodiscard]] T operator()(const vector<T>& x_vec) const { 
    return f(x_vec); 
  }

  /// @brief Evaluates the function at a specific 1D point.
  [[nodiscard]] T operator()(T x) const { 
    return f({x}); 
  }

  /// @brief Evaluates the function at a specific 2D point.
  [[nodiscard]] T operator()(T x, T y) const { 
    return f({x, y}); 
  }

  /// @brief Evaluates the function at a specific 3D point.
  [[nodiscard]] T operator()(T x, T y, T z) const { 
    return f({x, y, z}); 
  }

  /// @brief Evaluates the function on a vector of 1D values element-wise.
  [[nodiscard]] vector<T> map(const vector<T>& x_vec) const {
    vector<T> res(x_vec.size());
    for (size_t i = 0; i < x_vec.size(); ++i) {
      res[i] = f({x_vec[i]});
    }
    return res;
  }

  /// @brief Unary negation operator.
  func operator-() const {
    auto f_copy = this->f;
    return func([f_copy](const vector<T>& val) { return -f_copy(val); });
  }

  friend func operator+(const func& lhs, const func& rhs) {
    return func([lhs, rhs](const vector<T>& val) { return lhs(val) + rhs(val); });
  }

  friend func operator-(const func& lhs, const func& rhs) {
    return func([lhs, rhs](const vector<T>& val) { return lhs(val) - rhs(val); });
  }

  friend func operator*(const func& lhs, const func& rhs) {
    return func([lhs, rhs](const vector<T>& val) { return lhs(val) * rhs(val); });
  }

  friend func operator/(const func& lhs, const func& rhs) {
    return func([lhs, rhs](const vector<T>& val) { return lhs(val) / rhs(val); });
  }

  friend func operator^(const func& lhs, const func& rhs) {
    return func([lhs, rhs](const vector<T>& val) { return std::pow(lhs(val), rhs(val)); });
  }

  friend func operator+(const func& lhs, T rhs) {
    return func([lhs, rhs](const vector<T>& val) { return lhs(val) + rhs; });
  }

  friend func operator-(const func& lhs, T rhs) {
    return func([lhs, rhs](const vector<T>& val) { return lhs(val) - rhs; });
  }

  friend func operator*(const func& lhs, T rhs) {
    return func([lhs, rhs](const vector<T>& val) { return lhs(val) * rhs; });
  }

  friend func operator/(const func& lhs, T rhs) {
    return func([lhs, rhs](const vector<T>& val) { return lhs(val) / rhs; });
  }

  friend func operator^(const func& lhs, T rhs) {
    return func([lhs, rhs](const vector<T>& val) { return std::pow(lhs(val), rhs); });
  }

  friend func operator+(T lhs, const func& rhs) {
    return func([lhs, rhs](const vector<T>& val) { return lhs + rhs(val); });
  }

  friend func operator-(T lhs, const func& rhs) {
    return func([lhs, rhs](const vector<T>& val) { return lhs - rhs(val); });
  }

  friend func operator*(T lhs, const func& rhs) {
    return func([lhs, rhs](const vector<T>& val) { return lhs * rhs(val); });
  }

  friend func operator/(T lhs, const func& rhs) {
    return func([lhs, rhs](const vector<T>& val) { return lhs / rhs(val); });
  }

  friend func operator^(T lhs, const func& rhs) {
    return func([lhs, rhs](const vector<T>& val) { return std::pow(lhs, rhs(val)); });
  }

  /// @brief Calculates the definite integral of the 1D function using Adaptive Simpson's quadrature.
  T integral(T a, T b, T tol = T(1e-6)) const;
  
  /// @brief Calculates the definite integral of the 2D function over a rectangle.
  T integral_2d(T ax, T bx, T ay, T by, T tol = T(1e-4)) const;

  /// @brief Calculates the definite integral of the N-dimensional function using Monte Carlo.
  std::pair<T, T> integral_mc(const vector<T>& a, const vector<T>& b, size_t n_samples) const;
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

/// @} // Koniec grupy common

// Global coordinate functions to allow purely declarative mathematical expressions.
namespace funcv {
  inline const func<double> x = func<double>::x();
  inline const func<double> y = func<double>::y();
  inline const func<double> z = func<double>::z();
}


}  // namespace numc