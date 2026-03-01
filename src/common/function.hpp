#pragma once

#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {

/// @addtogroup structures
/// @{

/// @brief Class representing a mathematical function of a single variable, enabling lazy evaluation and the construction of expression trees.
/// @tparam T The underlying floating-point type (e.g., float, double, long double).
template <typename T = double>
class func {
  /// @brief Static assertion to ensure that the template parameter T is a floating-point type.
  static_assert(std::is_floating_point<T>::value, "numc::func mozna uzywac TYLKO z typami zmiennoprzecinkowymi!");

 private:
  /// @brief The internal callable object representing the mathematical expression.
  std::function<T(T)> f;

 public:
  /// @brief Constructs a function object from a generic callable (e.g., lambda, std::function).
  /// @tparam F The type of the callable.
  /// @param func The callable object to be wrapped.
  template <typename F>
  func(F func) : f(func) {}

  /// @brief Creates an identity function f(x) = x. This serves as the independent variable when building mathematical expressions.
  /// @return A function object representing the independent variable 'x'.
  static func x() {
    return func([](T val) { return val; });
  }

  /// @brief Evaluates the mathematical function at a specific value.
  /// @param x The value of the independent variable at which to evaluate the function.
  /// @return The result of the function evaluation.
  [[nodiscard]] T operator()(T x) const { return f(x); }

  /// @brief Unary negation operator. Creates a new function that returns the negated result of this function.
  /// @return A new function object representing -f(x).
  func operator-() const {
    auto f_copy = this->f;
    return func([f_copy](T val) { return -f_copy(val); });
  }

  /// @brief Addition of two functions.
  friend func operator+(const func& lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return lhs(val) + rhs(val); });
  }

  /// @brief Subtraction of two functions.
  friend func operator-(const func& lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return lhs(val) - rhs(val); });
  }

  /// @brief Multiplication of two functions.
  friend func operator*(const func& lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return lhs(val) * rhs(val); });
  }

  /// @brief Division of two functions.
  friend func operator/(const func& lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return lhs(val) / rhs(val); });
  }

  /// @brief Exponentiation of a function by another function (lhs^rhs).
  friend func operator^(const func& lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return std::pow(lhs(val), rhs(val)); });
  }

  /// @brief Addition of a function and a scalar.
  friend func operator+(const func& lhs, T rhs) {
    return func([lhs, rhs](T val) { return lhs(val) + rhs; });
  }

  /// @brief Subtraction of a scalar from a function.
  friend func operator-(const func& lhs, T rhs) {
    return func([lhs, rhs](T val) { return lhs(val) - rhs; });
  }

  /// @brief Multiplication of a function by a scalar.
  friend func operator*(const func& lhs, T rhs) {
    return func([lhs, rhs](T val) { return lhs(val) * rhs; });
  }

  /// @brief Division of a function by a scalar.
  friend func operator/(const func& lhs, T rhs) {
    return func([lhs, rhs](T val) { return lhs(val) / rhs; });
  }

  /// @brief Exponentiation of a function by a scalar.
  friend func operator^(const func& lhs, T rhs) {
    return func([lhs, rhs](T val) { return std::pow(lhs(val), rhs); });
  }

  /// @brief Addition of a scalar and a function.
  friend func operator+(T lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return lhs + rhs(val); });
  }

  /// @brief Subtraction of a function from a scalar.
  friend func operator-(T lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return lhs - rhs(val); });
  }

  /// @brief Multiplication of a scalar and a function.
  friend func operator*(T lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return lhs * rhs(val); });
  }

  /// @brief Division of a scalar by a function.
  friend func operator/(T lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return lhs / rhs(val); });
  }

  /// @brief Exponentiation of a scalar by a function.
  friend func operator^(T lhs, const func& rhs) {
    return func([lhs, rhs](T val) { return std::pow(lhs, rhs(val)); });
  }
};

/// @brief Computes the power of a function raised to another function. Alias for the ^ operator.
/// @tparam T The underlying floating-point type.
/// @param lhs The base function.
/// @param rhs The exponent function.
/// @return A new function representing the result.
template <typename T>
inline func<T> pow(const func<T>& lhs, const func<T>& rhs) {
  return lhs ^ rhs;
}

/// @brief Computes the power of a function raised to a scalar. Alias for the ^ operator.
/// @tparam T The underlying floating-point type.
/// @param lhs The base function.
/// @param rhs The scalar exponent.
/// @return A new function representing the result.
template <typename T>
inline func<T> pow(const func<T>& lhs, T rhs) {
  return lhs ^ rhs;
}

/// @brief Computes the power of a scalar raised to a function. Alias for the ^ operator.
/// @tparam T The underlying floating-point type.
/// @param lhs The scalar base.
/// @param rhs The exponent function.
/// @return A new function representing the result.
template <typename T>
inline func<T> pow(T lhs, const func<T>& rhs) {
  return lhs ^ rhs;
}

/// @} // Koniec grupy common

}  // namespace numc