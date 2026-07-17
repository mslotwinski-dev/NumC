#pragma once

#include "../inc.hpp"
#include "../utility/log.hpp"
#include "vector.hpp"
#include "tensor.hpp"

namespace numc {

/// @addtogroup structures
/// @{

/// @brief Dual number for forward-mode automatic differentiation.
/// Represents a value and its derivative: x + ε * dx, where ε^2 = 0.
template <typename T = double>
class dual {
  static_assert(std::is_floating_point<T>::value, "numc::dual supports floating-point types only!");

 private:
  T val;  // function value
  T der;  // derivative value

 public:
  /// @brief Constructs a dual number with a value and derivative.
  dual(T value = T(0.0), T derivative = T(0.0)) : val(value), der(derivative) {}

  /// @brief Returns the function value.
  T value() const { return val; }

  /// @brief Returns the derivative value.
  T deriv() const { return der; }

  /// @brief Creates a dual number representing the independent variable x.
  static dual variable(T x) { return dual(x, T(1.0)); }

  /// @brief Creates a dual number representing a constant.
  static dual constant(T c) { return dual(c, T(0.0)); }

  // Arithmetic operators
  dual operator+(const dual& rhs) const { return dual(val + rhs.val, der + rhs.der); }
  dual operator-(const dual& rhs) const { return dual(val - rhs.val, der - rhs.der); }
  dual operator*(const dual& rhs) const {
    return dual(val * rhs.val, val * rhs.der + der * rhs.val);
  }
  dual operator/(const dual& rhs) const {
    return dual(val / rhs.val, (der * rhs.val - val * rhs.der) / (rhs.val * rhs.val));
  }
  dual operator-() const { return dual(-val, -der); }

  // Scalar operators
  friend dual operator+(T lhs, const dual& rhs) { return dual(lhs + rhs.val, rhs.der); }
  friend dual operator-(T lhs, const dual& rhs) { return dual(lhs - rhs.val, -rhs.der); }
  friend dual operator*(T lhs, const dual& rhs) { return dual(lhs * rhs.val, lhs * rhs.der); }
  friend dual operator/(T lhs, const dual& rhs) {
    return dual(lhs / rhs.val, -lhs * rhs.der / (rhs.val * rhs.val));
  }
  dual operator+(T rhs) const { return dual(val + rhs, der); }
  dual operator-(T rhs) const { return dual(val - rhs, der); }
  dual operator*(T rhs) const { return dual(val * rhs, der * rhs); }
  dual operator/(T rhs) const { return dual(val / rhs, der / rhs); }

  // Comparison
  bool operator<(const dual& rhs) const { return val < rhs.val; }
  bool operator>(const dual& rhs) const { return val > rhs.val; }
  bool operator==(const dual& rhs) const { return val == rhs.val; }

  /// @brief Print
  friend std::ostream& operator<<(std::ostream& os, const dual& d) {
    os << d.val << " + " << d.der << "ε";
    return os;
  }
};

// Mathematical functions for dual numbers (chain rule applied automatically)

template <typename T>
dual<T> sin(const dual<T>& x) {
  return dual<T>(std::sin(x.value()), x.deriv() * std::cos(x.value()));
}

template <typename T>
dual<T> cos(const dual<T>& x) {
  return dual<T>(std::cos(x.value()), -x.deriv() * std::sin(x.value()));
}

template <typename T>
dual<T> tan(const dual<T>& x) {
  T c = std::cos(x.value());
  return dual<T>(std::tan(x.value()), x.deriv() / (c * c));
}

template <typename T>
dual<T> exp(const dual<T>& x) {
  T e = std::exp(x.value());
  return dual<T>(e, x.deriv() * e);
}

template <typename T>
dual<T> log(const dual<T>& x) {
  return dual<T>(std::log(x.value()), x.deriv() / x.value());
}

template <typename T>
dual<T> sqrt(const dual<T>& x) {
  T s = std::sqrt(x.value());
  return dual<T>(s, x.deriv() / (T(2.0) * s));
}

template <typename T>
dual<T> pow(const dual<T>& x, T n) {
  return dual<T>(std::pow(x.value(), n), x.deriv() * n * std::pow(x.value(), n - T(1.0)));
}

template <typename T>
dual<T> pow(const dual<T>& x, const dual<T>& y) {
  T v = std::pow(x.value(), y.value());
  T d = v * (y.deriv() * std::log(x.value()) + y.value() * x.deriv() / x.value());
  return dual<T>(v, d);
}

template <typename T>
dual<T> abs(const dual<T>& x) {
  return dual<T>(std::abs(x.value()), x.deriv() * (x.value() >= T(0.0) ? T(1.0) : T(-1.0)));
}

template <typename T>
dual<T> asin(const dual<T>& x) {
  return dual<T>(std::asin(x.value()), x.deriv() / std::sqrt(T(1.0) - x.value() * x.value()));
}

template <typename T>
dual<T> acos(const dual<T>& x) {
  return dual<T>(std::acos(x.value()), -x.deriv() / std::sqrt(T(1.0) - x.value() * x.value()));
}

template <typename T>
dual<T> atan(const dual<T>& x) {
  return dual<T>(std::atan(x.value()), x.deriv() / (T(1.0) + x.value() * x.value()));
}

template <typename T>
dual<T> sinh(const dual<T>& x) {
  return dual<T>(std::sinh(x.value()), x.deriv() * std::cosh(x.value()));
}

template <typename T>
dual<T> cosh(const dual<T>& x) {
  return dual<T>(std::cosh(x.value()), x.deriv() * std::sinh(x.value()));
}

template <typename T>
dual<T> tanh(const dual<T>& x) {
  T c = std::cosh(x.value());
  return dual<T>(std::tanh(x.value()), x.deriv() / (c * c));
}

/// @brief Computes the exact derivative of a function using automatic differentiation.
/// Usage: auto_diff([](auto x) { return sin(x) * x; }, 1.0);
template <typename F, typename T = double>
T auto_diff(F f, T x) {
  dual<T> result = f(dual<T>::variable(x));
  return result.deriv();
}

/// @brief Computes the numerical gradient of a vector function F: R^n -> R.
template <typename F, typename T = double>
vector<T> gradient(F F_func, const vector<T>& x, T h = T(1e-7)) {
  size_t n = x.size();
  vector<T> grad(n);
  for (size_t i = 0; i < n; ++i) {
    vector<T> x_plus = x, x_minus = x;
    x_plus[i] += h;
    x_minus[i] -= h;
    grad[i] = (F_func(x_plus) - F_func(x_minus)) / (T(2.0) * h);
  }
  return grad;
}

/// @brief Computes the numerical Jacobian matrix of F: R^n -> R^m.
template <typename F, typename T = double>
tensor<T> jacobian(F F_func, const vector<T>& x, T h = T(1e-7)) {
  size_t n = x.size();
  vector<T> F0 = F_func(x);
  size_t m = F0.size();
  tensor<T> J({m, n});
  for (size_t j = 0; j < n; ++j) {
    vector<T> x_plus = x, x_minus = x;
    x_plus[j] += h;
    x_minus[j] -= h;
    vector<T> Fp = F_func(x_plus), Fm = F_func(x_minus);
    for (size_t i = 0; i < m; ++i) {
      J(i, j) = (Fp[i] - Fm[i]) / (T(2.0) * h);
    }
  }
  return J;
}

/// @brief Computes the numerical Hessian matrix of F: R^n -> R.
template <typename F, typename T = double>
tensor<T> hessian(F F_func, const vector<T>& x, T h = T(1e-5)) {
  size_t n = x.size();
  tensor<T> H({n, n});
  T f0 = F_func(x);
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = i; j < n; ++j) {
      vector<T> xpp = x, xpm = x, xmp = x, xmm = x;
      xpp[i] += h; xpp[j] += h;
      xpm[i] += h; xpm[j] -= h;
      xmp[i] -= h; xmp[j] += h;
      xmm[i] -= h; xmm[j] -= h;
      H(i, j) = (F_func(xpp) - F_func(xpm) - F_func(xmp) + F_func(xmm)) / (T(4.0) * h * h);
      H(j, i) = H(i, j);
    }
  }
  return H;
}

/// @}

}  // namespace numc
