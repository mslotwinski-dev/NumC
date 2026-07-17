#pragma once

#include "../common/vector.hpp"
#include "../common/function.hpp"
#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {

/// @addtogroup structures
/// @{

/// @brief Class representing a polynomial with coefficients stored from lowest to highest degree.
/// p(x) = c[0] + c[1]*x + c[2]*x^2 + ... + c[n]*x^n
template <typename T = double>
class polynomial {
  static_assert(std::is_floating_point<T>::value, "numc::polynomial supports floating-point types only!");

 private:
  vector<T> coeffs;

 public:
  /// @brief Default constructor (zero polynomial).
  polynomial() : coeffs(1) { coeffs[0] = T(0.0); }

  /// @brief Construct from a vector of coefficients [c0, c1, c2, ...].
  explicit polynomial(const vector<T>& c) : coeffs(c) {}

  /// @brief Construct from an initializer list.
  polynomial(std::initializer_list<T> c) : coeffs(c) {}

  /// @brief Returns the degree of the polynomial.
  size_t degree() const { return coeffs.size() > 0 ? coeffs.size() - 1 : 0; }

  /// @brief Returns the coefficients.
  const vector<T>& coeff() const { return coeffs; }

  /// @brief Access coefficient by index.
  T& operator[](size_t i) { return coeffs[i]; }
  const T& operator[](size_t i) const { return coeffs[i]; }

  /// @brief Evaluates the polynomial at x using Horner's method.
  T operator()(T x) const {
    size_t n = coeffs.size();
    if (n == 0) return T(0.0);
    T result = coeffs[n - 1];
    for (int i = static_cast<int>(n) - 2; i >= 0; --i) {
      result = result * x + coeffs[i];
    }
    return result;
  }

  /// @brief Converts this polynomial to a numc::func<T> for use in analysis algorithms.
  numc::func<T> to_func() const {
    vector<T> c = coeffs;
    return numc::func<T>([c](T x) {
      size_t n = c.size();
      if (n == 0) return T(0.0);
      T result = c[n - 1];
      for (int i = static_cast<int>(n) - 2; i >= 0; --i) {
        result = result * x + c[i];
      }
      return result;
    });
  }

  /// @brief Addition of two polynomials.
  polynomial operator+(const polynomial& other) const {
    size_t maxN = std::max(coeffs.size(), other.coeffs.size());
    vector<T> result(maxN);
    for (size_t i = 0; i < coeffs.size(); ++i) result[i] += coeffs[i];
    for (size_t i = 0; i < other.coeffs.size(); ++i) result[i] += other.coeffs[i];
    return polynomial(result);
  }

  /// @brief Subtraction of two polynomials.
  polynomial operator-(const polynomial& other) const {
    size_t maxN = std::max(coeffs.size(), other.coeffs.size());
    vector<T> result(maxN);
    for (size_t i = 0; i < coeffs.size(); ++i) result[i] += coeffs[i];
    for (size_t i = 0; i < other.coeffs.size(); ++i) result[i] -= other.coeffs[i];
    return polynomial(result);
  }

  /// @brief Multiplication of two polynomials.
  polynomial operator*(const polynomial& other) const {
    size_t n = coeffs.size() + other.coeffs.size() - 1;
    vector<T> result(n);
    for (size_t i = 0; i < coeffs.size(); ++i)
      for (size_t j = 0; j < other.coeffs.size(); ++j)
        result[i + j] += coeffs[i] * other.coeffs[j];
    return polynomial(result);
  }

  /// @brief Scalar multiplication.
  polynomial operator*(T s) const {
    vector<T> result(coeffs.size());
    for (size_t i = 0; i < coeffs.size(); ++i) result[i] = coeffs[i] * s;
    return polynomial(result);
  }

  /// @brief Scalar multiplication (left).
  friend polynomial operator*(T s, const polynomial& p) { return p * s; }

  /// @brief Computes the derivative polynomial.
  polynomial derivative() const {
    if (coeffs.size() <= 1) return polynomial({T(0.0)});
    vector<T> result(coeffs.size() - 1);
    for (size_t i = 1; i < coeffs.size(); ++i) {
      result[i - 1] = coeffs[i] * T(i);
    }
    return polynomial(result);
  }

  /// @brief Computes the antiderivative polynomial (with C = 0).
  polynomial antiderivative(T C = T(0.0)) const {
    vector<T> result(coeffs.size() + 1);
    result[0] = C;
    for (size_t i = 0; i < coeffs.size(); ++i) {
      result[i + 1] = coeffs[i] / T(i + 1);
    }
    return polynomial(result);
  }

  /// @brief Computes the definite integral from a to b.
  T integrate(T a, T b) const {
    polynomial F = antiderivative();
    return F(b) - F(a);
  }

  /// @brief Returns a polynomial with roots shifted by delta: p(x + delta).
  polynomial shift(T delta) const {
    size_t n = coeffs.size();
    // Use Taylor expansion: p(x + delta) = sum p^(k)(delta)/k! * x^k
    vector<T> result(n);
    polynomial current = *this;
    T fact = T(1.0);
    for (size_t k = 0; k < n; ++k) {
      if (k > 0) fact *= T(k);
      result[k] = current(delta) / fact;
      current = current.derivative();
    }
    return polynomial(result);
  }

  /// @brief Pretty prints the polynomial.
  friend std::ostream& operator<<(std::ostream& os, const polynomial& p) {
    bool first = true;
    for (size_t i = p.coeffs.size(); i > 0; --i) {
      size_t idx = i - 1;
      T c = p.coeffs[idx];
      if (std::abs(c) < T(1e-15)) continue;
      if (!first && c > T(0.0)) os << " + ";
      else if (!first && c < T(0.0)) { os << " - "; c = -c; }
      else if (first && c < T(0.0)) { os << "-"; c = -c; }
      first = false;
      if (idx == 0) os << c;
      else if (idx == 1) {
        if (std::abs(c - T(1.0)) > T(1e-15)) os << c << "*";
        os << "x";
      } else {
        if (std::abs(c - T(1.0)) > T(1e-15)) os << c << "*";
        os << "x^" << idx;
      }
    }
    if (first) os << "0";
    return os;
  }
};

/// @}

}  // namespace numc
