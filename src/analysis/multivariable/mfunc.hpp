#pragma once

#include "../../common/vector.hpp"
#include <functional>

namespace numc {
namespace analysis {

/// @addtogroup analysis
/// @{

/// @brief Class representing a mathematical function of multiple variables.
template <typename T = double>
class mfunc {
 private:
  std::function<T(const vector<T>&)> f_;

 public:
  template <typename F>
  mfunc(F func) : f_(func) {}

  mfunc(const std::function<T(const vector<T>&)>& func) : f_(func) {}

  /// @brief Evaluates the multivariable function at a given point vector.
  T operator()(const vector<T>& x) const {
    return f_(x);
  }

  /// @brief Helper for two-variable function evaluation.
  T operator()(T x, T y) const {
    return f_({x, y});
  }

  /// @brief Helper for three-variable function evaluation.
  T operator()(T x, T y, T z) const {
    return f_({x, y, z});
  }

  /// @brief Calculates the definite integral of the 2D function over a rectangle using Adaptive Simpson's quadrature.
  T integral_2d(T ax, T bx, T ay, T by, T tol = T(1e-4)) const;

  /// @brief Calculates the definite integral of the N-dimensional function using Monte Carlo integration.
  /// @return A pair containing the estimated integral value and its estimated error.
  std::pair<T, T> integral_mc(const vector<T>& a, const vector<T>& b, size_t n_samples) const;
};

/// @}

}  // namespace analysis
}  // namespace numc
