#pragma once

#include "../../common/function.hpp"
#include "../../common/vector.hpp"
#include "../../optimalization/solvers/roots.hpp"
#include <stdexcept>

namespace numc {
namespace analysis {

/// @addtogroup analysis
/// @{

/// @brief Finds the root of an implicit function F(x, y) = 0 for a given x.
/// @param f The multivariable function F(x, y).
/// @param x The fixed value of x.
/// @param y_guess The initial guess for y.
/// @param tol Tolerance for the root finding.
/// @return The value of y such that F(x, y) = 0.
template <typename T = double>
T find_implicit_y(const func<T>& f, T x, T y_guess, T tol = T(1e-6)) {
  // We reduce F(x, y) to a 1D function g(y) = F(x, y)
  func<T> g = [&f, x](T y) {
    return f(x, y);
  };
  
  // Use Newton's method or secant. We'll use secant from analysis::optimalization::secant
  // Since we need two guesses for secant, we can use y_guess and y_guess + 1e-4
  T root = optimalization::secant<T>(g, y_guess, y_guess + T(1e-4), tol);
  return root;
}

/// @}

}  // namespace analysis
}  // namespace numc
