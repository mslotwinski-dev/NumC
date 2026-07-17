#pragma once

#include "../../common/function.hpp"
#include "../../common/vector.hpp"
#include "../../common/autodiff.hpp"
#include <cmath>
#include "../../utility/log.hpp"

namespace numc {
namespace analysis {

/// @addtogroup analysis
/// @{

/// @brief Finds the local minimum of a multivariable function using Gradient Descent.
template <typename T = double>
vector<T> find_minimum(const func<T>& f, const vector<T>& start_pt, T learning_rate = T(0.01), T tol = T(1e-6), int max_iter = 10000) {
  vector<T> x = start_pt;
  for (int iter = 0; iter < max_iter; ++iter) {
    // Wrap func to std::function
    std::function<T(const vector<T>&)> func_wrapper = [&f](const vector<T>& pt) { return f(pt); };
    vector<T> grad = numc::gradient(func_wrapper, x);
    T grad_norm = T(0.0);
    for (size_t i = 0; i < grad.size(); ++i) {
      grad_norm += grad[i] * grad[i];
    }
    grad_norm = std::sqrt(grad_norm);
    
    if (grad_norm < tol) break;

    for (size_t i = 0; i < x.size(); ++i) {
      x[i] -= learning_rate * grad[i];
    }
  }
  return x;
}

/// @}

}  // namespace analysis
}  // namespace numc
