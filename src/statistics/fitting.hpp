#pragma once

#include "../inc.hpp"
#include "../common/vector.hpp"
#include "../common/tensor.hpp"
#include "param.hpp"
#include "../optimalization/least_squares.hpp"
#include <vector>
#include <functional>

namespace numc {
namespace statistics {

/// @addtogroup statistics
/// @{

/// @brief Fits a model to data using Levenberg-Marquardt and updates the parameter objects directly.
/// @param model A function representing the model: y = f(x). It should capture the parameters by reference.
/// @param x_data The independent variable data.
/// @param y_data The dependent variable data.
/// @param params A list of pointers to the parameters to be optimized.
/// @param tol Tolerance for the LM algorithm.
/// @param max_iter Maximum iterations for the LM algorithm.
template <typename T = double, typename F>
void fit(F model, const vector<T>& x_data, const vector<T>& y_data, std::vector<param<T>*> params, T tol = T(1e-6), int max_iter = 1000) {
  if (x_data.size() != y_data.size()) {
    throw std::invalid_argument("x_data and y_data must have the same size.");
  }
  if (params.empty()) {
    throw std::invalid_argument("No parameters to fit.");
  }

  size_t n_params = params.size();
  size_t n_points = x_data.size();

  // Extract initial values
  vector<T> beta0(n_params);
  for (size_t i = 0; i < n_params; ++i) {
    beta0[i] = params[i]->value();
  }

  // Residuals function
  auto residuals = [&](const vector<T>& beta) {
    // Set parameters temporarily
    for (size_t i = 0; i < n_params; ++i) {
      params[i]->set_value(beta[i]);
    }
    vector<T> res(n_points);
    for (size_t i = 0; i < n_points; ++i) {
      res[i] = y_data[i] - model(x_data[i]);
    }
    return res;
  };

  // Numerical Jacobian function
  auto jacobian = [&](const vector<T>& beta) {
    tensor<T> J({n_points, n_params});
    T h = T(1e-5);
    for (size_t j = 0; j < n_params; ++j) {
      vector<T> beta_plus = beta;
      vector<T> beta_minus = beta;
      beta_plus[j] += h;
      beta_minus[j] -= h;

      vector<T> r_plus = residuals(beta_plus);
      vector<T> r_minus = residuals(beta_minus);

      for (size_t i = 0; i < n_points; ++i) {
        J(i, j) = (r_plus[i] - r_minus[i]) / (T(2.0) * h);
      }
    }
    // Restore beta
    residuals(beta);
    return J;
  };

  // Run Levenberg-Marquardt
  vector<T> best_beta = optimalization::levenberg_marquardt<T>(residuals, jacobian, beta0, tol, max_iter);

  // Set final values
  for (size_t i = 0; i < n_params; ++i) {
    params[i]->set_value(best_beta[i]);
    // Optionally compute and set errors from covariance matrix here!
  }
}

/// @}

}  // namespace statistics
}  // namespace numc
