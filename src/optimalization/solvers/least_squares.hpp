#pragma once

#include "../../common/vector.hpp"
#include "../../common/tensor.hpp"
#include "../../linear_algebra/solvers.hpp"
#include "../../inc.hpp"
#include "../../utility/log.hpp"

namespace numc {
namespace optimalization {

/// @addtogroup optimize
/// @{

/// @brief Gauss-Newton algorithm for non-linear least squares.
/// Minimizes S(beta) = sum( (y_i - f(x_i, beta))^2 ).
/// @tparam T Floating-point type.
/// @param residuals Function that returns the vector of residuals: r(beta) = y - f(x, beta).
/// @param jacobian Function that returns the Jacobian matrix J_ij = dr_i / dbeta_j.
/// @param beta0 Initial guess for the parameters.
/// @param tol Error tolerance.
/// @param max_iter Maximum number of iterations.
/// @return The optimized parameter vector.
template <typename T = double>
vector<T> gauss_newton(
    std::function<vector<T>(const vector<T>&)> residuals,
    std::function<tensor<T>(const vector<T>&)> jacobian,
    vector<T> beta0,
    T tol = T(1e-6),
    int max_iter = 100) {
  vector<T> beta = beta0;

  for (int iter = 0; iter < max_iter; ++iter) {
    vector<T> r = residuals(beta);
    tensor<T> J = jacobian(beta);

    // We need to solve (J^T * J) * delta_beta = J^T * (-r)
    tensor<T> JT = numc::linear_algebra::transpose(J);
    tensor<T> JTJ = numc::linear_algebra::matmul(JT, J);
    vector<T> JTr = numc::linear_algebra::matvec(JT, r);
    vector<T> neg_JTr(JTr.size());
    for(size_t i=0; i<JTr.size(); ++i) neg_JTr[i] = -JTr[i];

    vector<T> delta = numc::linear_algebra::gauss_elimin(JTJ, neg_JTr);
    beta = beta + delta;

    if (delta.norm() < tol) {
      return beta;
    }
  }

  Log::Warn("Gauss-Newton did not converge within the maximum iterations.");
  return beta;
}

/// @brief Levenberg-Marquardt algorithm for robust non-linear least squares.
/// Interpolates between Gauss-Newton and Gradient Descent.
/// @tparam T Floating-point type.
/// @param residuals Function that returns the vector of residuals: r(beta) = y - f(x, beta).
/// @param jacobian Function that returns the Jacobian matrix J_ij = dr_i / dbeta_j.
/// @param beta0 Initial guess for the parameters.
/// @param lambda0 Initial damping factor.
/// @param tol Error tolerance.
/// @param max_iter Maximum number of iterations.
/// @return The optimized parameter vector.
template <typename T = double>
vector<T> levenberg_marquardt(
    std::function<vector<T>(const vector<T>&)> residuals,
    std::function<tensor<T>(const vector<T>&)> jacobian,
    vector<T> beta0,
    T lambda0 = T(1e-3),
    T tol = T(1e-6),
    int max_iter = 100) {
  vector<T> beta = beta0;
  T lambda = lambda0;

  for (int iter = 0; iter < max_iter; ++iter) {
    vector<T> r = residuals(beta);
    T error = r * r; // sum of squared residuals
    tensor<T> J = jacobian(beta);

    tensor<T> JT = numc::linear_algebra::transpose(J);
    tensor<T> JTJ = numc::linear_algebra::matmul(JT, J);
    vector<T> JTr = numc::linear_algebra::matvec(JT, r);
    vector<T> neg_JTr(JTr.size());
    for(size_t i=0; i<JTr.size(); ++i) neg_JTr[i] = -JTr[i];

    // Add damping: JTJ + lambda * diag(JTJ)
    tensor<T> H = JTJ;
    for (size_t i = 0; i < H.shape()[0]; ++i) {
      H(i, i) += lambda * JTJ(i, i);
    }

    // Protect against zero diagonal
    for (size_t i = 0; i < H.shape()[0]; ++i) {
      if (H(i, i) == T(0.0)) H(i, i) = lambda;
    }

    vector<T> delta;
    try {
      delta = numc::linear_algebra::gauss_elimin(H, neg_JTr);
    } catch (...) {
      // If singular, increase damping and try again
      lambda *= T(10.0);
      continue;
    }

    if (delta.norm() < tol) {
      return beta + delta;
    }

    vector<T> beta_new = beta + delta;
    vector<T> r_new = residuals(beta_new);
    T error_new = r_new * r_new;

    if (error_new < error) {
      // Step accepted
      lambda /= T(10.0);
      beta = beta_new;
    } else {
      // Step rejected
      lambda *= T(10.0);
    }
  }

  Log::Warn("Levenberg-Marquardt did not converge within the maximum iterations.");
  return beta;
}

/// @}

}  // namespace optimalization
}  // namespace numc
