#pragma once

#include "stats.hpp"
#include "../common/tensor.hpp"
#include "../linear_algebra/matrix_ops.hpp"

namespace numc {
namespace statistics {

/// @addtogroup statistics
/// @{

/// @brief Propagate error through a linear transformation Y = T * X.
/// @param cov_X The covariance matrix of the original variables X.
/// @param transformation_matrix The transformation matrix T (Jacobian).
/// @return The covariance matrix of the new variables Y: Sigma_Y = T * Sigma_X * T^T.
template <typename T = double>
tensor<T> propagate_error(const tensor<T>& cov_X, const tensor<T>& transformation_matrix) {
  // Y = T * X  =>  Sigma_Y = T * Sigma_X * T^T
  // Ensure dimensions match
  if (cov_X.shape()[0] != cov_X.shape()[1] || cov_X.shape()[0] != transformation_matrix.shape()[1]) {
    throw std::invalid_argument("propagate_error: Invalid matrix dimensions.");
  }

  tensor<T> T_sigma = transformation_matrix * cov_X;
  tensor<T> T_transposed = numc::linear_algebra::transpose(transformation_matrix);
  return T_sigma * T_transposed;
}

/// @brief Compute the uncertainties (standard errors) from a covariance matrix.
/// Uncertainty is the square root of the diagonal elements.
template <typename T = double>
std::vector<T> extract_uncertainties(const tensor<T>& cov) {
  if (cov.shape()[0] != cov.shape()[1]) {
    throw std::invalid_argument("extract_uncertainties: Covariance matrix must be square.");
  }
  std::vector<T> u(cov.shape()[0]);
  for (size_t i = 0; i < u.size(); ++i) {
    u[i] = std::sqrt(std::max(T(0.0), cov(i, i)));
  }
  return u;
}

/// @}

}  // namespace statistics
}  // namespace numc
