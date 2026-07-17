#pragma once

#include "../common/vector.hpp"
#include "../common/tensor.hpp"
#include "../linear_algebra/eigen.hpp"
#include "../linear_algebra/matrix_ops.hpp"
#include "../statistics/stats.hpp"
#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {
namespace machine_learning {

/// @addtogroup machine_learning
/// @{

/// @brief Struct representing the result of Principal Component Analysis.
template <typename T = double>
struct PCAResult {
  tensor<T> components;     // Principal axes in feature space, representing the directions of maximum variance.
  vector<T> explained_variance; // The amount of variance explained by each of the selected components.
  std::vector<vector<T>> transformed_data; // Data transformed to the new subspace.
};

/// @brief Principal Component Analysis (PCA).
/// Performs linear dimensionality reduction using covariance matrix eigenvalue decomposition.
/// @tparam T Floating-point type.
/// @param data A vector of points (each point is a vector<T>), shape: (n_samples, n_features).
/// @param n_components Number of components to keep.
/// @return PCAResult containing components, variance, and transformed data.
template <typename T = double>
PCAResult<T> pca(
    const std::vector<vector<T>>& data,
    size_t n_components) {
  if (data.empty()) return {};
  size_t n_samples = data.size();
  size_t n_features = data[0].size();
  
  if (n_components > n_features) {
    Log::Warn("PCA: n_components cannot be greater than n_features. Setting n_components = n_features.");
    n_components = n_features;
  }

  // 1. Center the data (subtract mean)
  std::vector<vector<T>> centered_data = data;
  vector<T> means(n_features);
  for (size_t j = 0; j < n_features; ++j) {
    T sum = T(0.0);
    for (size_t i = 0; i < n_samples; ++i) {
      sum += data[i][j];
    }
    means[j] = sum / T(n_samples);
  }

  for (size_t i = 0; i < n_samples; ++i) {
    for (size_t j = 0; j < n_features; ++j) {
      centered_data[i][j] -= means[j];
    }
  }

  // 2. Compute Covariance Matrix (n_features x n_features)
  // To use existing cov_matrix which takes variables as vectors, we need to transpose the data format.
  // Actually, we can just build it directly here:
  tensor<T> cov({n_features, n_features});
  for (size_t j1 = 0; j1 < n_features; ++j1) {
    for (size_t j2 = j1; j2 < n_features; ++j2) {
      T sum = T(0.0);
      for (size_t i = 0; i < n_samples; ++i) {
        sum += centered_data[i][j1] * centered_data[i][j2];
      }
      cov(j1, j2) = sum / T(n_samples - 1);
      cov(j2, j1) = cov(j1, j2);
    }
  }

  // 3. Eigendecomposition of covariance matrix
  // Jacobi method works for symmetric positive semi-definite matrices
  auto [eigenvalues, eigenvectors] = numc::linear_algebra::jacobi(cov);

  // 4. Sort eigenvalues and eigenvectors in descending order
  std::vector<std::pair<T, vector<T>>> eigen_pairs(n_features);
  for (size_t j = 0; j < n_features; ++j) {
    vector<T> col(n_features);
    for (size_t i = 0; i < n_features; ++i) {
      col[i] = eigenvectors(i, j);
    }
    eigen_pairs[j] = {eigenvalues[j], col};
  }

  std::sort(eigen_pairs.begin(), eigen_pairs.end(),
      [](const std::pair<T, vector<T>>& a, const std::pair<T, vector<T>>& b) {
        return a.first > b.first;
      });

  // 5. Select top n_components
  tensor<T> components({n_components, n_features});
  vector<T> expl_var(n_components);
  for (size_t i = 0; i < n_components; ++i) {
    expl_var[i] = eigen_pairs[i].first;
    for (size_t j = 0; j < n_features; ++j) {
      components(i, j) = eigen_pairs[i].second[j];
    }
  }

  // 6. Transform the data (centered_data * components^T)
  std::vector<vector<T>> transformed(n_samples, vector<T>(n_components));
  for (size_t i = 0; i < n_samples; ++i) {
    for (size_t j = 0; j < n_components; ++j) {
      T sum = T(0.0);
      for (size_t k = 0; k < n_features; ++k) {
        sum += centered_data[i][k] * components(j, k);
      }
      transformed[i][j] = sum;
    }
  }

  return {components, expl_var, transformed};
}

/// @}

}  // namespace machine_learning
}  // namespace numc
