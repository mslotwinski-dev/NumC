#pragma once

#include "../common/vector.hpp"
#include "../common/tensor.hpp"
#include "../statistics/random.hpp"
#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {
namespace machine_learning {

/// @addtogroup machine_learning
/// @{

/// @brief Struct representing the result of a K-Means clustering.
template <typename T = double>
struct KMeansResult {
  std::vector<vector<T>> centroids;
  std::vector<size_t> labels;
  T inertia; // Sum of squared distances of samples to their closest cluster center
};

/// @brief K-Means clustering algorithm.
/// @tparam T Floating-point type.
/// @param data A vector of points (each point is a vector<T>).
/// @param k The number of clusters.
/// @param max_iter Maximum number of iterations.
/// @param tol Tolerance for early stopping (relative change in inertia).
/// @param seed Random seed for initial centroids selection.
/// @return KMeansResult containing centroids, labels, and inertia.
template <typename T = double>
KMeansResult<T> kmeans(
    const std::vector<vector<T>>& data,
    size_t k,
    int max_iter = 300,
    T tol = T(1e-4),
    unsigned seed = 42) {
  if (data.empty()) return {};
  size_t n_samples = data.size();
  size_t n_features = data[0].size();
  
  if (k > n_samples) {
    Log::Warn("K-Means: Number of clusters (k) is greater than number of samples. Setting k = n_samples.");
    k = n_samples;
  }

  // Random number generator
  std::mt19937 rng(seed);
  std::uniform_int_distribution<size_t> dist(0, n_samples - 1);

  // Initialize centroids (randomly select k points from data)
  std::vector<vector<T>> centroids(k);
  std::vector<size_t> chosen_indices;
  while (chosen_indices.size() < k) {
    size_t idx = dist(rng);
    if (std::find(chosen_indices.begin(), chosen_indices.end(), idx) == chosen_indices.end()) {
      chosen_indices.push_back(idx);
      centroids[chosen_indices.size() - 1] = data[idx];
    }
  }

  std::vector<size_t> labels(n_samples, 0);
  T best_inertia = std::numeric_limits<T>::max();

  for (int iter = 0; iter < max_iter; ++iter) {
    T current_inertia = T(0.0);
    std::vector<size_t> counts(k, 0);
    std::vector<vector<T>> new_centroids(k, vector<T>(n_features));

    // Assignment step
    for (size_t i = 0; i < n_samples; ++i) {
      T min_dist = std::numeric_limits<T>::max();
      size_t best_cluster = 0;
      for (size_t j = 0; j < k; ++j) {
        vector<T> diff = data[i] - centroids[j];
        T dist_sq = diff * diff;
        if (dist_sq < min_dist) {
          min_dist = dist_sq;
          best_cluster = j;
        }
      }
      labels[i] = best_cluster;
      current_inertia += min_dist;

      // Accumulate for update step
      new_centroids[best_cluster] = new_centroids[best_cluster] + data[i];
      counts[best_cluster]++;
    }

    // Update step
    for (size_t j = 0; j < k; ++j) {
      if (counts[j] > 0) {
        new_centroids[j] = new_centroids[j] / T(counts[j]);
      } else {
        // Handle empty cluster by reinitializing randomly
        new_centroids[j] = data[dist(rng)];
      }
    }

    centroids = new_centroids;

    // Check convergence
    if (std::abs(best_inertia - current_inertia) / best_inertia < tol) {
      best_inertia = current_inertia;
      break;
    }
    best_inertia = current_inertia;
  }

  return {centroids, labels, best_inertia};
}

/// @}

}  // namespace machine_learning
}  // namespace numc
