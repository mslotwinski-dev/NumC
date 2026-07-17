#pragma once

#include "../common/vector.hpp"
#include <vector>
#include <algorithm>
#include <queue>
#include <unordered_map>

namespace numc {
namespace machine_learning {

/// @addtogroup machine_learning
/// @{

/// @brief K-Nearest Neighbors Classifier.
template <typename T = double>
class KNNClassifier {
 private:
  std::vector<vector<T>> X_train;
  std::vector<T> y_train;
  size_t k_neighbors;

 public:
  KNNClassifier(size_t k = 3) : k_neighbors(k) {}

  void train(const std::vector<vector<T>>& X, const std::vector<T>& y) {
    if (X.size() != y.size()) return;
    X_train = X;
    y_train = y;
  }

  T predict(const vector<T>& x) const {
    if (X_train.empty()) return T(0.0);

    // pair<distance, class_label>
    std::vector<std::pair<T, T>> distances;
    distances.reserve(X_train.size());

    for (size_t i = 0; i < X_train.size(); ++i) {
      vector<T> diff = x - X_train[i];
      T dist = diff * diff; // Squared Euclidean distance
      distances.push_back({dist, y_train[i]});
    }

    // Partially sort to find the k smallest elements
    size_t k_eff = std::min(k_neighbors, distances.size());
    std::nth_element(distances.begin(), distances.begin() + k_eff, distances.end(),
                     [](const std::pair<T, T>& a, const std::pair<T, T>& b) {
                       return a.first < b.first;
                     });

    // Majority voting
    std::unordered_map<T, size_t> counts;
    for (size_t i = 0; i < k_eff; ++i) {
      counts[distances[i].second]++;
    }

    T best_class = T(0.0);
    size_t max_count = 0;
    for (const auto& [label, count] : counts) {
      if (count > max_count) {
        max_count = count;
        best_class = label;
      }
    }

    return best_class;
  }
};

/// @}

}  // namespace machine_learning
}  // namespace numc
