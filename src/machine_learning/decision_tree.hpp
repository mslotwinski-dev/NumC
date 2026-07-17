#pragma once

#include "../common/vector.hpp"
#include <vector>
#include <algorithm>
#include <memory>
#include <limits>
#include <cmath>

namespace numc {
namespace machine_learning {

/// @addtogroup machine_learning
/// @{

/// @brief Represents a node in the Decision Tree.
template <typename T = double>
struct TreeNode {
  bool is_leaf;
  T value; // Class label (classification) or mean value (regression) for leaf nodes
  size_t feature_idx;
  T threshold;
  std::unique_ptr<TreeNode<T>> left;
  std::unique_ptr<TreeNode<T>> right;

  TreeNode() : is_leaf(false), value(T(0.0)), feature_idx(0), threshold(T(0.0)) {}
};

/// @brief Basic Decision Tree Classifier using Gini Impurity.
template <typename T = double>
class DecisionTreeClassifier {
 private:
  std::unique_ptr<TreeNode<T>> root;
  size_t max_depth;
  size_t min_samples_split;

  // Calculates Gini Impurity for a set of labels
  T gini_impurity(const std::vector<T>& y) {
    if (y.empty()) return T(0.0);
    std::vector<T> sorted_y = y;
    std::sort(sorted_y.begin(), sorted_y.end());
    
    T impurity = T(1.0);
    size_t n = y.size();
    size_t count = 1;
    
    for (size_t i = 1; i <= n; ++i) {
      if (i == n || sorted_y[i] != sorted_y[i - 1]) {
        T p = static_cast<T>(count) / static_cast<T>(n);
        impurity -= p * p;
        count = 1;
      } else {
        count++;
      }
    }
    return impurity;
  }

  std::unique_ptr<TreeNode<T>> build_tree(const std::vector<vector<T>>& X, const std::vector<T>& y, size_t depth) {
    auto node = std::make_unique<TreeNode<T>>();
    
    // Check stopping criteria
    bool all_same = true;
    for (size_t i = 1; i < y.size(); ++i) {
      if (y[i] != y[0]) {
        all_same = false;
        break;
      }
    }

    if (depth >= max_depth || y.size() < min_samples_split || all_same) {
      node->is_leaf = true;
      // Assign the most common class
      if (!y.empty()) {
        std::vector<T> sorted_y = y;
        std::sort(sorted_y.begin(), sorted_y.end());
        T best_val = sorted_y[0];
        size_t max_count = 1, current_count = 1;
        for (size_t i = 1; i < sorted_y.size(); ++i) {
          if (sorted_y[i] == sorted_y[i - 1]) {
            current_count++;
          } else {
            if (current_count > max_count) {
              max_count = current_count;
              best_val = sorted_y[i - 1];
            }
            current_count = 1;
          }
        }
        if (current_count > max_count) best_val = sorted_y.back();
        node->value = best_val;
      }
      return node;
    }

    size_t n_features = X[0].size();
    T best_gini = std::numeric_limits<T>::max();
    size_t best_feature = 0;
    T best_threshold = T(0.0);

    // Find best split
    for (size_t j = 0; j < n_features; ++j) {
      std::vector<T> column(X.size());
      for (size_t i = 0; i < X.size(); ++i) column[i] = X[i][j];
      std::sort(column.begin(), column.end());

      for (size_t i = 1; i < column.size(); ++i) {
        if (column[i] == column[i - 1]) continue;
        T threshold = (column[i] + column[i - 1]) / T(2.0);

        std::vector<T> left_y, right_y;
        for (size_t k = 0; k < X.size(); ++k) {
          if (X[k][j] <= threshold) left_y.push_back(y[k]);
          else right_y.push_back(y[k]);
        }

        T p_left = static_cast<T>(left_y.size()) / static_cast<T>(y.size());
        T p_right = static_cast<T>(right_y.size()) / static_cast<T>(y.size());
        T gini = p_left * gini_impurity(left_y) + p_right * gini_impurity(right_y);

        if (gini < best_gini) {
          best_gini = gini;
          best_feature = j;
          best_threshold = threshold;
        }
      }
    }

    if (best_gini == std::numeric_limits<T>::max()) {
       node->is_leaf = true;
       node->value = y[0];
       return node;
    }

    node->feature_idx = best_feature;
    node->threshold = best_threshold;

    std::vector<vector<T>> left_X, right_X;
    std::vector<T> left_y, right_y;
    for (size_t i = 0; i < X.size(); ++i) {
      if (X[i][best_feature] <= best_threshold) {
        left_X.push_back(X[i]);
        left_y.push_back(y[i]);
      } else {
        right_X.push_back(X[i]);
        right_y.push_back(y[i]);
      }
    }

    node->left = build_tree(left_X, left_y, depth + 1);
    node->right = build_tree(right_X, right_y, depth + 1);

    return node;
  }

  T predict_node(const vector<T>& x, const TreeNode<T>* node) const {
    if (node->is_leaf) return node->value;
    if (x[node->feature_idx] <= node->threshold) {
      return predict_node(x, node->left.get());
    } else {
      return predict_node(x, node->right.get());
    }
  }

 public:
  DecisionTreeClassifier(size_t max_d = 10, size_t min_samples = 2) 
    : max_depth(max_d), min_samples_split(min_samples) {}

  void train(const std::vector<vector<T>>& X, const std::vector<T>& y) {
    if (X.empty() || X.size() != y.size()) return;
    root = build_tree(X, y, 0);
  }

  T predict(const vector<T>& x) const {
    if (!root) return T(0.0);
    return predict_node(x, root.get());
  }
};

/// @}

}  // namespace machine_learning
}  // namespace numc
