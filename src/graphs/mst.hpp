#pragma once

#include "pathfinding.hpp"
#include <algorithm>
#include <numeric>

namespace numc {
namespace graphs {

/// @addtogroup graphs
/// @{

/// @brief Disjoint Set Union (Union-Find) data structure used for Kruskal's algorithm.
class DSU {
 private:
  std::vector<int> parent;
  std::vector<int> rank;

 public:
  DSU(size_t n) {
    parent.resize(n);
    rank.resize(n, 0);
    std::iota(parent.begin(), parent.end(), 0);
  }

  int find(int i) {
    if (parent[i] == i)
      return i;
    return parent[i] = find(parent[i]); // Path compression
  }

  void unite(int i, int j) {
    int root_i = find(i);
    int root_j = find(j);
    if (root_i != root_j) {
      if (rank[root_i] < rank[root_j]) {
        parent[root_i] = root_j;
      } else if (rank[root_i] > rank[root_j]) {
        parent[root_j] = root_i;
      } else {
        parent[root_j] = root_i;
        rank[root_i]++;
      }
    }
  }
};

/// @brief Computes the Minimum Spanning Tree (MST) of a graph using Kruskal's Algorithm.
/// Returns the list of edges forming the MST and its total weight.
template <typename T = double>
std::pair<std::vector<std::tuple<size_t, size_t, T>>, T> kruskal_mst(const Graph<T>& graph) {
  size_t n = graph.size();
  std::vector<std::tuple<T, size_t, size_t>> all_edges;

  // Extract all unique edges from the graph
  for (size_t u = 0; u < n; ++u) {
    for (const auto& edge : graph.get_neighbors(u)) {
      size_t v = edge.to;
      if (u < v) { // Assuming undirected graph where both (u,v) and (v,u) exist
        all_edges.push_back({edge.weight, u, v});
      }
    }
  }

  // Sort edges by weight
  std::sort(all_edges.begin(), all_edges.end());

  DSU dsu(n);
  std::vector<std::tuple<size_t, size_t, T>> mst_edges;
  T total_weight = T(0.0);

  for (const auto& [w, u, v] : all_edges) {
    if (dsu.find(u) != dsu.find(v)) {
      dsu.unite(u, v);
      mst_edges.push_back({u, v, w});
      total_weight += w;
    }
  }

  return {mst_edges, total_weight};
}

/// @brief Computes the Minimum Spanning Tree (MST) of a graph using Prim's Algorithm.
/// Returns the list of edges forming the MST and its total weight.
template <typename T = double>
std::pair<std::vector<std::tuple<size_t, size_t, T>>, T> prim_mst(const Graph<T>& graph, size_t start = 0) {
  size_t n = graph.size();
  std::vector<bool> in_mst(n, false);
  std::vector<T> min_weight(n, std::numeric_limits<T>::infinity());
  std::vector<int> parent(n, -1);
  
  std::priority_queue<std::pair<T, size_t>, std::vector<std::pair<T, size_t>>, std::greater<std::pair<T, size_t>>> pq;

  min_weight[start] = T(0.0);
  pq.push({T(0.0), start});

  T total_weight = T(0.0);
  std::vector<std::tuple<size_t, size_t, T>> mst_edges;

  while (!pq.empty()) {
    auto [w, u] = pq.top();
    pq.pop();

    if (in_mst[u]) continue;

    in_mst[u] = true;
    total_weight += w;
    if (parent[u] != -1) {
      mst_edges.push_back({parent[u], u, w});
    }

    for (const auto& edge : graph.get_neighbors(u)) {
      size_t v = edge.to;
      T weight = edge.weight;

      if (!in_mst[v] && weight < min_weight[v]) {
        parent[v] = static_cast<int>(u);
        min_weight[v] = weight;
        pq.push({weight, v});
      }
    }
  }

  return {mst_edges, total_weight};
}

/// @brief Computes PageRank for each node in a directed graph.
/// @param damping The damping factor (usually 0.85).
/// @param tol Tolerance for convergence.
/// @param max_iter Maximum number of iterations.
template <typename T = double>
std::vector<T> pagerank(const Graph<T>& graph, T damping = T(0.85), T tol = T(1e-6), int max_iter = 100) {
  size_t n = graph.size();
  std::vector<T> pr(n, T(1.0) / T(n));
  
  // Calculate out-degrees
  std::vector<size_t> out_degree(n, 0);
  for (size_t u = 0; u < n; ++u) {
    out_degree[u] = graph.get_neighbors(u).size();
  }

  for (int iter = 0; iter < max_iter; ++iter) {
    std::vector<T> next_pr(n, (T(1.0) - damping) / T(n));
    T dangling_sum = T(0.0);

    // Sum up PR of dangling nodes
    for (size_t u = 0; u < n; ++u) {
      if (out_degree[u] == 0) {
        dangling_sum += pr[u];
      }
    }

    // Distribute PR
    for (size_t u = 0; u < n; ++u) {
      next_pr[u] += damping * (dangling_sum / T(n));
      for (const auto& edge : graph.get_neighbors(u)) {
        size_t v = edge.to;
        next_pr[v] += damping * (pr[u] / T(out_degree[u]));
      }
    }

    // Check convergence
    T diff = T(0.0);
    for (size_t i = 0; i < n; ++i) {
      diff += std::abs(next_pr[i] - pr[i]);
    }
    pr = next_pr;
    if (diff < tol) break;
  }

  return pr;
}

/// @}

}  // namespace graphs
}  // namespace numc
