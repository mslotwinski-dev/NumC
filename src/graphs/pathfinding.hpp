#pragma once

#include "../common/vector.hpp"
#include "../inc.hpp"
#include "../utility/log.hpp"
#include <queue>
#include <unordered_map>
#include <unordered_set>
#include <limits>

namespace numc {
namespace graphs {

/// @addtogroup graphs
/// @{

/// @brief Represents a weighted directed edge.
template <typename T = double>
struct Edge {
  size_t to;
  T weight;
  
  bool operator>(const Edge& other) const {
    return weight > other.weight;
  }
};

/// @brief Graph representation using an adjacency list.
template <typename T = double>
class Graph {
 private:
  size_t num_vertices;
  std::vector<std::vector<Edge<T>>> adj_list;

 public:
  Graph(size_t vertices) : num_vertices(vertices), adj_list(vertices) {}

  void add_edge(size_t from, size_t to, T weight = T(1.0), bool directed = false) {
    if (from >= num_vertices || to >= num_vertices) {
      Log::Warn("Graph::add_edge: Vertex index out of bounds.");
      return;
    }
    adj_list[from].push_back({to, weight});
    if (!directed) {
      adj_list[to].push_back({from, weight});
    }
  }

  size_t size() const { return num_vertices; }
  const std::vector<Edge<T>>& get_neighbors(size_t v) const { return adj_list[v]; }

  /// @brief Performs Breadth-First Search (BFS) starting from vertex `start`.
  /// Returns a vector of visited vertices in the order they were visited.
  std::vector<size_t> bfs(size_t start) const {
    std::vector<size_t> visited_order;
    std::vector<bool> visited(num_vertices, false);
    std::queue<size_t> q;

    q.push(start);
    visited[start] = true;

    while (!q.empty()) {
      size_t curr = q.front();
      q.pop();
      visited_order.push_back(curr);

      for (const auto& edge : adj_list[curr]) {
        if (!visited[edge.to]) {
          visited[edge.to] = true;
          q.push(edge.to);
        }
      }
    }
    return visited_order;
  }

  /// @brief Performs Depth-First Search (DFS) starting from vertex `start`.
  std::vector<size_t> dfs(size_t start) const {
    std::vector<size_t> visited_order;
    std::vector<bool> visited(num_vertices, false);

    std::function<void(size_t)> dfs_visit = [&](size_t v) {
      visited[v] = true;
      visited_order.push_back(v);
      for (const auto& edge : adj_list[v]) {
        if (!visited[edge.to]) {
          dfs_visit(edge.to);
        }
      }
    };

    dfs_visit(start);
    return visited_order;
  }

  /// @brief Finds shortest paths from `start` to all other vertices using Dijkstra's algorithm.
  /// Returns {distances, previous_vertices}. Unreachable vertices have infinite distance.
  std::pair<std::vector<T>, std::vector<int>> dijkstra(size_t start) const {
    std::vector<T> dist(num_vertices, std::numeric_limits<T>::infinity());
    std::vector<int> prev(num_vertices, -1);
    
    // priority_queue stores pairs of (distance, vertex), sorted by smallest distance first
    std::priority_queue<std::pair<T, size_t>, std::vector<std::pair<T, size_t>>, std::greater<std::pair<T, size_t>>> pq;

    dist[start] = T(0.0);
    pq.push({T(0.0), start});

    while (!pq.empty()) {
      auto [d, u] = pq.top();
      pq.pop();

      // If we found a shorter path previously, skip
      if (d > dist[u]) continue;

      for (const auto& edge : adj_list[u]) {
        size_t v = edge.to;
        T weight = edge.weight;

        if (dist[u] + weight < dist[v]) {
          dist[v] = dist[u] + weight;
          prev[v] = static_cast<int>(u);
          pq.push({dist[v], v});
        }
      }
    }

    return {dist, prev};
  }

  /// @brief Reconstructs the shortest path from the `prev` array returned by pathfinding algorithms.
  static std::vector<size_t> reconstruct_path(int target, const std::vector<int>& prev) {
    std::vector<size_t> path;
    for (int at = target; at != -1; at = prev[at]) {
      path.push_back(at);
    }
    std::reverse(path.begin(), path.end());
    if (path.size() == 1 && prev[path[0]] == -1) return {}; // Isolated or start node logic might need care, usually length 1 is just the start node.
    return path;
  }

  /// @brief A* Pathfinding Algorithm from `start` to `goal`.
  /// @param heuristic A function `h(v)` estimating the cost from vertex `v` to the `goal`.
  std::pair<std::vector<T>, std::vector<int>> a_star(size_t start, size_t goal, std::function<T(size_t)> heuristic) const {
    std::vector<T> g_score(num_vertices, std::numeric_limits<T>::infinity());
    std::vector<T> f_score(num_vertices, std::numeric_limits<T>::infinity());
    std::vector<int> prev(num_vertices, -1);
    
    std::priority_queue<std::pair<T, size_t>, std::vector<std::pair<T, size_t>>, std::greater<std::pair<T, size_t>>> pq;

    g_score[start] = T(0.0);
    f_score[start] = heuristic(start);
    pq.push({f_score[start], start});

    while (!pq.empty()) {
      auto [current_f, u] = pq.top();
      pq.pop();

      if (u == goal) break;

      // Skip outdated entries
      if (current_f > f_score[u]) continue;

      for (const auto& edge : adj_list[u]) {
        size_t v = edge.to;
        T tentative_g = g_score[u] + edge.weight;

        if (tentative_g < g_score[v]) {
          prev[v] = static_cast<int>(u);
          g_score[v] = tentative_g;
          f_score[v] = tentative_g + heuristic(v);
          pq.push({f_score[v], v});
        }
      }
    }

    return {g_score, prev};
  }

  /// @brief All-Pairs Shortest Path using Floyd-Warshall Algorithm.
  /// Returns a matrix D where D[i][j] is the shortest distance from i to j.
  std::vector<std::vector<T>> floyd_warshall() const {
    std::vector<std::vector<T>> dist(num_vertices, std::vector<T>(num_vertices, std::numeric_limits<T>::infinity()));
    
    for (size_t i = 0; i < num_vertices; ++i) {
      dist[i][i] = T(0.0);
      for (const auto& edge : adj_list[i]) {
        dist[i][edge.to] = std::min(dist[i][edge.to], edge.weight); // Take min in case of multigraph
      }
    }

    for (size_t k = 0; k < num_vertices; ++k) {
      for (size_t i = 0; i < num_vertices; ++i) {
        for (size_t j = 0; j < num_vertices; ++j) {
          if (dist[i][k] < std::numeric_limits<T>::infinity() && dist[k][j] < std::numeric_limits<T>::infinity()) {
            if (dist[i][k] + dist[k][j] < dist[i][j]) {
              dist[i][j] = dist[i][k] + dist[k][j];
            }
          }
        }
      }
    }
    return dist;
  }
};

/// @}

}  // namespace graphs
}  // namespace numc
