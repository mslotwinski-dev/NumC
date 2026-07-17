#pragma once

#include "../common/vector.hpp"
#include <vector>
#include <algorithm>
#include <cmath>

namespace numc {
namespace geometry {

/// @addtogroup geometry
/// @{

/// @brief Struct representing a 2D or 3D point (using vector).
template <typename T = double>
using Point = vector<T>;

/// @brief Computes the 2D cross product of OA and OB vectors, i.e., z-component of their 3D cross product.
/// Returns a positive value, if OAB makes a counter-clockwise turn,
/// negative for clockwise turn, and zero if the points are collinear.
template <typename T = double>
T cross_product_2d(const Point<T>& o, const Point<T>& a, const Point<T>& b) {
  return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0]);
}

/// @brief Computes the Convex Hull of a set of 2D points using the Monotone Chain algorithm.
/// Returns a sequence of points forming the hull in counter-clockwise order.
template <typename T = double>
std::vector<Point<T>> convex_hull(std::vector<Point<T>> pts) {
  size_t n = pts.size(), k = 0;
  if (n <= 3) return pts;
  
  std::vector<Point<T>> hull(2 * n);

  // Sort points lexicographically
  std::sort(pts.begin(), pts.end(), [](const Point<T>& a, const Point<T>& b) {
    return a[0] < b[0] || (a[0] == b[0] && a[1] < b[1]);
  });

  // Build lower hull
  for (size_t i = 0; i < n; ++i) {
    while (k >= 2 && cross_product_2d(hull[k - 2], hull[k - 1], pts[i]) <= 0) k--;
    hull[k++] = pts[i];
  }

  // Build upper hull
  for (size_t i = n - 1, t = k + 1; i > 0; --i) {
    while (k >= t && cross_product_2d(hull[k - 2], hull[k - 1], pts[i - 1]) <= 0) k--;
    hull[k++] = pts[i - 1];
  }

  hull.resize(k - 1);
  return hull;
}

/// @brief Checks if a 2D point is inside a polygon using the Ray Casting algorithm.
template <typename T = double>
bool point_in_polygon(const Point<T>& pt, const std::vector<Point<T>>& polygon) {
  bool inside = false;
  size_t n = polygon.size();
  for (size_t i = 0, j = n - 1; i < n; j = i++) {
    T xi = polygon[i][0], yi = polygon[i][1];
    T xj = polygon[j][0], yj = polygon[j][1];
    
    bool intersect = ((yi > pt[1]) != (yj > pt[1]))
        && (pt[0] < (xj - xi) * (pt[1] - yi) / (yj - yi) + xi);
    if (intersect) inside = !inside;
  }
  return inside;
}

/// @brief Computes the shortest distance from a point `pt` to a line segment `v`->`w`.
template <typename T = double>
T distance_point_to_segment(const Point<T>& pt, const Point<T>& v, const Point<T>& w) {
  T l2 = (w - v) * (w - v); // length squared
  if (l2 == T(0.0)) return (pt - v).norm(); // v == w case

  // Consider the line extending the segment, parameterized as v + t (w - v).
  // We find projection of point p onto the line.
  T t = std::max(T(0.0), std::min(T(1.0), ((pt - v) * (w - v)) / l2));
  Point<T> projection = v + t * (w - v);
  return (pt - projection).norm();
}

/// @}

}  // namespace geometry
}  // namespace numc
