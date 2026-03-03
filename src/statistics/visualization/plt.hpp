#pragma once

#include "../../common/function.hpp"
#include "../../common/vector.hpp"
#include "../../utility/log.hpp"

namespace numc {

/// @addtogroup visualization
/// @{

namespace plt {

/// @brief Quickly plots a single function.
/// @tparam T Floating-point type (deduced).
/// @tparam B1 Type of the lower bound (allows passing int instead of double).
/// @tparam B2 Type of the upper bound (allows passing int instead of double).
/// @param f The numc::func to draw.
/// @param x_min Start of the domain interval.
/// @param x_max End of the domain interval.
/// @param title Title of the plot.
/// @param xlabel X-axis label.
/// @param ylabel Y-axis label.
template <typename T, typename B1, typename B2>
inline void show(const numc::func<T>& f,
                 B1 x_min,
                 B2 x_max,
                 const std::string& title = "numc Plot",
                 const std::string& xlabel = "x",
                 const std::string& ylabel = "f(x)") {
  numc::plot<T> p(title);
  p.set_xlabel(xlabel).set_ylabel(ylabel);
  p.add(f, static_cast<T>(x_min), static_cast<T>(x_max));
  p.show();
}

/// @brief Quickly plots multiple functions on the same graph (e.g., {f, g}).
/// @tparam T Floating-point type (deduced).
/// @tparam B1 Type of the lower bound (allows passing int instead of double).
/// @tparam B2 Type of the upper bound (allows passing int instead of double).
/// @param funcs List of numc::func objects to draw.
/// @param x_min Start of the domain interval.
/// @param x_max End of the domain interval.
/// @param title Title of the plot.
/// @param xlabel X-axis label.
/// @param ylabel Y-axis label.
template <typename T, typename B1, typename B2>
inline void show(std::initializer_list<numc::func<T>> funcs,
                 B1 x_min,
                 B2 x_max,
                 const std::string& title = "numc Plot",
                 const std::string& xlabel = "x",
                 const std::string& ylabel = "f(x)") {
  numc::plot<T> p(title);
  p.set_xlabel(xlabel).set_ylabel(ylabel);
  for (const auto& f : funcs) {
    p.add(f, static_cast<T>(x_min), static_cast<T>(x_max));
  }
  p.show();
}

/// @brief Quickly creates a scatter plot from two data vectors.
/// @tparam T Floating-point type.
/// @param x Vector of X coordinates.
/// @param y Vector of Y coordinates.
/// @param title Title of the plot.
/// @param xlabel X-axis label.
/// @param ylabel Y-axis label.
template <typename T = double>
inline void scatter(const numc::vector<T>& x,
                    const numc::vector<T>& y,
                    const std::string& title = "Data Plot",
                    const std::string& xlabel = "x",
                    const std::string& ylabel = "y") {
  numc::plot<T> p(title);
  p.set_xlabel(xlabel).set_ylabel(ylabel);
  p.add_scatter(x, y, "Data");
  p.show();
}

}  // namespace plt

/// @} // End of visualization group

}  // namespace numc