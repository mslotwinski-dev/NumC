#pragma once

#include <initializer_list>
#include <string>

#include "../../common/function.hpp"
#include "../../common/vector.hpp"
#include "../../utility/log.hpp"

namespace numc {

/// @addtogroup visualization
/// @{

namespace plt {

/// @brief Szybko rysuje pojedynczą funkcję.
/// @tparam T Typ zmiennoprzecinkowy (dedukowany z funkcji).
/// @tparam B1 Typ początku przedziału (pozwala na podanie int zamiast double).
/// @tparam B2 Typ końca przedziału (pozwala na podanie int zamiast double).
/// @param f Funkcja numc::func do narysowania.
/// @param x_min Początek przedziału.
/// @param x_max Koniec przedziału.
/// @param title Tytuł wykresu.
/// @param xlabel Podpis osi X.
/// @param ylabel Podpis osi Y.
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

/// @brief Rysuje wiele funkcji na jednym wykresie (np. {f, g}).
/// @tparam T Typ zmiennoprzecinkowy (dedukowany z funkcji).
/// @tparam B1 Typ początku przedziału (pozwala na podanie int zamiast double).
/// @tparam B2 Typ końca przedziału (pozwala na podanie int zamiast double).
/// @param funcs Lista funkcji numc::func do narysowania.
/// @param x_min Początek przedziału.
/// @param x_max Koniec przedziału.
/// @param title Tytuł wykresu.
/// @param xlabel Podpis osi X.
/// @param ylabel Podpis osi Y.
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

/// @brief Szybko rysuje wykres punktowy (scatter) z dwóch wektorów.
/// @tparam T Typ zmiennoprzecinkowy.
/// @param x Wektor wartości X.
/// @param y Wektor wartości Y.
/// @param title Tytuł wykresu.
/// @param xlabel Podpis osi X.
/// @param ylabel Podpis osi Y.
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

/// @}

}  // namespace numc