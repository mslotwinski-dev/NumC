#pragma once

#include "../../common/function.hpp"
#include "../../common/vector.hpp"
#include "../../utility/log.hpp"

namespace numc {

/// @defgroup visualization Visualization
/// @ingroup statistics
/// @brief Tools for plotting mathematical functions and data using Plotly.js. Provides an easy interface to create interactive plots in the browser
/// with support for multiple functions, custom labels, and automatic handling of edge cases (NaN/Inf). Designed to work seamlessly with numc::func
/// and numc::vector for quick visualization of mathematical expressions and datasets.
/// @{

/// @brief Klasa do generowania interaktywnych wykresów w przeglądarce (via Plotly.js).
/// @tparam T Typ zmiennoprzecinkowy.
template <typename T = double>
class plot {
 private:
  struct Trace {
    std::vector<T> x;
    std::vector<T> y;
    std::string name;
    std::string mode;  // "lines" albo "markers"
  };

  std::vector<Trace> _traces;
  std::string _title;
  std::string _x_label;
  std::string _y_label;

  // Helper function to safely output floats to JS. Converts NaN/Inf to 'null'
  static std::string safe_js_number(T val) {
    if (std::isnan(val) || std::isinf(val)) {
      return "null";
    }
    std::ostringstream oss;
    oss << val;
    return oss.str();
  }

 public:
  /// @brief Konstruktor domyślny.
  plot(std::string title = "numc Plot") : _title(std::move(title)), _x_label("x"), _y_label("y") {}

  /// @brief Ustawia tytuł osi X.
  plot& set_xlabel(const std::string& label) {
    _x_label = label;
    return *this;
  }

  /// @brief Ustawia tytuł osi Y.
  plot& set_ylabel(const std::string& label) {
    _y_label = label;
    return *this;
  }

  /// @brief Dodaje funkcję matematyczną (func) do wykresu.
  /// @param f Funkcja numc::func do narysowania.
  /// @param x_min Początek przedziału.
  /// @param x_max Koniec przedziału.
  /// @param points Liczba punktów (rozdzielczość wykresu).
  /// @return Referencja do obiektu plot (pozwala na chainowanie).
  plot& add(const func<T>& f, T x_min, T x_max, size_t points = 1000) {
    if (points < 2 || x_max <= x_min) {
      Log::Error("Nieprawidlowe parametry przedzialu dla wykresu.");
      return *this;
    }

    Trace tr;
    tr.name = "f(x) = ";  // + f.to_string() Utilizing our expression string tracker!
    tr.mode = "lines";

    T step = (x_max - x_min) / T(points - 1);
    for (size_t i = 0; i < points; ++i) {
      T x_val = x_min + T(i) * step;
      tr.x.push_back(x_val);
      tr.y.push_back(f(x_val));  // Ewaluacja naszego drzewa (Lazy Evaluation) w locie!
    }

    _traces.push_back(std::move(tr));
    return *this;
  }

  /// @brief Dodaje surowe punkty (np. z numc::vector) jako wykres punktowy (scatter).
  plot& add_scatter(const numc::vector<T>& x, const numc::vector<T>& y, const std::string& name = "Data") {
    if (x.size() != y.size()) {
      Log::Error("Wektory X i Y musza miec ten sam rozmiar do narysowania wykresu!");
      return *this;
    }
    Trace tr;
    tr.name = name;
    tr.mode = "markers";
    for (size_t i = 0; i < x.size(); ++i) {
      tr.x.push_back(x[i]);
      tr.y.push_back(y[i]);
    }
    _traces.push_back(std::move(tr));
    return *this;
  }

  /// @brief Generuje plik HTML z wykresem i automatycznie otwiera go w przeglądarce.
  /// @param filename Nazwa pliku wyjściowego (domyślnie "numc_plot.html").
  void show(const std::string& filename = "numc_plot.html") const {
    std::ofstream file(filename);
    if (!file.is_open()) {
      Log::Error("Nie udalo sie otworzyc pliku do zapisu wykresu!");
      return;
    }

    // --- BUDOWANIE PLIKU HTML Z PLOTLY.JS ---
    file << "<!DOCTYPE html>\n<html>\n<head>\n";
    file << "  <meta charset=\"utf-8\">\n";  // <--- POPRAWKA 1: Naprawa 'dziwnego A'
    file << "  <script src=\"https://cdn.plot.ly/plotly-2.27.0.min.js\"></script>\n";
    file << "  <title>" << _title << "</title>\n";
    file << "</head>\n<body>\n";
    file << "  <div id=\"plot_div\" style=\"width:100%;height:100vh;\"></div>\n";
    file << "  <script>\n";

    // --- GENEROWANIE DANYCH ---
    std::vector<std::string> trace_names;
    for (size_t i = 0; i < _traces.size(); ++i) {
      const auto& tr = _traces[i];
      std::string t_name = "trace" + std::to_string(i);
      trace_names.push_back(t_name);

      file << "    var " << t_name << " = {\n";
      file << "      x: [";
      for (size_t j = 0; j < tr.x.size(); ++j) {
        file << safe_js_number(tr.x[j]) << (j + 1 < tr.x.size() ? "," : "");
      }
      file << "],\n";

      file << "      y: [";
      for (size_t j = 0; j < tr.y.size(); ++j) {
        file << safe_js_number(tr.y[j]) << (j + 1 < tr.y.size() ? "," : "");
      }
      file << "],\n";

      file << "      mode: '" << tr.mode << "',\n";
      file << "      name: '" << tr.name << "'\n";
      // file << "      hovertemplate: '(%{x:.3~f}, %{y:.3~f})<extra></extra>'\n";
      file << "    };\n";
    }

    // --- ŁĄCZENIE I RYSOWANIE ---
    file << "    var data = [";
    for (size_t i = 0; i < trace_names.size(); ++i) {
      file << trace_names[i] << (i + 1 < trace_names.size() ? ", " : "");
    }
    file << "];\n";

    file << "    var layout = {\n";
    file << "      title: '" << _title << "',\n";
    file << "      xaxis: { title: '" << _x_label << "' },\n";
    file << "      yaxis: { title: '" << _y_label << "' },\n";
    file << "      template: 'plotly_dark'\n";  // Ciemny motyw dla inżynierów!
    file << "    };\n";

    file << "    Plotly.newPlot('plot_div', data, layout);\n";
    file << "  </script>\n</body>\n</html>\n";
    file.close();

    // --- AUTOMATYCZNE OTWIERANIE W PRZEGLĄDARCE ---
#if defined(_WIN32)
    std::string cmd = "start \"\" \"" + filename + "\"";  // Fixed start command for paths with spaces
#elif defined(__APPLE__)
    std::string cmd = "open \"" + filename + "\"";
#else
    std::string cmd = "xdg-open \"" + filename + "\"";
#endif
    std::system(cmd.c_str());
  }
};

/// @} // Koniec grupy visualization

}  // namespace numc