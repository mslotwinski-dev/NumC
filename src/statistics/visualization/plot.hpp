#pragma once

#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <sstream>
#include <string>
#include <vector>

#include "../../common/function.hpp"
#include "../../common/vector.hpp"
#include "../../utility/log.hpp"
#include "../histogram.hpp"
#include "plot_renderer.hpp"
#include "plot_window.hpp"

namespace numc {

/// @defgroup visualization Visualization
/// @ingroup statistics
/// @brief Tools for plotting mathematical functions and data. Provides an easy interface to create
/// publication-quality SVG plots (default) or interactive HTML plots (via Plotly.js).
/// Supports multiple functions, custom labels, scatter data, and automatic handling of edge cases (NaN/Inf).
/// Designed to work seamlessly with numc::func and numc::vector for quick visualization of mathematical expressions and datasets.
/// @{

/// @brief Backend selection for plot rendering.
enum class PlotBackend {
  WINDOW, ///< Default. Opens a native GUI window with the plot.
  SVG,    ///< Generates a beautiful matplotlib-style SVG file.
  WEB     ///< Generates an interactive HTML file using Plotly.js.
};

/// @brief A class for generating publication-quality plots.
/// @tparam T Underlying floating-point type.
template <typename T = double>
class plot {
 private:
  struct Trace {
    std::vector<T> x;
    std::vector<T> y;
    std::string name;
    std::string mode;  // "lines" or "markers"
    std::string color;
  };

  static constexpr const char* TRACE_COLORS[] = {
      "#58a6ff", "#f78166", "#3fb950", "#d2a8ff",
      "#f0883e", "#79c0ff", "#ff7b72", "#7ee787",
  };
  static constexpr size_t TRACE_COLOR_COUNT = 8;

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
  /// @brief Default constructor.
  /// @param title Title of the plot.
  plot(std::string title = "numc Plot") : _title(std::move(title)), _x_label("x"), _y_label("y") {}

  /// @brief Sets the plot title.
  plot& set_title(const std::string& title) {
    _title = title;
    return *this;
  }

  /// @brief Sets the label for the X axis.
  plot& set_xlabel(const std::string& label) {
    _x_label = label;
    return *this;
  }

  /// @brief Sets the label for the Y axis.
  plot& set_ylabel(const std::string& label) {
    _y_label = label;
    return *this;
  }

  /// @brief Dodaje funkcję matematyczną (func) do wykresu.
  /// @param f Funkcja numc::func do narysowania.
  /// @param x_min Początek przedziału.
  /// @param x_max Koniec przedziału.
  /// @param name Niestandardowa nazwa do legendy (domyślnie pusta, co wygeneruje "y1", "y2" itd.).
  /// @param points Liczba punktów (rozdzielczość wykresu).
  /// @return Referencja do obiektu plot (pozwala na chainowanie).
  plot& add(const func<T>& f, T x_min, T x_max, const std::string& name = "", size_t points = 1000) {
    if (points < 2 || x_max <= x_min) {
      Log::Error("Nieprawidlowe parametry przedzialu dla wykresu.");
      return *this;
    }

    Trace tr;
    // --- NOWA, CZYSTA LOGIKA NAZEWNICTWA ---
    // Jeśli nie podano nazwy, generujemy y1, y2, y3 na podstawie liczby już dodanych wykresów.
    tr.name = name.empty() ? ("y" + std::to_string(_traces.size() + 1)) : name;
    tr.mode = "lines";
    tr.color = TRACE_COLORS[_traces.size() % TRACE_COLOR_COUNT];

    T step = (x_max - x_min) / T(points - 1);
    for (size_t i = 0; i < points; ++i) {
      T x_val = x_min + T(i) * step;
      tr.x.push_back(x_val);
      tr.y.push_back(f(x_val));  // Lazy Evaluation w locie
    }

    _traces.push_back(std::move(tr));
    return *this;
  }

  /// @brief Adds raw data points (e.g., from numc::vector) as a scatter plot.
  /// @param x Vector of X coordinates.
  /// @param y Vector of Y coordinates.
  /// @param name The name of the dataset for the legend.
  plot& add_scatter(const numc::vector<T>& x, const numc::vector<T>& y, const std::string& name = "Data") {
    if (x.size() != y.size()) {
      Log::Error("X and Y vectors must have the same size to plot!");
      return *this;
    }
    Trace tr;
    tr.name = name;
    tr.mode = "markers";
    tr.color = TRACE_COLORS[_traces.size() % TRACE_COLOR_COUNT];
    for (size_t i = 0; i < x.size(); ++i) {
      tr.x.push_back(x[i]);
      tr.y.push_back(y[i]);
    }
    _traces.push_back(std::move(tr));
    return *this;
  }

  /// @brief Alias for add_scatter — same beautiful chain, shorter name.
  plot& add(const numc::vector<T>& x, const numc::vector<T>& y, const std::string& name = "Data") {
    return add_scatter(x, y, name);
  }

  /// @brief Dodaje histogram 1D do wykresu.
  /// @param hist Histogram do narysowania.
  /// @param name Nazwa w legendzie.
  plot& add(const statistics::Histogram1D<T>& hist, const std::string& name = "Histogram") {
    vector<T> centers = hist.get_bin_centers();
    vector<T> contents = hist.get_bin_contents();
    // Plotly supports "bar" traces, but our Trace struct currently assumes "lines" or "markers".
    // We can use lines to draw the outline, or just add_scatter for now, but drawing it as lines looks better for PDFs.
    Trace tr;
    tr.name = name;
    tr.mode = "lines";
    tr.color = TRACE_COLORS[_traces.size() % TRACE_COLOR_COUNT];
    // For a step-like histogram look:
    T hw = (centers[1] - centers[0]) / T(2.0);
    for (size_t i = 0; i < centers.size(); ++i) {
      tr.x.push_back(centers[i] - hw); tr.y.push_back(contents[i]);
      tr.x.push_back(centers[i] + hw); tr.y.push_back(contents[i]);
    }
    _traces.push_back(std::move(tr));
    return *this;
  }

  /// @brief Generates a publication-quality plot.
  /// @param backend The rendering backend (default: WINDOW).
  /// @param theme The visual theme (default: DARK).
  /// @param filename Output filename for SVG/WEB backends (optional).
  void show(PlotBackend backend = PlotBackend::WINDOW, 
            PlotTheme theme = PlotTheme::DARK,
            std::string filename = "") const {
    if (backend == PlotBackend::WEB) {
      _show_web(filename);
    } else if (backend == PlotBackend::SVG) {
      _show_svg(filename, theme);
    } else {
      _show_window(theme);
    }
  }

  /// @brief Generates a plot using the specified theme (defaults to WINDOW backend).
  void show(PlotTheme theme, std::string filename = "") const {
    show(PlotBackend::WINDOW, theme, filename);
  }

  /// @brief Generates an SVG plot with the specified theme.
  void show_svg(PlotTheme theme = PlotTheme::DARK, std::string filename = "") const {
    _show_svg(filename, theme);
  }

  /// @brief Overload for string backend name and theme
  void show(const std::string& backend_str, PlotTheme theme = PlotTheme::DARK, std::string filename = "") const {
    if (backend_str == "web" || backend_str == "WEB" || backend_str == "html" || backend_str == "plotly") {
      _show_web(filename);
    } else if (backend_str == "svg" || backend_str == "SVG") {
      _show_svg(filename, theme);
    } else {
      _show_window(theme);
    }
  }



  /// @brief Explicit method to show the interactive web plot (Plotly.js).
  void show_web(std::string filename = "") const {
    _show_web(filename);
  }

 private:
  // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  // Native Window Backend — GDI+ renderer
  // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  void _show_window(PlotTheme theme) const {
#if defined(_WIN32)
    win_detail::PlotWindowData data;
    data.title = _title;
    data.xlabel = _x_label;
    data.ylabel = _y_label;
    data.theme = theme;

    for (const auto& tr : _traces) {
      svg_detail::TraceData<T> td;
      td.x = tr.x;
      td.y = tr.y;
      td.name = tr.name;
      td.mode = tr.mode;
      td.color = tr.color;
      data.traces.push_back(std::move(td));
    }

    win_detail::show_window(data);
#else
    // Fallback on non-Windows
    Log::Warning("Native window rendering is only supported on Windows. Falling back to SVG.");
    _show_svg("", theme);
#endif
  }

  // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  // SVG Backend — matplotlib-style publication-quality vector graphics
  // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  void _show_svg(std::string filename, PlotTheme theme) const {
    // Build TraceData from internal Trace format
    std::vector<svg_detail::TraceData<T>> trace_data;
    for (const auto& tr : _traces) {
      svg_detail::TraceData<T> td;
      td.x = tr.x;
      td.y = tr.y;
      td.name = tr.name;
      td.mode = tr.mode;
      td.color = tr.color;
      trace_data.push_back(std::move(td));
    }

    svg_detail::SvgPlotRenderer<T> renderer(_title, _x_label, _y_label, trace_data, theme);

    // Determine output path
    std::filesystem::path filepath;
    if (filename.empty()) {
      static int svg_counter = 0;
      svg_counter++;
      std::time_t now = std::time(nullptr);
      std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "numc";
      if (!std::filesystem::exists(temp_dir)) {
        std::filesystem::create_directories(temp_dir);
      }
      std::string temp_filename = "numc_plot_" + std::to_string(now) + "_" + std::to_string(svg_counter) + ".svg";
      filepath = temp_dir / temp_filename;
    } else {
      filepath = filename;
    }

    // Render and write
    std::ofstream file(filepath);
    if (!file.is_open()) {
      Log::Error("Failed to open file for writing the SVG plot!");
      return;
    }
    file << renderer.render();
    file.close();

    // Open in system viewer
    std::string path_str = filepath.string();
#if defined(_WIN32)
    std::string cmd = "start \"\" \"" + path_str + "\"";
#elif defined(__APPLE__)
    std::string cmd = "open \"" + path_str + "\"";
#else
    std::string cmd = "xdg-open \"" + path_str + "\"";
#endif
    std::system(cmd.c_str());
  }

  // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  // Web Backend — interactive Plotly.js HTML (original implementation)
  // ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━
  void _show_web(std::string filename) const {
    std::filesystem::path filepath;

    std::time_t now = std::time(nullptr);
    std::string current_time = std::to_string(now);

    // --- HANDLE TEMP FILE GENERATION ---
    if (filename.empty()) {
      static int plot_counter = 0;
      plot_counter++;

      // Zbudowanie poprawnej ścieżki (użycie operatora '/' dodaje odpowiedni separator systemu)
      std::filesystem::path temp_dir = std::filesystem::temp_directory_path() / "numc";

      // Upewnienie się, że folder 'numc' istnieje w katalogu tymczasowym
      if (!std::filesystem::exists(temp_dir)) {
          std::filesystem::create_directories(temp_dir);
      }

      // Zbudowanie pełnej ścieżki do pliku
      std::string temp_filename = "numc_plot_" + current_time + ".html";
      filepath = temp_dir / temp_filename;
    } else {
      filepath = filename;
    }

    std::ofstream file(filepath);
    if (!file.is_open()) {
      Log::Error("Failed to open file for writing the plot!");
      return;
    }

    // --- BUILD HTML FILE WITH PLOTLY.JS (DARK MODE) ---
    file << "<!DOCTYPE html>\n<html>\n<head>\n";
    file << "  <meta charset=\"utf-8\">\n";
    file << "  <script src=\"https://cdn.plot.ly/plotly-2.27.0.min.js\"></script>\n";
    file << "  <title>" << _title << "</title>\n";
    // set icon
    file << "  <link rel=\"icon\" type=\"image/png\" "
            "href=\"https://raw.githubusercontent.com/mslotwinski-dev/NumC/refs/heads/main/public/readme_icon.png\">\n";
    file << "</head>\n";
    file << "<body style=\"background-color: #121212; margin: 0; padding: 0;\">\n";
    file << "  <div id=\"plot_div\" style=\"width:100%;height:100vh;\"></div>\n";
    file << "  <script>\n";

    // --- GENERATE DATA ARRAYS ---
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
      file << "      name: '" << tr.name << "',\n";
      if (!tr.color.empty()) {
        file << "      line: { color: '" << tr.color << "' },\n";
        file << "      marker: { color: '" << tr.color << "', size: 8 },\n";
      }
      file << "      hovertemplate: '(%{x:.3~f}, %{y:.3~f})<extra></extra>'\n";
      file << "    };\n";
    }

    // --- ASSEMBLE DATA AND LAYOUT ---
    file << "    var data = [";
    for (size_t i = 0; i < trace_names.size(); ++i) {
      file << trace_names[i] << (i + 1 < trace_names.size() ? ", " : "");
    }
    file << "];\n";

    file << "    var layout = {\n";
    file << "      title: { text: '" << _title << "', font: { color: '#ffffff' } },\n";
    file << "      paper_bgcolor: '#121212',\n";
    file << "      plot_bgcolor: '#121212',\n";
    file << "      xaxis: { \n";
    file << "        title: { text: '" << _x_label << "', font: { color: '#ffffff' } },\n";
    file << "        tickfont: { color: '#aaaaaa' },\n";
    file << "        gridcolor: '#333333',\n";
    file << "        zerolinecolor: '#666666'\n";
    file << "      },\n";
    file << "      yaxis: { \n";
    file << "        title: { text: '" << _y_label << "', font: { color: '#ffffff' } },\n";
    file << "        tickfont: { color: '#aaaaaa' },\n";
    file << "        gridcolor: '#333333',\n";
    file << "        zerolinecolor: '#666666'\n";
    file << "      },\n";
    file << "      legend: { font: { color: '#ffffff' } },\n";
    file << "      margin: { t: 60, b: 60, l: 60, r: 60 }\n";
    file << "    };\n";

    file << "    Plotly.newPlot('plot_div', data, layout);\n";
    file << "  </script>\n</body>\n</html>\n";
    file.close();

    // --- AUTOMATICALLY OPEN IN BROWSER ---
    std::string path_str = filepath.string();
#if defined(_WIN32)
    std::string cmd = "start \"\" \"" + path_str + "\"";
#elif defined(__APPLE__)
    std::string cmd = "open \"" + path_str + "\"";
#else
    std::string cmd = "xdg-open \"" + path_str + "\"";
#endif
    std::system(cmd.c_str());
  }
};

/// @} // End of visualization group

}  // namespace numc