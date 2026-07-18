#pragma once

#include <cmath>
#include <ctime>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstdlib>

namespace numc {

// ─────────────────────────────────────────────────────────────────────────────
// PlotTheme — light or dark color scheme
// ─────────────────────────────────────────────────────────────────────────────

enum class PlotTheme { DARK, LIGHT };

namespace svg_detail {

// ─────────────────────────────────────────────────────────────────────────────
// PlotStyle — visual configuration for the matplotlib-inspired renderer
// ─────────────────────────────────────────────────────────────────────────────

struct PlotStyle {
  // Canvas
  int width  = 1100;
  int height = 680;

  // Margins (left, right, top, bottom)
  int margin_left   = 90;
  int margin_right  = 30;
  int margin_top    = 65;
  int margin_bottom = 70;

  // Colors (populated by theme)
  std::string bg_color;
  std::string plot_bg_top;
  std::string plot_bg_bottom;
  std::string grid_color;
  std::string grid_zero_color;
  std::string axis_color;
  std::string tick_color;
  std::string label_color;
  std::string title_color;
  std::string legend_bg;
  std::string legend_border;
  std::string legend_text;

  // Trace palette
  static constexpr int PALETTE_SIZE = 8;
  const char* palette[PALETTE_SIZE] = {};

  // Fonts
  std::string font_family = "'Inter', 'Segoe UI', 'Helvetica Neue', Arial, sans-serif";
  int title_font_size  = 20;
  int label_font_size  = 14;
  int tick_font_size   = 12;
  int legend_font_size = 13;

  // Line rendering
  double line_width = 2.4;
  double scatter_radius = 4.0;

  // Grid
  double grid_width = 0.7;
  double grid_dash_array = 4.0;

  const char* trace_color(size_t index) const {
    return palette[index % PALETTE_SIZE];
  }

  int plot_width()  const { return width - margin_left - margin_right; }
  int plot_height() const { return height - margin_top - margin_bottom; }
  int plot_x()      const { return margin_left; }
  int plot_y()      const { return margin_top; }

  // ── Factory methods ──

  static PlotStyle dark() {
    PlotStyle s;
    s.bg_color         = "#0d1117";
    s.plot_bg_top      = "#141a24";
    s.plot_bg_bottom   = "#1a2233";
    s.grid_color       = "#21262d";
    s.grid_zero_color  = "#3b4252";
    s.axis_color       = "#484f58";
    s.tick_color        = "#8b949e";
    s.label_color      = "#c9d1d9";
    s.title_color      = "#f0f6fc";
    s.legend_bg        = "#161b22";
    s.legend_border    = "#30363d";
    s.legend_text      = "#c9d1d9";
    s.palette[0] = "#58a6ff";  // Electric blue
    s.palette[1] = "#f78166";  // Warm coral
    s.palette[2] = "#3fb950";  // Fresh green
    s.palette[3] = "#d2a8ff";  // Soft purple
    s.palette[4] = "#f0883e";  // Amber
    s.palette[5] = "#79c0ff";  // Sky blue
    s.palette[6] = "#ff7b72";  // Salmon red
    s.palette[7] = "#7ee787";  // Mint green
    return s;
  }

  static PlotStyle light() {
    PlotStyle s;
    s.bg_color         = "#ffffff";
    s.plot_bg_top      = "#f8f9fb";
    s.plot_bg_bottom   = "#f0f2f5";
    s.grid_color       = "#e1e4e8";
    s.grid_zero_color  = "#9ca3af";
    s.axis_color       = "#6b7280";
    s.tick_color        = "#4b5563";
    s.label_color      = "#1f2937";
    s.title_color      = "#111827";
    s.legend_bg        = "#ffffff";
    s.legend_border    = "#d1d5db";
    s.legend_text      = "#1f2937";
    s.palette[0] = "#2563eb";  // Royal blue
    s.palette[1] = "#dc2626";  // Red
    s.palette[2] = "#16a34a";  // Green
    s.palette[3] = "#9333ea";  // Purple
    s.palette[4] = "#ea580c";  // Orange
    s.palette[5] = "#0891b2";  // Teal
    s.palette[6] = "#e11d48";  // Rose
    s.palette[7] = "#65a30d";  // Lime
    return s;
  }

  static PlotStyle from_theme(PlotTheme theme) {
    return theme == PlotTheme::LIGHT ? light() : dark();
  }
};

// ─────────────────────────────────────────────────────────────────────────────
// Nice Numbers Algorithm — produces clean, human-readable tick values
// (Based on Paul Heckbert's "Nice Numbers for Graph Labels")
// ─────────────────────────────────────────────────────────────────────────────

inline double nice_num(double x, bool round_flag) {
  if (x <= 0.0) return 1.0;
  int exp_val = static_cast<int>(std::floor(std::log10(x)));
  double frac = x / std::pow(10.0, exp_val);
  double nice;

  if (round_flag) {
    if (frac < 1.5)      nice = 1.0;
    else if (frac < 3.0) nice = 2.0;
    else if (frac < 7.0) nice = 5.0;
    else                  nice = 10.0;
  } else {
    if (frac <= 1.0)      nice = 1.0;
    else if (frac <= 2.0) nice = 2.0;
    else if (frac <= 5.0) nice = 5.0;
    else                  nice = 10.0;
  }
  return nice * std::pow(10.0, exp_val);
}

struct TickResult {
  double tick_min;
  double tick_max;
  double tick_spacing;
  std::vector<double> ticks;
};

/// @brief Computes clean, human-readable tick positions for an axis.
/// The key change from the original Heckbert algorithm: we compute tick_spacing
/// directly from the raw range (without rounding the range first via nice_num),
/// which prevents the axis from expanding beyond the data boundaries.
inline TickResult compute_ticks(double data_min, double data_max, int max_ticks = 7) {
  TickResult result;

  if (std::abs(data_max - data_min) < 1e-15) {
    data_min -= 1.0;
    data_max += 1.0;
  }

  double range = data_max - data_min;
  result.tick_spacing = nice_num(range / (max_ticks - 1), true);
  result.tick_min = std::floor(data_min / result.tick_spacing) * result.tick_spacing;
  result.tick_max = std::ceil(data_max / result.tick_spacing) * result.tick_spacing;

  for (double t = result.tick_min; t <= result.tick_max + 0.5 * result.tick_spacing; t += result.tick_spacing) {
    result.ticks.push_back(t);
  }

  return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// Format number for axis labels — remove trailing zeros, handle -0
// ─────────────────────────────────────────────────────────────────────────────

inline std::string format_tick(double val) {
  if (std::abs(val) < 1e-12) val = 0.0;

  std::ostringstream oss;
  if (std::abs(val) >= 1000.0 || (std::abs(val) < 0.01 && val != 0.0)) {
    oss << std::setprecision(3) << val;
  } else {
    oss << std::fixed;
    if (std::abs(val) < 0.1) oss << std::setprecision(3);
    else if (std::abs(val) < 10.0) oss << std::setprecision(2);
    else oss << std::setprecision(1);
    oss << val;
  }

  std::string s = oss.str();
  if (s.find('.') != std::string::npos) {
    size_t last = s.find_last_not_of('0');
    if (last != std::string::npos && s[last] == '.') last--;
    s = s.substr(0, last + 1);
  }

  return s;
}

// ─────────────────────────────────────────────────────────────────────────────
// XML/SVG escape utility
// ─────────────────────────────────────────────────────────────────────────────

inline std::string xml_escape(const std::string& s) {
  std::string out;
  out.reserve(s.size());
  for (char c : s) {
    switch (c) {
      case '&':  out += "&amp;";  break;
      case '<':  out += "&lt;";   break;
      case '>':  out += "&gt;";   break;
      case '"':  out += "&quot;"; break;
      case '\'': out += "&apos;"; break;
      default:   out += c;
    }
  }
  return out;
}

// ─────────────────────────────────────────────────────────────────────────────
// TraceData — data passed to the renderer
// ─────────────────────────────────────────────────────────────────────────────

template <typename T>
struct TraceData {
  std::vector<T> x;
  std::vector<T> y;
  std::string name;
  std::string mode;   // "lines", "markers"
  std::string color;
};

// ─────────────────────────────────────────────────────────────────────────────
// compute_data_range — raw data range without artificial padding
// ─────────────────────────────────────────────────────────────────────────────

template <typename T>
void compute_data_range(const std::vector<TraceData<T>>& traces,
                        T& x_min, T& x_max, T& y_min, T& y_max) {
  x_min = std::numeric_limits<T>::max();
  x_max = std::numeric_limits<T>::lowest();
  y_min = std::numeric_limits<T>::max();
  y_max = std::numeric_limits<T>::lowest();

  for (const auto& tr : traces) {
    for (size_t i = 0; i < tr.x.size(); ++i) {
      T xv = tr.x[i];
      T yv = tr.y[i];
      if (std::isnan(static_cast<double>(xv)) || std::isinf(static_cast<double>(xv)) ||
          std::isnan(static_cast<double>(yv)) || std::isinf(static_cast<double>(yv))) {
        continue;
      }
      if (xv < x_min) x_min = xv;
      if (xv > x_max) x_max = xv;
      if (yv < y_min) y_min = yv;
      if (yv > y_max) y_max = yv;
    }
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// SvgPlotRenderer — the SVG rendering engine
// ─────────────────────────────────────────────────────────────────────────────

template <typename T>
class SvgPlotRenderer {
 public:
  SvgPlotRenderer(const std::string& title,
                  const std::string& xlabel,
                  const std::string& ylabel,
                  const std::vector<TraceData<T>>& traces,
                  PlotTheme theme = PlotTheme::DARK)
      : _title(title), _xlabel(xlabel), _ylabel(ylabel), _traces(traces), _theme(theme) {}

  std::string render() const {
    PlotStyle s = PlotStyle::from_theme(_theme);
    std::ostringstream svg;

    // Compute data ranges
    T x_min_data, x_max_data, y_min_data, y_max_data;
    compute_data_range(_traces, x_min_data, x_max_data, y_min_data, y_max_data);

    // Compute nice ticks
    auto x_ticks = compute_ticks(static_cast<double>(x_min_data), static_cast<double>(x_max_data), 7);
    auto y_ticks = compute_ticks(static_cast<double>(y_min_data), static_cast<double>(y_max_data), 7);

    double xmin = x_ticks.tick_min;
    double xmax = x_ticks.tick_max;
    double ymin = y_ticks.tick_min;
    double ymax = y_ticks.tick_max;

    int pw = s.plot_width();
    int ph = s.plot_height();
    int px = s.plot_x();
    int py = s.plot_y();

    auto map_x = [&](double val) -> double {
      return px + (val - xmin) / (xmax - xmin) * pw;
    };
    auto map_y = [&](double val) -> double {
      return py + ph - (val - ymin) / (ymax - ymin) * ph;
    };

    // ── SVG header ──
    svg << "<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n";
    svg << "<svg xmlns=\"http://www.w3.org/2000/svg\" "
        << "width=\"" << s.width << "\" height=\"" << s.height << "\" "
        << "viewBox=\"0 0 " << s.width << " " << s.height << "\" "
        << "style=\"font-family:" << s.font_family << ";\">\n";

    // ── Defs: gradients and filters ──
    svg << "<defs>\n";

    svg << "  <linearGradient id=\"plotBg\" x1=\"0\" y1=\"0\" x2=\"0\" y2=\"1\">\n";
    svg << "    <stop offset=\"0%\" stop-color=\"" << s.plot_bg_top << "\"/>\n";
    svg << "    <stop offset=\"100%\" stop-color=\"" << s.plot_bg_bottom << "\"/>\n";
    svg << "  </linearGradient>\n";

    // Glow filters for scatter points
    for (size_t i = 0; i < _traces.size(); ++i) {
      if (_traces[i].mode == "markers") {
        std::string color = _traces[i].color.empty() ? s.trace_color(i) : _traces[i].color;
        svg << "  <filter id=\"glow" << i << "\" x=\"-50%\" y=\"-50%\" width=\"200%\" height=\"200%\">\n";
        svg << "    <feGaussianBlur stdDeviation=\"2.0\" result=\"blur\"/>\n";
        svg << "    <feFlood flood-color=\"" << color << "\" flood-opacity=\"0.3\" result=\"color\"/>\n";
        svg << "    <feComposite in=\"color\" in2=\"blur\" operator=\"in\" result=\"shadow\"/>\n";
        svg << "    <feMerge>\n";
        svg << "      <feMergeNode in=\"shadow\"/>\n";
        svg << "      <feMergeNode in=\"SourceGraphic\"/>\n";
        svg << "    </feMerge>\n";
        svg << "  </filter>\n";
      }
    }

    // Line glow filters
    for (size_t i = 0; i < _traces.size(); ++i) {
      if (_traces[i].mode == "lines") {
        std::string color = _traces[i].color.empty() ? s.trace_color(i) : _traces[i].color;
        svg << "  <filter id=\"lineGlow" << i << "\" x=\"-5%\" y=\"-20%\" width=\"110%\" height=\"140%\">\n";
        svg << "    <feGaussianBlur in=\"SourceGraphic\" stdDeviation=\"2.5\" result=\"blur\"/>\n";
        svg << "    <feFlood flood-color=\"" << color << "\" flood-opacity=\"0.2\" result=\"color\"/>\n";
        svg << "    <feComposite in=\"color\" in2=\"blur\" operator=\"in\" result=\"glow\"/>\n";
        svg << "    <feMerge>\n";
        svg << "      <feMergeNode in=\"glow\"/>\n";
        svg << "      <feMergeNode in=\"SourceGraphic\"/>\n";
        svg << "    </feMerge>\n";
        svg << "  </filter>\n";
      }
    }

    svg << "  <clipPath id=\"plotClip\">\n";
    svg << "    <rect x=\"" << px << "\" y=\"" << py
        << "\" width=\"" << pw << "\" height=\"" << ph << "\"/>\n";
    svg << "  </clipPath>\n";

    svg << "  <filter id=\"legendShadow\" x=\"-5%\" y=\"-5%\" width=\"115%\" height=\"115%\">\n";
    svg << "    <feDropShadow dx=\"0\" dy=\"2\" stdDeviation=\"4\" flood-color=\"#000000\" flood-opacity=\"0.4\"/>\n";
    svg << "  </filter>\n";

    svg << "</defs>\n";

    // ── Background ──
    svg << "<rect width=\"" << s.width << "\" height=\"" << s.height
        << "\" rx=\"8\" fill=\"" << s.bg_color << "\"/>\n";

    // ── Plot area background ──
    svg << "<rect x=\"" << px << "\" y=\"" << py
        << "\" width=\"" << pw << "\" height=\"" << ph
        << "\" rx=\"4\" fill=\"url(#plotBg)\"/>\n";

    // ── Grid lines ──
    svg << "<g stroke=\"" << s.grid_color << "\" stroke-width=\"" << s.grid_width
        << "\" stroke-dasharray=\"" << s.grid_dash_array << "\">\n";

    for (double yt : y_ticks.ticks) {
      double sy = map_y(yt);
      if (sy >= py && sy <= py + ph) {
        bool is_zero = std::abs(yt) < 1e-12;
        if (is_zero) {
          svg << "  <line x1=\"" << px << "\" y1=\"" << sy
              << "\" x2=\"" << px + pw << "\" y2=\"" << sy
              << "\" stroke=\"" << s.grid_zero_color
              << "\" stroke-width=\"1.2\" stroke-dasharray=\"none\"/>\n";
        } else {
          svg << "  <line x1=\"" << px << "\" y1=\"" << sy
              << "\" x2=\"" << px + pw << "\" y2=\"" << sy << "\"/>\n";
        }
      }
    }

    for (double xt : x_ticks.ticks) {
      double sx = map_x(xt);
      if (sx >= px && sx <= px + pw) {
        bool is_zero = std::abs(xt) < 1e-12;
        if (is_zero) {
          svg << "  <line x1=\"" << sx << "\" y1=\"" << py
              << "\" x2=\"" << sx << "\" y2=\"" << py + ph
              << "\" stroke=\"" << s.grid_zero_color
              << "\" stroke-width=\"1.2\" stroke-dasharray=\"none\"/>\n";
        } else {
          svg << "  <line x1=\"" << sx << "\" y1=\"" << py
              << "\" x2=\"" << sx << "\" y2=\"" << py + ph << "\"/>\n";
        }
      }
    }

    svg << "</g>\n";

    // ── Axis border ──
    svg << "<rect x=\"" << px << "\" y=\"" << py
        << "\" width=\"" << pw << "\" height=\"" << ph
        << "\" rx=\"4\" fill=\"none\" stroke=\"" << s.axis_color << "\" stroke-width=\"1.2\"/>\n";

    // ── X axis ticks and labels ──
    svg << "<g text-anchor=\"middle\" fill=\"" << s.tick_color
        << "\" font-size=\"" << s.tick_font_size << "\">\n";
    for (double xt : x_ticks.ticks) {
      double sx = map_x(xt);
      if (sx >= px - 1 && sx <= px + pw + 1) {
        svg << "  <line x1=\"" << sx << "\" y1=\"" << py + ph
            << "\" x2=\"" << sx << "\" y2=\"" << py + ph + 6
            << "\" stroke=\"" << s.axis_color << "\" stroke-width=\"1.2\"/>\n";
        svg << "  <text x=\"" << sx << "\" y=\"" << py + ph + 22 << "\">"
            << xml_escape(format_tick(xt)) << "</text>\n";
      }
    }
    svg << "</g>\n";

    // ── Y axis ticks and labels ──
    svg << "<g text-anchor=\"end\" fill=\"" << s.tick_color
        << "\" font-size=\"" << s.tick_font_size << "\">\n";
    for (double yt : y_ticks.ticks) {
      double sy = map_y(yt);
      if (sy >= py - 1 && sy <= py + ph + 1) {
        svg << "  <line x1=\"" << px - 6 << "\" y1=\"" << sy
            << "\" x2=\"" << px << "\" y2=\"" << sy
            << "\" stroke=\"" << s.axis_color << "\" stroke-width=\"1.2\"/>\n";
        svg << "  <text x=\"" << px - 10 << "\" y=\"" << sy + 4 << "\">"
            << xml_escape(format_tick(yt)) << "</text>\n";
      }
    }
    svg << "</g>\n";

    // ── Data traces (clipped) ──
    svg << "<g clip-path=\"url(#plotClip)\">\n";

    for (size_t i = 0; i < _traces.size(); ++i) {
      const auto& tr = _traces[i];
      std::string color = tr.color.empty() ? s.trace_color(i) : tr.color;

      if (tr.mode == "lines") {
        svg << "  <polyline fill=\"none\" stroke=\"" << color
            << "\" stroke-width=\"" << s.line_width
            << "\" stroke-linecap=\"round\" stroke-linejoin=\"round\""
            << " filter=\"url(#lineGlow" << i << ")\""
            << " points=\"";

        bool first = true;
        for (size_t j = 0; j < tr.x.size(); ++j) {
          T xv = tr.x[j];
          T yv = tr.y[j];
          if (std::isnan(static_cast<double>(xv)) || std::isinf(static_cast<double>(xv)) ||
              std::isnan(static_cast<double>(yv)) || std::isinf(static_cast<double>(yv))) {
            continue;
          }
          double sx = map_x(static_cast<double>(xv));
          double sy = map_y(static_cast<double>(yv));
          if (!first) svg << " ";
          svg << std::fixed << std::setprecision(2) << sx << "," << sy;
          first = false;
        }

        svg << "\"/>\n";

      } else if (tr.mode == "markers") {
        svg << "  <g filter=\"url(#glow" << i << ")\">\n";
        for (size_t j = 0; j < tr.x.size(); ++j) {
          T xv = tr.x[j];
          T yv = tr.y[j];
          if (std::isnan(static_cast<double>(xv)) || std::isinf(static_cast<double>(xv)) ||
              std::isnan(static_cast<double>(yv)) || std::isinf(static_cast<double>(yv))) {
            continue;
          }
          double sx = map_x(static_cast<double>(xv));
          double sy = map_y(static_cast<double>(yv));

          svg << "    <circle cx=\"" << std::fixed << std::setprecision(2) << sx
              << "\" cy=\"" << sy << "\" r=\"" << s.scatter_radius + 1.5
              << "\" fill=\"none\" stroke=\"" << color
              << "\" stroke-width=\"1.5\" opacity=\"0.5\"/>\n";
          svg << "    <circle cx=\"" << sx << "\" cy=\"" << sy
              << "\" r=\"" << s.scatter_radius
              << "\" fill=\"" << color << "\" stroke=\"" << s.bg_color
              << "\" stroke-width=\"1.8\"/>\n";
        }
        svg << "  </g>\n";
      }
    }

    svg << "</g>\n";

    // ── Title ──
    svg << "<text x=\"" << s.width / 2 << "\" y=\"" << s.margin_top - 22
        << "\" text-anchor=\"middle\" fill=\"" << s.title_color
        << "\" font-size=\"" << s.title_font_size
        << "\" font-weight=\"600\" letter-spacing=\"0.5\">"
        << xml_escape(_title) << "</text>\n";

    // ── X axis label ──
    svg << "<text x=\"" << px + pw / 2 << "\" y=\"" << s.height - 14
        << "\" text-anchor=\"middle\" fill=\"" << s.label_color
        << "\" font-size=\"" << s.label_font_size << "\">"
        << xml_escape(_xlabel) << "</text>\n";

    // ── Y axis label (rotated) ──
    svg << "<text x=\"" << 22 << "\" y=\"" << py + ph / 2
        << "\" text-anchor=\"middle\" fill=\"" << s.label_color
        << "\" font-size=\"" << s.label_font_size
        << "\" transform=\"rotate(-90, 22, " << py + ph / 2 << ")\">"
        << xml_escape(_ylabel) << "</text>\n";

    // ── Legend ──
    if (!_traces.empty()) {
      render_legend(svg, s);
    }

    // ── Watermark ──
    svg << "<text x=\"" << s.width - 12 << "\" y=\"" << s.height - 8
        << "\" text-anchor=\"end\" fill=\"" << s.axis_color
        << "\" font-size=\"10\" opacity=\"0.5\">NumC</text>\n";

    svg << "</svg>\n";
    return svg.str();
  }

  void render_to_file(const std::string& path) const {
    std::ofstream file(path);
    if (!file.is_open()) return;
    file << render();
    file.close();
  }

 private:
  std::string _title;
  std::string _xlabel;
  std::string _ylabel;
  std::vector<TraceData<T>> _traces;
  PlotTheme _theme;

  void render_legend(std::ostringstream& svg, const PlotStyle& s) const {
    int legend_x = s.plot_x() + s.plot_width() - 10;
    int legend_y = s.plot_y() + 12;
    int item_height = 24;
    int padding_x = 14;
    int padding_y = 10;

    int max_text_width = 0;
    for (const auto& tr : _traces) {
      int w = static_cast<int>(tr.name.size()) * 8;
      if (w > max_text_width) max_text_width = w;
    }

    int box_width = max_text_width + 50;
    int box_height = static_cast<int>(_traces.size()) * item_height + padding_y * 2;
    int box_x = legend_x - box_width - 4;
    int box_y = legend_y;

    svg << "<g filter=\"url(#legendShadow)\">\n";
    svg << "  <rect x=\"" << box_x << "\" y=\"" << box_y
        << "\" width=\"" << box_width << "\" height=\"" << box_height
        << "\" rx=\"6\" fill=\"" << s.legend_bg
        << "\" stroke=\"" << s.legend_border << "\" stroke-width=\"1\" opacity=\"0.92\"/>\n";
    svg << "</g>\n";

    svg << "<g font-size=\"" << s.legend_font_size << "\" fill=\"" << s.legend_text << "\">\n";
    for (size_t i = 0; i < _traces.size(); ++i) {
      const auto& tr = _traces[i];
      std::string color = tr.color.empty() ? s.trace_color(i) : tr.color;
      int iy = box_y + padding_y + static_cast<int>(i) * item_height + 16;

      if (tr.mode == "markers") {
        svg << "  <circle cx=\"" << box_x + padding_x + 6 << "\" cy=\"" << iy - 4
            << "\" r=\"5\" fill=\"" << color << "\"/>\n";
      } else {
        svg << "  <line x1=\"" << box_x + padding_x << "\" y1=\"" << iy - 4
            << "\" x2=\"" << box_x + padding_x + 18 << "\" y2=\"" << iy - 4
            << "\" stroke=\"" << color << "\" stroke-width=\"2.5\" stroke-linecap=\"round\"/>\n";
      }

      svg << "  <text x=\"" << box_x + padding_x + 26 << "\" y=\"" << iy << "\">"
          << xml_escape(tr.name) << "</text>\n";
    }
    svg << "</g>\n";
  }
};

}  // namespace svg_detail
}  // namespace numc
