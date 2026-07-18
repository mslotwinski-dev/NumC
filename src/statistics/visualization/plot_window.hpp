#pragma once

// ─────────────────────────────────────────────────────────────────────────────
// NumC Plot Window — Native Win32 + GDI+ rendering for matplotlib-quality plots
//
// On Windows: Creates a resizable, anti-aliased graphical window.
// On other platforms: Falls back to opening SVG in the default browser.
// ─────────────────────────────────────────────────────────────────────────────

#include "plot_renderer.hpp"

#if defined(_WIN32)

#ifndef NOMINMAX
#define NOMINMAX
#endif

#include <windows.h>

#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

#include <objidl.h>
#include <gdiplus.h>

#include <algorithm>
#include <cstdio>
#include <string>
#include <vector>

namespace numc {
namespace win_detail {

// ─────────────────────────────────────────────────────────────────────────────
// Utility: UTF-8 → UTF-16 conversion
// ─────────────────────────────────────────────────────────────────────────────

inline std::wstring to_wide(const std::string& str) {
  if (str.empty()) return L"";
  int size = MultiByteToWideChar(CP_UTF8, 0, str.c_str(), -1, nullptr, 0);
  if (size <= 0) return L"";
  std::wstring result(size - 1, 0);
  MultiByteToWideChar(CP_UTF8, 0, str.c_str(), -1, &result[0], size);
  return result;
}

// ─────────────────────────────────────────────────────────────────────────────
// Utility: hex color string → Gdiplus::Color
// ─────────────────────────────────────────────────────────────────────────────

inline Gdiplus::Color hex_to_color(const char* hex, BYTE alpha = 255) {
  unsigned int r = 0, g = 0, b = 0;
  if (hex && hex[0] == '#') {
    std::sscanf(hex + 1, "%02x%02x%02x", &r, &g, &b);
  }
  return Gdiplus::Color(alpha, (BYTE)r, (BYTE)g, (BYTE)b);
}

inline Gdiplus::Color hex_to_color(const std::string& hex, BYTE alpha = 255) {
  return hex_to_color(hex.c_str(), alpha);
}

// ─────────────────────────────────────────────────────────────────────────────
// PlotWindowData — all data needed for rendering, stored in GWLP_USERDATA
// ─────────────────────────────────────────────────────────────────────────────

struct PlotWindowData {
  std::string title;
  std::string xlabel;
  std::string ylabel;
  std::vector<svg_detail::TraceData<double>> traces;
  PlotTheme theme;
};

// ─────────────────────────────────────────────────────────────────────────────
// render_plot — the core GDI+ rendering routine
// ─────────────────────────────────────────────────────────────────────────────

inline void render_plot(Gdiplus::Graphics& g, int win_w, int win_h,
                        const PlotWindowData& data) {
  using namespace Gdiplus;
  using namespace svg_detail;

  g.SetSmoothingMode(SmoothingModeHighQuality);
  g.SetTextRenderingHint(TextRenderingHintClearTypeGridFit);
  g.SetPixelOffsetMode(PixelOffsetModeHighQuality);

  PlotStyle s = PlotStyle::from_theme(data.theme);
  s.width  = win_w;
  s.height = win_h;

  // Adjust margins for window size
  s.margin_left   = std::max(70, win_w / 14);
  s.margin_right  = std::max(20, win_w / 40);
  s.margin_top    = std::max(50, win_h / 12);
  s.margin_bottom = std::max(55, win_h / 10);

  int pw = s.plot_width();
  int ph = s.plot_height();
  int px = s.plot_x();
  int py = s.plot_y();

  if (pw < 50 || ph < 50) return;  // Window too small

  // ── Data ranges and ticks ──
  double x_min_d, x_max_d, y_min_d, y_max_d;
  compute_data_range(data.traces, x_min_d, x_max_d, y_min_d, y_max_d);

  auto x_ticks = compute_ticks(x_min_d, x_max_d, 7);
  auto y_ticks = compute_ticks(y_min_d, y_max_d, 7);

  double xmin = x_ticks.tick_min;
  double xmax = x_ticks.tick_max;
  double ymin = y_ticks.tick_min;
  double ymax = y_ticks.tick_max;

  auto map_x = [&](double val) -> float {
    return (float)(px + (val - xmin) / (xmax - xmin) * pw);
  };
  auto map_y = [&](double val) -> float {
    return (float)(py + ph - (val - ymin) / (ymax - ymin) * ph);
  };

  // ── Fonts ──
  FontFamily fontFamily(L"Segoe UI");
  Font titleFont(&fontFamily, (float)std::max(14, win_h / 40), FontStyleBold, UnitPixel);
  Font labelFont(&fontFamily, (float)std::max(11, win_h / 55), FontStyleRegular, UnitPixel);
  Font tickFont(&fontFamily, (float)std::max(10, win_h / 65), FontStyleRegular, UnitPixel);
  Font legendFont(&fontFamily, (float)std::max(10, win_h / 60), FontStyleRegular, UnitPixel);
  Font watermarkFont(&fontFamily, 9.0f, FontStyleRegular, UnitPixel);

  // ── Background ──
  SolidBrush bgBrush(hex_to_color(s.bg_color));
  g.FillRectangle(&bgBrush, 0, 0, win_w, win_h);

  // ── Plot area background (gradient) ──
  LinearGradientBrush plotBg(
    Point(px, py), Point(px, py + ph),
    hex_to_color(s.plot_bg_top), hex_to_color(s.plot_bg_bottom));
  g.FillRectangle(&plotBg, px, py, pw, ph);

  // ── Grid lines ──
  {
    Pen gridPen(hex_to_color(s.grid_color), (float)s.grid_width);
    REAL dashPattern[] = {6.0f, 4.0f};
    gridPen.SetDashPattern(dashPattern, 2);

    Pen zeroPen(hex_to_color(s.grid_zero_color), 1.2f);

    // Horizontal grid (Y ticks)
    for (double yt : y_ticks.ticks) {
      float sy = map_y(yt);
      if (sy >= py && sy <= py + ph) {
        bool is_zero = std::abs(yt) < 1e-12;
        g.DrawLine(is_zero ? &zeroPen : &gridPen,
                   (float)px, sy, (float)(px + pw), sy);
      }
    }

    // Vertical grid (X ticks)
    for (double xt : x_ticks.ticks) {
      float sx = map_x(xt);
      if (sx >= px && sx <= px + pw) {
        bool is_zero = std::abs(xt) < 1e-12;
        g.DrawLine(is_zero ? &zeroPen : &gridPen,
                   sx, (float)py, sx, (float)(py + ph));
      }
    }
  }

  // ── Axis border ──
  Pen axisPen(hex_to_color(s.axis_color), 1.2f);
  g.DrawRectangle(&axisPen, px, py, pw, ph);

  // ── Tick marks and labels ──
  {
    Pen tickPen(hex_to_color(s.axis_color), 1.2f);
    SolidBrush tickBrush(hex_to_color(s.tick_color));

    StringFormat sfCenter;
    sfCenter.SetAlignment(StringAlignmentCenter);
    sfCenter.SetLineAlignment(StringAlignmentNear);

    StringFormat sfRight;
    sfRight.SetAlignment(StringAlignmentFar);
    sfRight.SetLineAlignment(StringAlignmentCenter);

    // X ticks
    for (double xt : x_ticks.ticks) {
      float sx = map_x(xt);
      if (sx >= px - 1 && sx <= px + pw + 1) {
        g.DrawLine(&tickPen, sx, (float)(py + ph), sx, (float)(py + ph + 5));
        std::wstring label = to_wide(format_tick(xt));
        RectF labelRect(sx - 40, (float)(py + ph + 7), 80, 20);
        g.DrawString(label.c_str(), -1, &tickFont, labelRect, &sfCenter, &tickBrush);
      }
    }

    // Y ticks
    for (double yt : y_ticks.ticks) {
      float sy = map_y(yt);
      if (sy >= py - 1 && sy <= py + ph + 1) {
        g.DrawLine(&tickPen, (float)(px - 5), sy, (float)px, sy);
        std::wstring label = to_wide(format_tick(yt));
        RectF labelRect((float)(px - 70), sy - 8, 60, 16);
        g.DrawString(label.c_str(), -1, &tickFont, labelRect, &sfRight, &tickBrush);
      }
    }
  }

  // ── Data traces (clipped to plot area) ──
  {
    Region clipRegion(Rect(px, py, pw, ph));
    g.SetClip(&clipRegion);

    for (size_t i = 0; i < data.traces.size(); ++i) {
      const auto& tr = data.traces[i];
      Color traceColor = hex_to_color(
        tr.color.empty() ? s.trace_color(i) : tr.color.c_str());

      if (tr.mode == "lines") {
        // Build PointF array, splitting at NaN/Inf
        std::vector<PointF> segment;

        auto flush_segment = [&]() {
          if (segment.size() >= 2) {
            // Glow effect: draw wider semi-transparent line first
            Color glowColor(40, traceColor.GetRed(), traceColor.GetGreen(), traceColor.GetBlue());
            Pen glowPen(glowColor, (float)s.line_width + 5.0f);
            glowPen.SetLineCap(LineCapRound, LineCapRound, DashCapRound);
            glowPen.SetLineJoin(LineJoinRound);
            g.DrawLines(&glowPen, segment.data(), (int)segment.size());

            // Main line
            Pen linePen(traceColor, (float)s.line_width);
            linePen.SetLineCap(LineCapRound, LineCapRound, DashCapRound);
            linePen.SetLineJoin(LineJoinRound);
            g.DrawLines(&linePen, segment.data(), (int)segment.size());
          }
          segment.clear();
        };

        for (size_t j = 0; j < tr.x.size(); ++j) {
          double xv = tr.x[j];
          double yv = tr.y[j];
          if (std::isnan(xv) || std::isinf(xv) || std::isnan(yv) || std::isinf(yv)) {
            flush_segment();
            continue;
          }
          segment.push_back(PointF(map_x(xv), map_y(yv)));
        }
        flush_segment();

      } else if (tr.mode == "markers") {
        float r = (float)s.scatter_radius;

        for (size_t j = 0; j < tr.x.size(); ++j) {
          double xv = tr.x[j];
          double yv = tr.y[j];
          if (std::isnan(xv) || std::isinf(xv) || std::isnan(yv) || std::isinf(yv)) {
            continue;
          }
          float sx = map_x(xv);
          float sy = map_y(yv);

          // Outer glow
          Color glowColor(40, traceColor.GetRed(), traceColor.GetGreen(), traceColor.GetBlue());
          SolidBrush glowBrush(glowColor);
          g.FillEllipse(&glowBrush, sx - r * 1.5f, sy - r * 1.5f, r * 3, r * 3);

          // Outer ring
          Pen ringPen(Color(100, traceColor.GetRed(), traceColor.GetGreen(), traceColor.GetBlue()), 1.2f);
          g.DrawEllipse(&ringPen, sx - r - 1.0f, sy - r - 1.0f, (r + 1.0f) * 2, (r + 1.0f) * 2);

          // Inner solid circle with border
          SolidBrush fillBrush(traceColor);
          Pen borderPen(hex_to_color(s.bg_color), 1.8f);
          g.FillEllipse(&fillBrush, sx - r, sy - r, r * 2, r * 2);
          g.DrawEllipse(&borderPen, sx - r, sy - r, r * 2, r * 2);
        }
      }
    }

    g.ResetClip();
  }

  // ── Title ──
  {
    SolidBrush titleBrush(hex_to_color(s.title_color));
    std::wstring title_w = to_wide(data.title);
    StringFormat sf;
    sf.SetAlignment(StringAlignmentCenter);
    sf.SetLineAlignment(StringAlignmentFar);
    RectF titleRect(0, 0, (float)win_w, (float)(py - 6));
    g.DrawString(title_w.c_str(), -1, &titleFont, titleRect, &sf, &titleBrush);
  }

  // ── X axis label ──
  {
    SolidBrush labelBrush(hex_to_color(s.label_color));
    std::wstring xlabel_w = to_wide(data.xlabel);
    StringFormat sf;
    sf.SetAlignment(StringAlignmentCenter);
    sf.SetLineAlignment(StringAlignmentFar);
    RectF labelRect((float)px, (float)(py + ph + 28), (float)pw, 28);
    g.DrawString(xlabel_w.c_str(), -1, &labelFont, labelRect, &sf, &labelBrush);
  }

  // ── Y axis label (rotated) ──
  {
    SolidBrush labelBrush(hex_to_color(s.label_color));
    std::wstring ylabel_w = to_wide(data.ylabel);
    StringFormat sf;
    sf.SetAlignment(StringAlignmentCenter);
    sf.SetLineAlignment(StringAlignmentCenter);

    // Rotate around the label center
    float cx = 18.0f;
    float cy = (float)(py + ph / 2);

    g.TranslateTransform(cx, cy);
    g.RotateTransform(-90.0f);
    RectF labelRect(-200, -10, 400, 20);
    g.DrawString(ylabel_w.c_str(), -1, &labelFont, labelRect, &sf, &labelBrush);
    g.ResetTransform();
  }

  // ── Legend ──
  if (!data.traces.empty()) {
    int item_h = std::max(18, win_h / 35);
    int pad_x = 12;
    int pad_y = 8;
    int icon_w = 20;
    int gap = 8;

    // Measure text widths
    float max_text_w = 0;
    for (const auto& tr : data.traces) {
      std::wstring name_w = to_wide(tr.name);
      RectF measureRect;
      g.MeasureString(name_w.c_str(), -1, &legendFont, PointF(0, 0), &measureRect);
      if (measureRect.Width > max_text_w) max_text_w = measureRect.Width;
    }

    int box_w = (int)(pad_x * 2 + icon_w + gap + max_text_w + 8);
    int box_h = (int)(pad_y * 2 + data.traces.size() * item_h);
    int box_x = px + pw - box_w - 10;
    int box_y = py + 10;

    // Legend background with shadow
    SolidBrush shadowBrush(Color(50, 0, 0, 0));
    g.FillRectangle(&shadowBrush, box_x + 2, box_y + 3, box_w, box_h);

    SolidBrush legendBg(hex_to_color(s.legend_bg, 235));
    Pen legendBorder(hex_to_color(s.legend_border), 1.0f);
    g.FillRectangle(&legendBg, box_x, box_y, box_w, box_h);
    g.DrawRectangle(&legendBorder, box_x, box_y, box_w, box_h);

    SolidBrush legendTextBrush(hex_to_color(s.legend_text));
    StringFormat sfLeft;
    sfLeft.SetAlignment(StringAlignmentNear);
    sfLeft.SetLineAlignment(StringAlignmentCenter);

    for (size_t i = 0; i < data.traces.size(); ++i) {
      const auto& tr = data.traces[i];
      Color color = hex_to_color(
        tr.color.empty() ? s.trace_color(i) : tr.color.c_str());

      int iy = box_y + pad_y + (int)i * item_h + item_h / 2;

      if (tr.mode == "markers") {
        SolidBrush dotBrush(color);
        float dotR = 4.5f;
        g.FillEllipse(&dotBrush,
                      (float)(box_x + pad_x + icon_w / 2) - dotR,
                      (float)iy - dotR, dotR * 2, dotR * 2);
      } else {
        Pen linePen(color, 2.5f);
        linePen.SetLineCap(LineCapRound, LineCapRound, DashCapRound);
        g.DrawLine(&linePen,
                   (float)(box_x + pad_x), (float)iy,
                   (float)(box_x + pad_x + icon_w), (float)iy);
      }

      std::wstring name_w = to_wide(tr.name);
      RectF textRect((float)(box_x + pad_x + icon_w + gap),
                     (float)(iy - item_h / 2),
                     max_text_w + 8, (float)item_h);
      g.DrawString(name_w.c_str(), -1, &legendFont, textRect, &sfLeft, &legendTextBrush);
    }
  }

  // ── Watermark ──
  {
    SolidBrush wmBrush(hex_to_color(s.axis_color, 80));
    std::wstring wm = L"NumC";
    StringFormat sf;
    sf.SetAlignment(StringAlignmentFar);
    sf.SetLineAlignment(StringAlignmentFar);
    RectF wmRect(0, 0, (float)(win_w - 8), (float)(win_h - 4));
    g.DrawString(wm.c_str(), -1, &watermarkFont, wmRect, &sf, &wmBrush);
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// Window procedure
// ─────────────────────────────────────────────────────────────────────────────

inline LRESULT CALLBACK PlotWndProc(HWND hwnd, UINT msg, WPARAM wParam, LPARAM lParam) {
  switch (msg) {
    case WM_PAINT: {
      PAINTSTRUCT ps;
      HDC hdc = BeginPaint(hwnd, &ps);

      RECT rc;
      GetClientRect(hwnd, &rc);
      int w = rc.right - rc.left;
      int h = rc.bottom - rc.top;

      // Double buffering
      HDC bufferDC = CreateCompatibleDC(hdc);
      HBITMAP bufferBmp = CreateCompatibleBitmap(hdc, w, h);
      HGDIOBJ oldBmp = SelectObject(bufferDC, bufferBmp);

      {
        Gdiplus::Graphics g(bufferDC);
        PlotWindowData* data = reinterpret_cast<PlotWindowData*>(
            GetWindowLongPtr(hwnd, GWLP_USERDATA));
        if (data) {
          render_plot(g, w, h, *data);
        }
      }

      BitBlt(hdc, 0, 0, w, h, bufferDC, 0, 0, SRCCOPY);

      SelectObject(bufferDC, oldBmp);
      DeleteObject(bufferBmp);
      DeleteDC(bufferDC);
      EndPaint(hwnd, &ps);
      return 0;
    }

    case WM_SIZE:
      InvalidateRect(hwnd, nullptr, FALSE);
      return 0;

    case WM_KEYDOWN:
      if (wParam == VK_ESCAPE) {
        DestroyWindow(hwnd);
      }
      return 0;

    case WM_DESTROY:
      PostQuitMessage(0);
      return 0;

    default:
      return DefWindowProcW(hwnd, msg, wParam, lParam);
  }
}

// ─────────────────────────────────────────────────────────────────────────────
// show_window — creates the native window and runs the message loop
// ─────────────────────────────────────────────────────────────────────────────

inline void show_window(PlotWindowData& data) {
  // Initialize GDI+
  Gdiplus::GdiplusStartupInput gdiplusStartupInput;
  ULONG_PTR gdiplusToken;
  Gdiplus::GdiplusStartup(&gdiplusToken, &gdiplusStartupInput, nullptr);

  // Register window class
  static bool class_registered = false;
  static const wchar_t* CLASS_NAME = L"NumC_PlotWindow";

  if (!class_registered) {
    WNDCLASSEXW wc = {};
    wc.cbSize        = sizeof(WNDCLASSEXW);
    wc.lpfnWndProc   = PlotWndProc;
    wc.hInstance      = GetModuleHandle(nullptr);
    wc.hCursor       = LoadCursor(nullptr, IDC_ARROW);
    wc.hbrBackground = nullptr;
    wc.lpszClassName  = CLASS_NAME;
    wc.style         = CS_HREDRAW | CS_VREDRAW;

    // Try loading the icon from a few common relative paths
    HICON hIcon = (HICON)LoadImageW(nullptr, L"src/icon.ico", IMAGE_ICON, 0, 0, LR_LOADFROMFILE | LR_DEFAULTSIZE);
    if (!hIcon) {
      hIcon = (HICON)LoadImageW(nullptr, L"../../src/icon.ico", IMAGE_ICON, 0, 0, LR_LOADFROMFILE | LR_DEFAULTSIZE);
    }
    if (!hIcon) {
      hIcon = (HICON)LoadImageW(nullptr, L"../src/icon.ico", IMAGE_ICON, 0, 0, LR_LOADFROMFILE | LR_DEFAULTSIZE);
    }

    if (hIcon) {
      wc.hIcon = hIcon;
      wc.hIconSm = hIcon;
    } else {
      wc.hIcon = LoadIcon(nullptr, IDI_APPLICATION);
      wc.hIconSm = LoadIcon(nullptr, IDI_APPLICATION);
    }

    RegisterClassExW(&wc);
    class_registered = true;
  }

  // Calculate centered window position
  int screen_w = GetSystemMetrics(SM_CXSCREEN);
  int screen_h = GetSystemMetrics(SM_CYSCREEN);
  int win_w = std::min(1100, screen_w - 100);
  int win_h = std::min(680, screen_h - 100);
  int win_x = (screen_w - win_w) / 2;
  int win_y = (screen_h - win_h) / 2;

  // Window title
  std::wstring title_w = to_wide(data.title) + L" \u2014 NumC";

  // Create window
  HWND hwnd = CreateWindowExW(
    0,
    CLASS_NAME,
    title_w.c_str(),
    WS_OVERLAPPEDWINDOW,
    win_x, win_y, win_w, win_h,
    nullptr, nullptr, GetModuleHandle(nullptr), nullptr);

  if (!hwnd) {
    Gdiplus::GdiplusShutdown(gdiplusToken);
    return;
  }

  // Store plot data in window
  SetWindowLongPtr(hwnd, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(&data));

  ShowWindow(hwnd, SW_SHOW);
  UpdateWindow(hwnd);

  // Message loop (blocking until window is closed)
  MSG msg;
  while (GetMessage(&msg, nullptr, 0, 0)) {
    TranslateMessage(&msg);
    DispatchMessage(&msg);
  }

  Gdiplus::GdiplusShutdown(gdiplusToken);
}

}  // namespace win_detail
}  // namespace numc

#endif  // _WIN32
