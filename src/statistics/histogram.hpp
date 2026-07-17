#pragma once

#include "../common/vector.hpp"
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <iomanip>

namespace numc {
namespace statistics {

/// @addtogroup statistics
/// @{

/// @brief 1D Histogram class similar to ROOT's TH1D.
template <typename T = double>
class Histogram1D {
 private:
  size_t num_bins;
  T x_min, x_max;
  T bin_width;
  std::vector<T> bins;
  size_t underflow, overflow;

 public:
  Histogram1D(size_t n_bins, T xmin, T xmax) 
    : num_bins(n_bins), x_min(xmin), x_max(xmax), underflow(0), overflow(0) {
    if (n_bins == 0 || xmin >= xmax) {
      throw std::invalid_argument("Invalid histogram parameters");
    }
    bins.resize(n_bins, T(0.0));
    bin_width = (x_max - x_min) / T(num_bins);
  }

  Histogram1D(size_t n_bins, T xmin, T xmax, const vector<T>& samples)
    : Histogram1D(n_bins, xmin, xmax) {
    for (T val : samples) {
      fill(val);
    }
  }

  void fill(T value, T weight = T(1.0)) {
    if (value < x_min) underflow++;
    else if (value >= x_max) overflow++;
    else {
      size_t idx = static_cast<size_t>((value - x_min) / bin_width);
      if (idx >= num_bins) idx = num_bins - 1; // Safeguard against floating point imprecision
      bins[idx] += weight;
    }
  }

  void set_bin_content(size_t bin, T content) {
    if (bin < num_bins) bins[bin] = content;
  }

  T get_bin_content(size_t bin) const {
    if (bin < num_bins) return bins[bin];
    return T(0.0);
  }

  T get_bin_center(size_t bin) const {
    if (bin >= num_bins) return T(0.0);
    return x_min + (static_cast<T>(bin) + T(0.5)) * bin_width;
  }

  size_t get_num_bins() const { return num_bins; }

  vector<T> get_bin_centers() const {
    vector<T> res(num_bins);
    for (size_t i = 0; i < num_bins; ++i) res[i] = get_bin_center(i);
    return res;
  }

  vector<T> get_bin_contents() const {
    vector<T> res(num_bins);
    for (size_t i = 0; i < num_bins; ++i) res[i] = bins[i];
    return res;
  }

  T get_integral() const {
    T sum = T(0.0);
    for (T val : bins) sum += val;
    return sum * bin_width;
  }

  T get_entries() const {
    T sum = T(0.0);
    for (T val : bins) sum += val;
    return sum;
  }

  void scale(T factor) {
    for (T& val : bins) val *= factor;
  }

  /// @brief Normalizes the histogram such that the integral over all bins equals 1.
  void normalize() {
    T area = get_integral();
    if (area > T(0.0)) scale(T(1.0) / area);
  }

  T get_mean() const {
    T sum = T(0.0);
    T entries = get_entries();
    if (entries == T(0.0)) return T(0.0);
    for (size_t i = 0; i < num_bins; ++i) {
      sum += bins[i] * get_bin_center(i);
    }
    return sum / entries;
  }

  T get_variance() const {
    T sum_sq = T(0.0);
    T mean = get_mean();
    T entries = get_entries();
    if (entries <= T(1.0)) return T(0.0);
    for (size_t i = 0; i < num_bins; ++i) {
      T diff = get_bin_center(i) - mean;
      sum_sq += bins[i] * diff * diff;
    }
    return sum_sq / (entries - T(1.0)); // Sample variance based on binned data approximation
  }

  T get_std_dev() const {
    return std::sqrt(get_variance());
  }

  /// @brief Returns the Cumulative Distribution Function (CDF) as a new histogram.
  Histogram1D<T> get_cdf() const {
    Histogram1D<T> cdf(num_bins, x_min, x_max);
    T cumsum = T(0.0);
    for (size_t i = 0; i < num_bins; ++i) {
      cumsum += bins[i];
      cdf.set_bin_content(i, cumsum * bin_width); 
      // Technically ROOT's integral just sums the bins, but normalized CDF requires bin_width handling.
      // Usually CDF ends at 1.0. If histogram is normalized, cumsum * bin_width ends at 1.0.
    }
    return cdf;
  }

  void print() const {
    std::cout << "Histogram1D [" << x_min << ", " << x_max << "] with " << num_bins << " bins.\n";
    T max_val = T(0.0);
    for (T val : bins) max_val = std::max(max_val, val);
    
    for (size_t i = 0; i < num_bins; ++i) {
      std::cout << std::fixed << std::setprecision(2) << get_bin_center(i) << " | ";
      int stars = max_val > T(0.0) ? static_cast<int>((bins[i] / max_val) * 50.0) : 0;
      for (int s = 0; s < stars; ++s) std::cout << "*";
      std::cout << " (" << bins[i] << ")\n";
    }
    std::cout << "Underflow: " << underflow << ", Overflow: " << overflow << "\n";
  }

  // Operators
  Histogram1D operator+(const Histogram1D& other) const {
    if (num_bins != other.num_bins || x_min != other.x_min || x_max != other.x_max)
      throw std::invalid_argument("Histograms must have same binning to add");
    Histogram1D res(num_bins, x_min, x_max);
    for (size_t i = 0; i < num_bins; ++i) res.bins[i] = bins[i] + other.bins[i];
    return res;
  }

  Histogram1D operator-(const Histogram1D& other) const {
    if (num_bins != other.num_bins || x_min != other.x_min || x_max != other.x_max)
      throw std::invalid_argument("Histograms must have same binning to subtract");
    Histogram1D res(num_bins, x_min, x_max);
    for (size_t i = 0; i < num_bins; ++i) res.bins[i] = bins[i] - other.bins[i];
    return res;
  }

  Histogram1D operator*(const Histogram1D& other) const {
    if (num_bins != other.num_bins || x_min != other.x_min || x_max != other.x_max)
      throw std::invalid_argument("Histograms must have same binning to multiply");
    Histogram1D res(num_bins, x_min, x_max);
    for (size_t i = 0; i < num_bins; ++i) res.bins[i] = bins[i] * other.bins[i];
    return res;
  }

  Histogram1D operator/(const Histogram1D& other) const {
    if (num_bins != other.num_bins || x_min != other.x_min || x_max != other.x_max)
      throw std::invalid_argument("Histograms must have same binning to divide");
    Histogram1D res(num_bins, x_min, x_max);
    for (size_t i = 0; i < num_bins; ++i) {
      if (other.bins[i] != T(0.0)) res.bins[i] = bins[i] / other.bins[i];
    }
    return res;
  }
};

/// @brief 2D Histogram class similar to ROOT's TH2D.
template <typename T = double>
class Histogram2D {
 private:
  size_t nx, ny;
  T x_min, x_max, y_min, y_max;
  T bx, by;
  std::vector<std::vector<T>> bins;

 public:
  Histogram2D(size_t nx_bins, T xmin, T xmax, size_t ny_bins, T ymin, T ymax)
    : nx(nx_bins), ny(ny_bins), x_min(xmin), x_max(xmax), y_min(ymin), y_max(ymax) {
    if (nx == 0 || ny == 0 || xmin >= xmax || ymin >= ymax) {
      throw std::invalid_argument("Invalid Histogram2D parameters");
    }
    bins.assign(nx, std::vector<T>(ny, T(0.0)));
    bx = (x_max - x_min) / T(nx);
    by = (y_max - y_min) / T(ny);
  }

  Histogram2D(size_t nx_bins, T xmin, T xmax, size_t ny_bins, T ymin, T ymax, const std::pair<vector<T>, vector<T>>& samples)
    : Histogram2D(nx_bins, xmin, xmax, ny_bins, ymin, ymax) {
    size_t n = std::min(samples.first.size(), samples.second.size());
    for (size_t i = 0; i < n; ++i) {
      fill(samples.first[i], samples.second[i]);
    }
  }

  void fill(T x, T y, T weight = T(1.0)) {
    if (x >= x_min && x < x_max && y >= y_min && y < y_max) {
      size_t ix = static_cast<size_t>((x - x_min) / bx);
      size_t iy = static_cast<size_t>((y - y_min) / by);
      if (ix >= nx) ix = nx - 1;
      if (iy >= ny) iy = ny - 1;
      bins[ix][iy] += weight;
    }
  }

  T get_bin_center_x(size_t ix) const { return x_min + (static_cast<T>(ix) + T(0.5)) * bx; }
  T get_bin_center_y(size_t iy) const { return y_min + (static_cast<T>(iy) + T(0.5)) * by; }
  T get_bin_content(size_t ix, size_t iy) const { return bins[ix][iy]; }
  void set_bin_content(size_t ix, size_t iy, T content) { bins[ix][iy] = content; }
  size_t get_nx() const { return nx; }
  size_t get_ny() const { return ny; }

  T get_integral() const {
    T sum = T(0.0);
    for (size_t ix = 0; ix < nx; ++ix) {
      for (size_t iy = 0; iy < ny; ++iy) {
        sum += bins[ix][iy];
      }
    }
    return sum * bx * by;
  }

  T get_entries() const {
    T sum = T(0.0);
    for (size_t ix = 0; ix < nx; ++ix) {
      for (size_t iy = 0; iy < ny; ++iy) {
        sum += bins[ix][iy];
      }
    }
    return sum;
  }

  void scale(T factor) {
    for (size_t ix = 0; ix < nx; ++ix) {
      for (size_t iy = 0; iy < ny; ++iy) bins[ix][iy] *= factor;
    }
  }

  void normalize() {
    T area = get_integral();
    if (area > T(0.0)) scale(T(1.0) / area);
  }

  Histogram1D<T> projection_x() const {
    Histogram1D<T> h(nx, x_min, x_max);
    for (size_t ix = 0; ix < nx; ++ix) {
      T sum = T(0.0);
      for (size_t iy = 0; iy < ny; ++iy) {
        sum += bins[ix][iy];
      }
      h.set_bin_content(ix, sum * by); 
    }
    return h;
  }

  Histogram1D<T> projection_y() const {
    Histogram1D<T> h(ny, y_min, y_max);
    for (size_t iy = 0; iy < ny; ++iy) {
      T sum = T(0.0);
      for (size_t ix = 0; ix < nx; ++ix) {
        sum += bins[ix][iy];
      }
      h.set_bin_content(iy, sum * bx);
    }
    return h;
  }

  Histogram2D<T> get_cdf() const {
    Histogram2D<T> cdf(nx, x_min, x_max, ny, y_min, y_max);
    T cumsum = T(0.0);
    // Standard 2D CDF builds up from bottom left to top right
    for (size_t ix = 0; ix < nx; ++ix) {
      for (size_t iy = 0; iy < ny; ++iy) {
        T sum = T(0.0);
        for(size_t i=0; i<=ix; ++i) {
          for(size_t j=0; j<=iy; ++j) {
            sum += bins[i][j];
          }
        }
        cdf.set_bin_content(ix, iy, sum * bx * by);
      }
    }
    return cdf;
  }

  T mean_x() const {
    T sum = T(0.0);
    T entries = get_entries();
    if (entries == T(0.0)) return T(0.0);
    for (size_t ix = 0; ix < nx; ++ix) {
      for (size_t iy = 0; iy < ny; ++iy) {
        sum += bins[ix][iy] * get_bin_center_x(ix);
      }
    }
    return sum / entries;
  }

  T mean_y() const {
    T sum = T(0.0);
    T entries = get_entries();
    if (entries == T(0.0)) return T(0.0);
    for (size_t ix = 0; ix < nx; ++ix) {
      for (size_t iy = 0; iy < ny; ++iy) {
        sum += bins[ix][iy] * get_bin_center_y(iy);
      }
    }
    return sum / entries;
  }

  T variance_x() const {
    T sum_sq = T(0.0);
    T mx = mean_x();
    T entries = get_entries();
    if (entries <= T(1.0)) return T(0.0);
    for (size_t ix = 0; ix < nx; ++ix) {
      for (size_t iy = 0; iy < ny; ++iy) {
        T dx = get_bin_center_x(ix) - mx;
        sum_sq += bins[ix][iy] * dx * dx;
      }
    }
    return sum_sq / (entries - T(1.0));
  }

  T variance_y() const {
    T sum_sq = T(0.0);
    T my = mean_y();
    T entries = get_entries();
    if (entries <= T(1.0)) return T(0.0);
    for (size_t ix = 0; ix < nx; ++ix) {
      for (size_t iy = 0; iy < ny; ++iy) {
        T dy = get_bin_center_y(iy) - my;
        sum_sq += bins[ix][iy] * dy * dy;
      }
    }
    return sum_sq / (entries - T(1.0));
  }

  T covariance() const {
    T sum_xy = T(0.0);
    T mx = mean_x(), my = mean_y();
    T entries = get_entries();
    if (entries <= T(1.0)) return T(0.0);
    for (size_t ix = 0; ix < nx; ++ix) {
      for (size_t iy = 0; iy < ny; ++iy) {
        T dx = get_bin_center_x(ix) - mx;
        T dy = get_bin_center_y(iy) - my;
        sum_xy += bins[ix][iy] * dx * dy;
      }
    }
    return sum_xy / (entries - T(1.0));
  }

  T correlation() const {
    T vx = variance_x();
    T vy = variance_y();
    if (vx == T(0.0) || vy == T(0.0)) return T(0.0);
    return covariance() / std::sqrt(vx * vy);
  }

  // Operators
  Histogram2D operator+(const Histogram2D& other) const {
    if (nx != other.nx || ny != other.ny || x_min != other.x_min || x_max != other.x_max || y_min != other.y_min || y_max != other.y_max)
      throw std::invalid_argument("Histograms must have same binning to add");
    Histogram2D res(nx, x_min, x_max, ny, y_min, y_max);
    for (size_t i = 0; i < nx; ++i)
      for (size_t j = 0; j < ny; ++j) res.bins[i][j] = bins[i][j] + other.bins[i][j];
    return res;
  }

  Histogram2D operator-(const Histogram2D& other) const {
    if (nx != other.nx || ny != other.ny || x_min != other.x_min || x_max != other.x_max || y_min != other.y_min || y_max != other.y_max)
      throw std::invalid_argument("Histograms must have same binning to subtract");
    Histogram2D res(nx, x_min, x_max, ny, y_min, y_max);
    for (size_t i = 0; i < nx; ++i)
      for (size_t j = 0; j < ny; ++j) res.bins[i][j] = bins[i][j] - other.bins[i][j];
    return res;
  }

  Histogram2D operator*(const Histogram2D& other) const {
    if (nx != other.nx || ny != other.ny || x_min != other.x_min || x_max != other.x_max || y_min != other.y_min || y_max != other.y_max)
      throw std::invalid_argument("Histograms must have same binning to multiply");
    Histogram2D res(nx, x_min, x_max, ny, y_min, y_max);
    for (size_t i = 0; i < nx; ++i)
      for (size_t j = 0; j < ny; ++j) res.bins[i][j] = bins[i][j] * other.bins[i][j];
    return res;
  }

  Histogram2D operator/(const Histogram2D& other) const {
    if (nx != other.nx || ny != other.ny || x_min != other.x_min || x_max != other.x_max || y_min != other.y_min || y_max != other.y_max)
      throw std::invalid_argument("Histograms must have same binning to divide");
    Histogram2D res(nx, x_min, x_max, ny, y_min, y_max);
    for (size_t i = 0; i < nx; ++i)
      for (size_t j = 0; j < ny; ++j) {
        if (other.bins[i][j] != T(0.0)) res.bins[i][j] = bins[i][j] / other.bins[i][j];
      }
    return res;
  }
};

/// @}

}  // namespace statistics
}  // namespace numc
