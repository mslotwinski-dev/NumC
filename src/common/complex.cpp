#include "complex.hpp"

namespace numc {

// Properties

complex::complex(double real, double imag) : x(real), y(imag) {}

complex::~complex() {}

double complex::abs() const { return std::sqrt(x * x + y * y); }

double complex::arg() const { return std::atan2(y, x); }

complex complex::round(int places) {
  double factor = std::pow(10.0, places);
  x = std::round(x * factor) / factor;
  y = std::round(y * factor) / factor;

  return *this;
}

// Arithmetics

complex complex::operator+(const complex& other) const { return complex(x + other.x, y + other.y); }
complex complex::operator-(const complex& other) const { return complex(x - other.x, y - other.y); }
complex complex::operator*(const complex& other) const { return complex(x * other.x - y * other.y, x * other.y + y * other.x); }
complex complex::operator/(const complex& other) const {
  return complex((x * other.x + y * other.y) / other.abs(), (y * other.x - x * other.y) / other.abs());
}

complex& complex::operator+=(const complex& other) {
  x += other.x;
  y += other.y;
  return *this;
}

complex& complex::operator-=(const complex& other) {
  x -= other.x;
  y -= other.y;
  return *this;
}

complex& complex::operator*=(const complex& other) {
  double new_x = x * other.x - y * other.y;
  double new_y = x * other.y + y * other.x;
  x = new_x;
  y = new_y;
  return *this;
}

complex& complex::operator/=(const complex& other) {
  double denom = other.abs();
  x = (x * other.x + y * other.y) / denom;
  y = (y * other.x - x * other.y) / denom;
  return *this;
}

complex complex::operator~() const { return complex(x, -y); }

complex complex::operator^(const complex& other) const {
  double r = this->abs();
  double theta = this->arg();

  double new_r = std::pow(r, other.x) * std::exp(-other.y * theta);
  double new_theta = other.x * theta + other.y * std::log(r);

  return complex(new_r * std::cos(new_theta), new_r * std::sin(new_theta));
}

complex complex::conjugate() const { return ~(*this); }
complex complex::pow(const complex& other) const { return (*this) ^ other; }
std::vector<complex> complex::roots(int n) const {
  if (n <= 0) {
    Log::Error("Number of roots must be positive.");
    return {};
  }

  std::vector<complex> result;
  result.reserve(n);

  double r = this->abs();
  double theta = this->arg();

  double root_r = std::pow(r, 1.0 / n);
  for (int k = 0; k < n; ++k) {
    double root_theta = (theta + 2 * M_PI * k) / n;
    result[k] = complex(root_r * std::cos(root_theta), root_r * std::sin(root_theta));
  }

  return result;
}

// Logic

bool complex::operator==(const complex& other) const { return x == other.x && y == other.y; }
bool complex::operator!=(const complex& other) const { return !(*this == other); }

bool complex::operator<(const complex& other) const {
  Log::Warn(
      "Comparison of complex numbers is not well-defined. Comparing by magnitude. If you want to compare by real or imaginary part, use Re() or Im() "
      "methods.");
  return this->abs() < other.abs();
}

bool complex::operator>(const complex& other) const {
  Log::Warn(
      "Comparison of complex numbers is not well-defined. Comparing by magnitude. If you want to compare by real or imaginary part, use real() or "
      "imag() methods or Re() or Im() functions.");
  return this->arg() > other.arg();
}

bool complex::operator<=(const complex& other) const {
  if (*this == other) return true;
  return *this < other;
}

bool complex::operator>=(const complex& other) const {
  if (*this == other) return true;
  return *this > other;
}

// Display

std::ostream& operator<<(std::ostream& os, const complex& c) {
  double re = (std::abs(c.x) < 1e-14) ? 0.0 : c.x;
  double im = (std::abs(c.y) < 1e-14) ? 0.0 : c.y;

  if (re != 0) {
    os << re;
  }

  if (im > 0 && re != 0) {
    os << " + ";
  } else if (im < 0 && re != 0) {
    os << " - ";
  } else if (im < 0 && re == 0) {
    os << "-";
  }

  if (im != 0) {
    if (std::abs(im) != 1) {
      os << std::abs(im);
    }
    os << "i";
  }

  return os;
}

// Functions

double Re(const complex& c) { return c.real(); }
double Im(const complex& c) { return c.imag(); }

}  // namespace numc