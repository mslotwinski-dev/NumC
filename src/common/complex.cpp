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

complex complex::operator*(const complex& other) const {
  return complex(x * other.x - y * other.y, x * other.y + y * other.x);
}

complex complex::operator/(const complex& other) const {
  return complex((x * other.x + y * other.y) / other.abs(), (y * other.x - x * other.y) / other.abs());
}

complex complex::operator!() const { return complex(x, -y); }

complex complex::operator^(const complex& other) const {
  double r = this->abs();
  double theta = this->arg();

  double new_r = std::pow(r, other.x) * std::exp(-other.y * theta);
  double new_theta = other.x * theta + other.y * std::log(r);

  return complex(new_r * std::cos(new_theta), new_r * std::sin(new_theta));
}

// Display

std::ostream& operator<<(std::ostream& os, const complex& c) {
  if (c.real() != 0) {
    os << c.real();
  }

  if (c.imag() > 0 && c.real() != 0) {
    os << " + ";
  } else if (c.imag() < 0 && c.real() != 0) {
    os << " - ";
  } else if (c.imag() < 0 && c.real() == 0) {
    os << "-";
  }

  if (c.imag() != 0) {
    if (std::abs(c.imag()) != 1) {
      os << std::abs(c.imag());
    }
    os << "i";
  }

  return os;
}

}  // namespace numc