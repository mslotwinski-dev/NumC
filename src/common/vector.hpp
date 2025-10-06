#pragma once

#include "../inc.hpp"

namespace numc {

class vector {
  double* data;
  size_t size;

 public:
  vector();

  explicit vector(size_t size);

  vector(std::initializer_list<double> values);

  //  // vector(const vector& other);
  //  // vector(vector&& other) noexcept;
  //  // vector& operator=(const vector& other);
  //  // vector& operator=(vector&& other) noexcept;

  ~vector();

  vector operator+(const vector& other) const;
  vector operator-(const vector& other) const;
  double operator*(const vector& other) const;
  vector operator^(const vector& other) const;

  vector operator+=(const vector& other) const;
  vector operator-=(const vector& other) const;
  vector operator^=(const vector& other) const;

  // vector operator+(double scalar) const;
  // vector operator-(double scalar) const;

  vector operator*(double scalar) const;
  vector operator/(double scalar) const;

  double& operator[](size_t index);
  const double& operator[](size_t index) const;

  inline size_t size() const { return size; }

  // double norm() const;

  // vector round(int places = 0) const;

  friend std::ostream& operator<<(std::ostream& os, const vector& v);
};

std::ostream& operator<<(std::ostream& os, const vector& v);

}  // namespace numc
