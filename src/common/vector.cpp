#include "vector.hpp"

namespace numc {

// Structural

vector::vector() : data(nullptr), _size(0) {}

vector::vector(size_t size) : data(new double[size]()), _size(size) {}

vector::vector(std::initializer_list<double> values) : data(new double[values.size()]), _size(values.size()) {
  std::copy(values.begin(), values.end(), data);
}

vector::~vector() { delete[] data; }

double& vector::operator[](size_t index) { return this->data[index]; };

const double& vector::operator[](size_t index) const { return this->data[index]; };

// Arithmetic

vector vector::operator+(const vector& other) const {
  if (this->size() != other.size()) {
    Log::Warn("Addition requires vectors of the same size. Filling missing values with 0.");
  }

  vector blank(this->size() > other.size() ? this->size() : other.size());
  for (size_t i = 0; i < blank.size(); ++i) {
    blank[i] = (*this)[i] + other[i];
  }
  return blank;
}

vector vector::operator-(const vector& other) const {
  if (this->size() != other.size()) {
    Log::Warn("Subtraction requires vectors of the same size. Filling missing values with 0.");
  }

  vector blank(this->size() > other.size() ? this->size() : other.size());
  for (size_t i = 0; i < blank.size(); ++i) {
    blank[i] = (*this)[i] - other[i];
  }
  return blank;
}

double vector::operator*(const vector& other) const {
  if (this->size() != other.size()) {
    Log::Error("Dot product requires vectors of the same size.");
    return 0.0;
  }

  double result = 0.0;
  for (size_t i = 0; i < this->size(); ++i) {
    result += (*this)[i] * other[i];
  }
  return result;
}

vector vector::operator^(const vector& other) const {
  if (this->size() != 3 || other.size() != 3) {
    Log::Error("Cross product requires two vectors of size 3.");
    return vector();
  }

  double result_1 = (*this)[1] * other[2] - (*this)[2] * other[1];
  double result_2 = (*this)[2] * other[0] - (*this)[0] * other[2];
  double result_3 = (*this)[0] * other[1] - (*this)[1] * other[0];

  return {result_1, result_2, result_3};
}

// Display

std::ostream& operator<<(std::ostream& os, const vector& v) {
  os << "(";
  for (size_t i = 0; i < v.size(); ++i) {
    os << v.data[i];
    if (i + 1 < v.size()) os << ", ";
  }
  os << ")";
  return os;
}

}  // namespace numc
