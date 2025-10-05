  #include "vector.hpp"

  namespace numc {

  vector::vector() : data(nullptr), size(0) {}

  vector::vector(size_t size) : data(new double[size]()), size(size) {}

  vector::vector(std::initializer_list<double> values) : data(new double[values.size()]), size(values.size()) {
    std::copy(values.begin(), values.end(), data);
  }

  vector::~vector() { delete[] data; }

  std::ostream& operator<<(std::ostream& os, const vector& v) {
    os << "(";
    for (size_t i = 0; i < v.size; ++i) {
      os << v.data[i];
      if (i + 1 < v.size) os << ", ";
    }
    os << ")";
    return os;
  }

  }  // namespace numc
