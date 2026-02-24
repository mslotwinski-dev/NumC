#pragma once

#include "../inc.hpp"
#include "../utility/log.hpp"
#include "vector.hpp"

namespace numc {

class Tensor {
  Vector data;

 public:
  Tensor();
  ~Tensor();
};

}  // namespace numc