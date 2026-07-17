#pragma once

#include "../../common/function.hpp"
#include "../../common/vector.hpp"

namespace numc {
namespace analysis {

/// @brief Represents a mathematical vector field: a vector of functions.
/// Inherits from std::vector<func<T>> but adds an evaluation operator 
/// that evaluates all functions at a given point and returns a numc::vector<T>.
template <typename T = double>
class vector_field : public std::vector<func<T>> {
public:
  using std::vector<func<T>>::vector;
  
  vector_field(const std::vector<func<T>>& base) : std::vector<func<T>>(base) {}

  /// @brief Evaluates the vector field at a given ND point.
  vector<T> operator()(const vector<T>& pt) const {
    vector<T> res(this->size());
    for (size_t i = 0; i < this->size(); ++i) {
      res[i] = (*this)[i](pt);
    }
    return res;
  }
};

/// @brief The Nabla (Del) operator.
/// Computes gradient, divergence, and curl of functions and vector fields.
template <typename T = double>
class nabla_op {
public:
  /// @brief Gradient of a scalar field. Returns a 3D vector field.
  vector_field<T> operator()(const func<T>& f) const {
    std::vector<func<T>> grad(3);
    for (size_t i = 0; i < 3; ++i) {
      grad[i] = func<T>([f, i](const vector<T>& v) {
        vector<T> v_plus = v;
        if (v_plus.size() <= i) v_plus.resize(i + 1, T(0.0));
        v_plus[i] += T(1e-7);

        vector<T> v_minus = v;
        if (v_minus.size() <= i) v_minus.resize(i + 1, T(0.0));
        v_minus[i] -= T(1e-7);

        return (f(v_plus) - f(v_minus)) / T(2e-7);
      });
    }
    return vector_field<T>(grad);
  }

  func<T> operator*(const std::vector<func<T>>& F) const {
    return func<T>([F](const vector<T>& v) {
      T div = T(0.0);
      for (size_t i = 0; i < F.size(); ++i) {
        vector<T> v_plus = v;
        if (v_plus.size() <= i) v_plus.resize(i + 1, T(0.0));
        v_plus[i] += T(1e-7);

        vector<T> v_minus = v;
        if (v_minus.size() <= i) v_minus.resize(i + 1, T(0.0));
        v_minus[i] -= T(1e-7);

        div += (F[i](v_plus) - F[i](v_minus)) / T(2e-7);
      }
      return div;
    });
  }

  /// @brief Curl (Rot) of a vector field (nabla ^ F). Returns a 3D vector field.
  vector_field<T> operator^(const std::vector<func<T>>& F) const {
    std::vector<func<T>> rot(3);
    if (F.size() < 3) {
      Log::Warn("Curl is only defined for 3D vector fields. Missing dimensions are treated as 0.");
    }
    
    // (dFz/dy - dFy/dz)
    rot[0] = func<T>([F](const vector<T>& v) {
      if (F.size() < 3) return T(0.0);
      vector<T> vpy = v, vmy = v, vpz = v, vmz = v;
      if (vpy.size() <= 1) vpy.resize(2, T(0.0)); if (vmy.size() <= 1) vmy.resize(2, T(0.0));
      if (vpz.size() <= 2) vpz.resize(3, T(0.0)); if (vmz.size() <= 2) vmz.resize(3, T(0.0));
      vpy[1] += T(1e-7); vmy[1] -= T(1e-7); 
      vpz[2] += T(1e-7); vmz[2] -= T(1e-7);
      return (F[2](vpy) - F[2](vmy)) / T(2e-7) - (F[1](vpz) - F[1](vmz)) / T(2e-7);
    });

    // (dFx/dz - dFz/dx)
    rot[1] = func<T>([F](const vector<T>& v) {
      if (F.size() < 3) return T(0.0);
      vector<T> vpx = v, vmx = v, vpz = v, vmz = v;
      if (vpx.size() <= 0) vpx.resize(1, T(0.0)); if (vmx.size() <= 0) vmx.resize(1, T(0.0));
      if (vpz.size() <= 2) vpz.resize(3, T(0.0)); if (vmz.size() <= 2) vmz.resize(3, T(0.0));
      vpz[2] += T(1e-7); vmz[2] -= T(1e-7); 
      vpx[0] += T(1e-7); vmx[0] -= T(1e-7);
      return (F[0](vpz) - F[0](vmz)) / T(2e-7) - (F[2](vpx) - F[2](vmx)) / T(2e-7);
    });

    // (dFy/dx - dFx/dy)
    rot[2] = func<T>([F](const vector<T>& v) {
      if (F.size() < 2) return T(0.0);
      vector<T> vpx = v, vmx = v, vpy = v, vmy = v;
      if (vpx.size() <= 0) vpx.resize(1, T(0.0)); if (vmx.size() <= 0) vmx.resize(1, T(0.0));
      if (vpy.size() <= 1) vpy.resize(2, T(0.0)); if (vmy.size() <= 1) vmy.resize(2, T(0.0));
      vpx[0] += T(1e-7); vmx[0] -= T(1e-7); 
      vpy[1] += T(1e-7); vmy[1] -= T(1e-7);
      return (F[1](vpx) - F[1](vmx)) / T(2e-7) - (F[0](vpy) - F[0](vmy)) / T(2e-7);
    });

    return vector_field<T>(rot);
  }
};

// Global nabla operator
inline const nabla_op<double> nabla;

} // namespace analysis

// Provide it globally in numc as well
using analysis::nabla;
using analysis::vector_field;

} // namespace numc
