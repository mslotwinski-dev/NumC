#pragma once

#include "../../common/constants.hpp"
#include "../../common/function.hpp"
#include "../../inc.hpp"
#include "../../utility/log.hpp"

namespace numc {

template <typename T>
numc::func<T> heaviside(const numc::func<T>& f) {
  return numc::func<T>([f](T x) {
    T val = f(x);
    if (val > T(0.0)) return T(1.0);
    if (val < T(0.0)) return T(0.0);
    return T(0.5);
  });
}

template <typename T>
numc::func<T> dirac(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return (f(x) == T(0.0)) ? T(1.0) : T(0.0); });
}

template <typename T>
numc::func<T> ramp(const numc::func<T>& f) {
  return numc::func<T>([f](T x) {
    T val = f(x);
    return (val > T(0.0)) ? val : T(0.0);
  });
}

template <typename T>
numc::func<T> sinc(const numc::func<T>& f) {
  return numc::func<T>([f](T x) {
    T val = f(x);
    if (std::abs(val) < T(1e-9)) {
      return T(1.0);
    }
    return std::sin(val) / val;
  });
}

template <typename T>
numc::func<T> rect(const numc::func<T>& f) {
  return numc::func<T>([f](T x) {
    T abs_val = std::abs(f(x));
    if (abs_val < T(0.5)) return T(1.0);
    if (abs_val == T(0.5)) return T(0.5);
    return T(0.0);
  });
}

template <typename T>
numc::func<T> tri(const numc::func<T>& f) {
  return numc::func<T>([f](T x) {
    T abs_val = std::abs(f(x));
    return (abs_val < T(1.0)) ? (T(1.0) - abs_val) : T(0.0);
  });
}

template <typename T>
numc::func<T> sawtooth(const numc::func<T>& f) {
  return numc::func<T>([f](T x) {
    T val = f(x);
    T phase = val / (T(2.0) * T(numc::PI));
    return T(2.0) * (phase - std::floor(phase + T(0.5)));
  });
}

template <typename T>
numc::func<T> square(const numc::func<T>& f, T duty_cycle = T(0.5)) {
  T safe_duty = duty_cycle;
  if (safe_duty < T(0.0)) safe_duty = T(0.0);
  if (safe_duty > T(1.0)) safe_duty = T(1.0);

  return numc::func<T>([f, safe_duty](T x) {
    T val = f(x);
    T phase = val / (T(2.0) * T(numc::PI));
    phase = phase - std::floor(phase);

    return (phase < safe_duty) ? T(1.0) : T(-1.0);
  });
}

}  // namespace numc