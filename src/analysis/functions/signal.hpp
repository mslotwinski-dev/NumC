#pragma once

#include <cmath>

#include "../../common/constants.hpp"
#include "../../common/function.hpp"
#include "../../inc.hpp"
#include "../../utility/log.hpp"

namespace numc {

/// @defgroup signal Signal Processing
/// @ingroup functions
/// @brief Step functions, impulses, windowing functions, and periodic waveform generators.
/// @{

/// @brief Computes the Heaviside step function, $H(x)$.
/// Returns 1.0 for strictly positive values, 0.0 for strictly negative values, and 0.5 at exactly 0.0.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing $H(f(x))$.
template <typename T>
numc::func<T> heaviside(const numc::func<T>& f) {
  return numc::func<T>([f](T x) {
    T val = f(x);
    if (val > T(0.0)) return T(1.0);
    if (val < T(0.0)) return T(0.0);
    return T(0.5);
  });
}

/// @brief Computes the discrete Dirac delta (impulse) function, $\delta(x)$.
/// Returns 1.0 if the input is exactly 0.0, and 0.0 everywhere else.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing $\delta(f(x))$.
template <typename T>
numc::func<T> dirac(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return (f(x) == T(0.0)) ? T(1.0) : T(0.0); });
}

/// @brief Computes the Ramp function, $R(x) = \max(0, x)$.
/// Linearly increases for positive values and remains 0.0 for negative values.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing the ramp function applied to $f(x)$.
template <typename T>
numc::func<T> ramp(const numc::func<T>& f) {
  return numc::func<T>([f](T x) {
    T val = f(x);
    return (val > T(0.0)) ? val : T(0.0);
  });
}

/// @brief Computes the unnormalized sinc function, $\operatorname{sinc}(x) = \sin(x)/x$.
/// Handles the singularity at $x = 0$ safely by returning 1.0 for values extremely close to zero.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing $\operatorname{sinc}(f(x))$.
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

/// @brief Computes the rectangular (boxcar) function, $\operatorname{rect}(x)$.
/// Returns 1.0 inside the interval (-0.5, 0.5), 0.5 at the boundaries ($\pm0.5$), and 0.0 outside.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing $\operatorname{rect}(f(x))$.
template <typename T>
numc::func<T> rect(const numc::func<T>& f) {
  return numc::func<T>([f](T x) {
    T abs_val = std::abs(f(x));
    if (abs_val < T(0.5)) return T(1.0);
    if (abs_val == T(0.5)) return T(0.5);
    return T(0.0);
  });
}

/// @brief Computes the triangular function, $\operatorname{tri}(x)$.
/// Returns $1 - |x|$ for $|x| < 1$, and 0.0 otherwise.
/// @tparam T Floating-point type.
/// @param f The input expression.
/// @return An expression representing $\operatorname{tri}(f(x))$.
template <typename T>
numc::func<T> tri(const numc::func<T>& f) {
  return numc::func<T>([f](T x) {
    T abs_val = std::abs(f(x));
    return (abs_val < T(1.0)) ? (T(1.0) - abs_val) : T(0.0);
  });
}

/// @brief Generates a periodic sawtooth wave with a period of $2\pi$.
/// The amplitude ranges from -1.0 to 1.0.
/// @tparam T Floating-point type.
/// @param f The input expression (usually representing time or phase).
/// @return An expression representing the sawtooth wave evaluated at $f(x)$.
template <typename T>
numc::func<T> sawtooth(const numc::func<T>& f) {
  return numc::func<T>([f](T x) {
    T val = f(x);
    T phase = val / (T(2.0) * T(numc::PI));
    return T(2.0) * (phase - std::floor(phase + T(0.5)));
  });
}

/// @brief Generates a periodic square wave with a period of $2\pi$ and adjustable duty cycle.
/// The amplitude alternates between 1.0 and -1.0.
/// @tparam T Floating-point type.
/// @param f The input expression (usually representing time or phase).
/// @param duty_cycle The fraction of the period where the signal is positive (clamped between 0.0 and 1.0). Default is 0.5 (50%).
/// @return An expression representing the square wave evaluated at $f(x)$.
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

/// @} // End of signal group

}  // namespace numc