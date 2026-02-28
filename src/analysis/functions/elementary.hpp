#pragma once

#include "../../common/function.hpp"
#include "../../inc.hpp"
#include "../../utility/log.hpp"

namespace numc {

template <typename T>
numc::func<T> exp(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::exp(f(x)); });
}

template <typename T>
numc::func<T> ln(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::log(f(x)); });
}

template <typename T>
numc::func<T> log2(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::log2(f(x)); });
}

template <typename T>
numc::func<T> log10(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::log10(f(x)); });
}

template <typename T>
numc::func<T> log(T base, const numc::func<T>& f) {
  return numc::func<T>([base, f](T x) { return std::log(f(x)) / std::log(base); });
}

template <typename T>
numc::func<T> sin(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::sin(f(x)); });
}

template <typename T>
numc::func<T> cos(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::cos(f(x)); });
}

template <typename T>
numc::func<T> tan(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::tan(f(x)); });
}

template <typename T>
numc::func<T> cot(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return T(1.0) / std::tan(f(x)); });
}

template <typename T>
numc::func<T> asin(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::asin(f(x)); });
}

template <typename T>
numc::func<T> acos(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::acos(f(x)); });
}

template <typename T>
numc::func<T> atan(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::atan(f(x)); });
}

template <typename T>
numc::func<T> acot(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::atan2(T(1.0), f(x)); });
}

template <typename T>
numc::func<T> sinh(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::sinh(f(x)); });
}

template <typename T>
numc::func<T> cosh(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::cosh(f(x)); });
}

template <typename T>
numc::func<T> tanh(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::tanh(f(x)); });
}

template <typename T>
numc::func<T> coth(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return T(1.0) / std::tanh(f(x)); });
}

template <typename T>
numc::func<T> asinh(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::asinh(f(x)); });
}

template <typename T>
numc::func<T> acosh(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::acosh(f(x)); });
}

template <typename T>
numc::func<T> atanh(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::atanh(f(x)); });
}

template <typename T>
numc::func<T> acoth(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::atanh(T(1.0) / f(x)); });
}

template <typename T>
numc::func<T> sqrt(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::sqrt(f(x)); });
}

template <typename T>
numc::func<T> cbrt(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::cbrt(f(x)); });
}

template <typename T>
numc::func<T> root(T degree, const numc::func<T>& f) {
  return numc::func<T>([degree, f](T x) { return std::pow(f(x), T(1.0) / degree); });
}

template <typename T>
numc::func<T> abs(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::abs(f(x)); });
}

template <typename T>
numc::func<T> sign(const numc::func<T>& f) {
  return numc::func<T>([f](T x) {
    T val = f(x);
    if (val > T(0)) return T(1);
    if (val < T(0)) return T(-1);
    return T(0);
  });
}

template <typename T>
numc::func<T> floor(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::floor(f(x)); });
}

template <typename T>
numc::func<T> ceil(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::ceil(f(x)); });
}

template <typename T>
numc::func<T> round(const numc::func<T>& f, int places = 0) {
  return numc::func<T>([f, places](T x) {
    T factor = std::pow(T(10.0), T(places));
    return std::round(f(x) * factor) / factor;
  });
}

}  // namespace numc