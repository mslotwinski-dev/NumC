#pragma once

#include "../../common/function.hpp"
#include "../../inc.hpp"
#include "../../utility/log.hpp"

namespace numc {

template <typename T>
numc::func<T> erf(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::erf(f(x)); });
}

template <typename T>
numc::func<T> erfc(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::erfc(f(x)); });
}

template <typename T>
numc::func<T> tgamma(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::tgamma(f(x)); });
}

template <typename T>
numc::func<T> lgamma(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::lgamma(f(x)); });
}

template <typename T>
numc::func<T> beta(const numc::func<T>& f1, const numc::func<T>& f2) {
  return numc::func<T>([f1, f2](T x) { return std::beta(f1(x), f2(x)); });
}

template <typename T>
numc::func<T> comp_ellint_1(const numc::func<T>& k) {
  return numc::func<T>([k](T x) { return std::comp_ellint_1(k(x)); });
}

template <typename T>
numc::func<T> comp_ellint_2(const numc::func<T>& k) {
  return numc::func<T>([k](T x) { return std::comp_ellint_2(k(x)); });
}

template <typename T>
numc::func<T> comp_ellint_3(T nu, const numc::func<T>& k) {
  return numc::func<T>([nu, k](T x) { return std::comp_ellint_3(nu, k(x)); });
}

template <typename T>
numc::func<T> ellint_1(const numc::func<T>& k, const numc::func<T>& phi) {
  return numc::func<T>([k, phi](T x) { return std::ellint_1(k(x), phi(x)); });
}

template <typename T>
numc::func<T> ellint_2(const numc::func<T>& k, const numc::func<T>& phi) {
  return numc::func<T>([k, phi](T x) { return std::ellint_2(k(x), phi(x)); });
}

template <typename T>
numc::func<T> ellint_3(T nu, const numc::func<T>& k, const numc::func<T>& phi) {
  return numc::func<T>([nu, k, phi](T x) { return std::ellint_3(nu, k(x), phi(x)); });
}

template <typename T>
numc::func<T> expint(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::expint(f(x)); });
}

template <typename T>
numc::func<T> riemann_zeta(const numc::func<T>& f) {
  return numc::func<T>([f](T x) { return std::riemann_zeta(f(x)); });
}

}  // namespace numc