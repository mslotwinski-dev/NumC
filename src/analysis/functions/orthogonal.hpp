#pragma once

#include "../../common/function.hpp"
#include "../../inc.hpp"
#include "../../utility/log.hpp"

namespace numc {

template <typename T>
numc::func<T> bessel_j(T nu, const numc::func<T>& f) {
  return numc::func<T>([nu, f](T x) { return std::cyl_bessel_j(nu, f(x)); });
}

template <typename T>
numc::func<T> bessel_y(T nu, const numc::func<T>& f) {
  // W standardzie C++ funkcja Y (Neumanna) nazywa się cyl_neumann
  return numc::func<T>([nu, f](T x) { return std::cyl_neumann(nu, f(x)); });
}

template <typename T>
numc::func<T> bessel_i(T nu, const numc::func<T>& f) {
  return numc::func<T>([nu, f](T x) { return std::cyl_bessel_i(nu, f(x)); });
}

template <typename T>
numc::func<T> bessel_k(T nu, const numc::func<T>& f) {
  return numc::func<T>([nu, f](T x) { return std::cyl_bessel_k(nu, f(x)); });
}

template <typename T>
numc::func<T> sph_bessel_j(unsigned n, const numc::func<T>& f) {
  return numc::func<T>([n, f](T x) { return std::sph_bessel(n, f(x)); });
}

template <typename T>
numc::func<T> sph_bessel_y(unsigned n, const numc::func<T>& f) {
  return numc::func<T>([n, f](T x) { return std::sph_neumann(n, f(x)); });
}

template <typename T>
numc::func<T> sph_legendre(unsigned l, unsigned m, const numc::func<T>& theta) {
  return numc::func<T>([l, m, theta](T x) { return std::sph_legendre(l, m, theta(x)); });
}

template <typename T>
numc::func<T> legendre(unsigned n, const numc::func<T>& f) {
  return numc::func<T>([n, f](T x) { return std::legendre(n, f(x)); });
}

template <typename T>
numc::func<T> assoc_legendre(unsigned l, unsigned m, const numc::func<T>& f) {
  return numc::func<T>([l, m, f](T x) { return std::assoc_legendre(l, m, f(x)); });
}

template <typename T>
numc::func<T> hermite(unsigned n, const numc::func<T>& f) {
  return numc::func<T>([n, f](T x) { return std::hermite(n, f(x)); });
}

template <typename T>
numc::func<T> laguerre(unsigned n, const numc::func<T>& f) {
  return numc::func<T>([n, f](T x) { return std::laguerre(n, f(x)); });
}

template <typename T>
numc::func<T> assoc_laguerre(unsigned n, unsigned m, const numc::func<T>& f) {
  return numc::func<T>([n, m, f](T x) { return std::assoc_laguerre(n, m, f(x)); });
}

template <typename T>
numc::func<T> chebyshev(unsigned n, const numc::func<T>& f) {
  return numc::func<T>([n, f](T x) {
    T val = f(x);
    if (n == 0) return T(1.0);
    if (n == 1) return val;

    T t0 = T(1.0), t1 = val, tn = T(0.0);
    for (unsigned i = 2; i <= n; ++i) {
      tn = T(2.0) * val * t1 - t0;
      t0 = t1;
      t1 = tn;
    }
    return tn;
  });
}

template <typename T>
numc::func<T> chebyshev_u(unsigned n, const numc::func<T>& f) {
  return numc::func<T>([n, f](T x) {
    T val = f(x);
    if (n == 0) return T(1.0);
    if (n == 1) return T(2.0) * val;

    T u0 = T(1.0), u1 = T(2.0) * val, un = T(0.0);
    for (unsigned i = 2; i <= n; ++i) {
      un = T(2.0) * val * u1 - u0;
      u0 = u1;
      u1 = un;
    }
    return un;
  });
}

template <typename T>
numc::func<T> gegenbauer(unsigned n, T alpha, const numc::func<T>& f) {
  return numc::func<T>([n, alpha, f](T x) {
    T val = f(x);
    if (n == 0) return T(1.0);
    if (n == 1) return T(2.0) * alpha * val;

    T c0 = T(1.0), c1 = T(2.0) * alpha * val, cn = T(0.0);
    for (unsigned i = 2; i <= n; ++i) {
      cn = (T(2.0) * val * (T(i) + alpha - T(1.0)) * c1 - (T(i) + T(2.0) * alpha - T(2.0)) * c0) / T(i);
      c0 = c1;
      c1 = cn;
    }
    return cn;
  });
}

template <typename T>
numc::func<T> jacobi(unsigned n, T alpha, T beta, const numc::func<T>& f) {
  return numc::func<T>([n, alpha, beta, f](T x) {
    T val = f(x);
    if (n == 0) return T(1.0);

    T p0 = T(1.0);
    T p1 = T(0.5) * (alpha - beta + (alpha + beta + T(2.0)) * val);
    if (n == 1) return p1;

    T pn = T(0.0);
    for (unsigned i = 2; i <= n; ++i) {
      T Ti = T(i);
      T c1 = T(2.0) * Ti * (Ti + alpha + beta) * (T(2.0) * Ti + alpha + beta - T(2.0));
      T c2 = (T(2.0) * Ti + alpha + beta - T(1.0)) *
             ((T(2.0) * Ti + alpha + beta) * (T(2.0) * Ti + alpha + beta - T(2.0)) * val + alpha * alpha - beta * beta);
      T c3 = T(2.0) * (Ti + alpha - T(1.0)) * (Ti + beta - T(1.0)) * (T(2.0) * Ti + alpha + beta);

      pn = (c2 * p1 - c3 * p0) / c1;
      p0 = p1;
      p1 = pn;
    }
    return pn;
  });
}

}  // namespace numc