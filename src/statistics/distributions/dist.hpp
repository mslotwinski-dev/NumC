#pragma once

#include "../../common/constants.hpp"
#include "../../common/function.hpp"
#include "../../inc.hpp"
#include "../../utility/log.hpp"

namespace numc {

template <typename T>
numc::func<T> norm_pdf(const numc::func<T>& x, T mu = T(0.0), T sigma = T(1.0)) {
  return numc::func<T>([x, mu, sigma](T t) {
    T val = x(t);
    return std::exp(-T(0.5) * std::pow((val - mu) / sigma, T(2.0))) / (sigma * std::sqrt(T(2.0) * T(PI)));
  });
}

template <typename T>
numc::func<T> norm_cdf(const numc::func<T>& x, T mu = T(0.0), T sigma = T(1.0)) {
  return numc::func<T>([x, mu, sigma](T t) {
    T val = x(t);
    return T(0.5) * (T(1.0) + std::erf((val - mu) / (sigma * T(SQRT2))));
  });
}

template <typename T>
numc::func<T> uniform_pdf(const numc::func<T>& x, T a = T(0.0), T b = T(1.0)) {
  return numc::func<T>([x, a, b](T t) {
    T val = x(t);
    if (val < a || val > b) return T(0.0);
    return T(1.0) / (b - a);
  });
}

template <typename T>
numc::func<T> uniform_cdf(const numc::func<T>& x, T a = T(0.0), T b = T(1.0)) {
  return numc::func<T>([x, a, b](T t) {
    T val = x(t);
    if (val < a) return T(0.0);
    if (val > b) return T(1.0);
    return (val - a) / (b - a);
  });
}

template <typename T>
numc::func<T> exp_pdf(const numc::func<T>& x, T lambda = T(1.0)) {
  return numc::func<T>([x, lambda](T t) {
    T val = x(t);
    return (val >= T(0.0)) ? lambda * std::exp(-lambda * val) : T(0.0);
  });
}

template <typename T>
numc::func<T> exp_cdf(const numc::func<T>& x, T lambda = T(1.0)) {
  return numc::func<T>([x, lambda](T t) {
    T val = x(t);
    return (val >= T(0.0)) ? T(1.0) - std::exp(-lambda * val) : T(0.0);
  });
}

template <typename T>
numc::func<T> weibull_pdf(const numc::func<T>& x, T shape, T scale = T(1.0)) {
  return numc::func<T>([x, shape, scale](T t) {
    T val = x(t);
    if (val < T(0.0)) return T(0.0);
    return (shape / scale) * std::pow(val / scale, shape - T(1.0)) * std::exp(-std::pow(val / scale, shape));
  });
}

template <typename T>
numc::func<T> weibull_cdf(const numc::func<T>& x, T shape, T scale = T(1.0)) {
  return numc::func<T>([x, shape, scale](T t) {
    T val = x(t);
    if (val < T(0.0)) return T(0.0);
    return T(1.0) - std::exp(-std::pow(val / scale, shape));
  });
}

template <typename T>
numc::func<T> rayleigh_pdf(const numc::func<T>& x, T sigma = T(1.0)) {
  return numc::func<T>([x, sigma](T t) {
    T val = x(t);
    if (val < T(0.0)) return T(0.0);
    return (val / (sigma * sigma)) * std::exp(-(val * val) / (T(2.0) * sigma * sigma));
  });
}

template <typename T>
numc::func<T> rayleigh_cdf(const numc::func<T>& x, T sigma = T(1.0)) {
  return numc::func<T>([x, sigma](T t) {
    T val = x(t);
    if (val < T(0.0)) return T(0.0);
    return T(1.0) - std::exp(-(val * val) / (T(2.0) * sigma * sigma));
  });
}

template <typename T>
numc::func<T> lognorm_pdf(const numc::func<T>& x, T mu = T(0.0), T sigma = T(1.0)) {
  return numc::func<T>([x, mu, sigma](T t) {
    T val = x(t);
    if (val <= T(0.0)) return T(0.0);
    return std::exp(-T(0.5) * std::pow((std::log(val) - mu) / sigma, T(2.0))) / (val * sigma * std::sqrt(T(2.0) * T(PI)));
  });
}

template <typename T>
numc::func<T> lognorm_cdf(const numc::func<T>& x, T mu = T(0.0), T sigma = T(1.0)) {
  return numc::func<T>([x, mu, sigma](T t) {
    T val = x(t);
    if (val <= T(0.0)) return T(0.0);
    return T(0.5) * (T(1.0) + std::erf((std::log(val) - mu) / (sigma * T(SQRT2))));
  });
}

template <typename T>
numc::func<T> laplace_pdf(const numc::func<T>& x, T mu = T(0.0), T b = T(1.0)) {
  return numc::func<T>([x, mu, b](T t) {
    T val = x(t);
    return std::exp(-std::abs(val - mu) / b) / (T(2.0) * b);
  });
}

template <typename T>
numc::func<T> laplace_cdf(const numc::func<T>& x, T mu = T(0.0), T b = T(1.0)) {
  return numc::func<T>([x, mu, b](T t) {
    T val = x(t);
    if (val < mu) return T(0.5) * std::exp((val - mu) / b);
    return T(1.0) - T(0.5) * std::exp(-(val - mu) / b);
  });
}

template <typename T>
numc::func<T> cauchy_pdf(const numc::func<T>& x, T x0 = T(0.0), T gamma = T(1.0)) {
  return numc::func<T>([x, x0, gamma](T t) {
    T val = x(t);
    return T(1.0) / (T(PI) * gamma * (T(1.0) + std::pow((val - x0) / gamma, T(2.0))));
  });
}

template <typename T>
numc::func<T> cauchy_cdf(const numc::func<T>& x, T x0 = T(0.0), T gamma = T(1.0)) {
  return numc::func<T>([x, x0, gamma](T t) {
    T val = x(t);
    return T(1.0) / T(PI) * std::atan((val - x0) / gamma) + T(0.5);
  });
}

template <typename T>
numc::func<T> pareto_pdf(const numc::func<T>& x, T xm, T alpha) {
  return numc::func<T>([x, xm, alpha](T t) {
    T val = x(t);
    if (val < xm) return T(0.0);
    return alpha * std::pow(xm, alpha) / std::pow(val, alpha + T(1.0));
  });
}

template <typename T>
numc::func<T> pareto_cdf(const numc::func<T>& x, T xm, T alpha) {
  return numc::func<T>([x, xm, alpha](T t) {
    T val = x(t);
    if (val < xm) return T(0.0);
    return T(1.0) - std::pow(xm / val, alpha);
  });
}

template <typename T>
numc::func<T> student_t_pdf(const numc::func<T>& x, T nu) {
  return numc::func<T>([x, nu](T t) {
    T val = x(t);
    T c = std::tgamma((nu + T(1.0)) / T(2.0)) / (std::sqrt(nu * T(PI)) * std::tgamma(nu / T(2.0)));
    return c * std::pow(T(1.0) + (val * val) / nu, -(nu + T(1.0)) / T(2.0));
  });
}

template <typename T>
numc::func<T> student_t_cdf(const numc::func<T>& x, T nu) {
  return numc::func<T>([x](T t) {
    Log::Error("student_t_cdf wymaga niekompletnej funkcji Beta. Zwracam 0.");
    return T(0.0);
  });
}

template <typename T>
numc::func<T> chisq_pdf(const numc::func<T>& x, T k) {
  return numc::func<T>([x, k](T t) {
    T val = x(t);
    if (val <= T(0.0)) return T(0.0);
    return std::exp((k / T(2.0) - T(1.0)) * std::log(val) - val / T(2.0) - (k / T(2.0)) * std::log(T(2.0)) - std::lgamma(k / T(2.0)));
  });
}

template <typename T>
numc::func<T> chisq_cdf(const numc::func<T>& x, T k) {
  return numc::func<T>([x](T t) {
    Log::Error("chisq_cdf wymaga niekompletnej funkcji Gamma. Zwracam 0.");
    return T(0.0);
  });
}

template <typename T>
numc::func<T> f_pdf(const numc::func<T>& x, T d1, T d2) {
  return numc::func<T>([x, d1, d2](T t) {
    T val = x(t);
    if (val <= T(0.0)) return T(0.0);
    T log_num = (d1 / T(2.0)) * std::log(d1 * val) + (d2 / T(2.0)) * std::log(d2) - ((d1 + d2) / T(2.0)) * std::log(d1 * val + d2);
    return std::exp(log_num) / (val * std::beta(d1 / T(2.0), d2 / T(2.0)));
  });
}

template <typename T>
numc::func<T> f_cdf(const numc::func<T>& x, T d1, T d2) {
  return numc::func<T>([x](T t) {
    Log::Error("f_cdf wymaga niekompletnej funkcji Beta. Zwracam 0.");
    return T(0.0);
  });
}

template <typename T>
numc::func<T> gamma_pdf(const numc::func<T>& x, T shape, T scale = T(1.0)) {
  return numc::func<T>([x, shape, scale](T t) {
    T val = x(t);
    if (val <= T(0.0)) return T(0.0);
    return std::exp((shape - T(1.0)) * std::log(val) - val / scale - shape * std::log(scale) - std::lgamma(shape));
  });
}

template <typename T>
numc::func<T> gamma_cdf(const numc::func<T>& x, T shape, T scale = T(1.0)) {
  return numc::func<T>([x](T t) {
    Log::Error("gamma_cdf wymaga niekompletnej funkcji Gamma. Zwracam 0.");
    return T(0.0);
  });
}

template <typename T>
numc::func<T> beta_pdf(const numc::func<T>& x, T alpha, T beta) {
  return numc::func<T>([x, alpha, beta](T t) {
    T val = x(t);
    if (val < T(0.0) || val > T(1.0)) return T(0.0);
    return std::pow(val, alpha - T(1.0)) * std::pow(T(1.0) - val, beta - T(1.0)) / std::beta(alpha, beta);
  });
}

template <typename T>
numc::func<T> beta_cdf(const numc::func<T>& x, T alpha, T beta) {
  return numc::func<T>([x](T t) {
    Log::Error("beta_cdf wymaga niekompletnej funkcji Beta. Zwracam 0.");
    return T(0.0);
  });
}

}  // namespace numc