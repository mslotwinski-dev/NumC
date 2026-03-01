#pragma once

#include "../../common/constants.hpp"
#include "../../common/function.hpp"
#include "../../inc.hpp"
#include "../../utility/log.hpp"

namespace numc {

/// @defgroup statistics Statistics
/// @brief Statistical analysis module, data processing, probability distributions, and hypothesis testing.

/// @defgroup distributions Probability Distributions
/// @ingroup statistics
/// @brief Probability Density Functions (PDF) and Cumulative Distribution Functions (CDF) for continuous random variables.
/// @{

/// @brief Computes the Probability Density Function (PDF) of the Normal (Gaussian) distribution.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param mu The mean (expectation) of the distribution, $\mu$. Default is 0.0.
/// @param sigma The standard deviation of the distribution, $\sigma$. Default is 1.0.
/// @return An expression representing the Normal PDF.
template <typename T>
numc::func<T> norm_pdf(const numc::func<T>& x, T mu = T(0.0), T sigma = T(1.0)) {
  return numc::func<T>([x, mu, sigma](T t) {
    T val = x(t);
    return std::exp(-T(0.5) * std::pow((val - mu) / sigma, T(2.0))) / (sigma * std::sqrt(T(2.0) * T(numc::PI)));
  });
}

/// @brief Computes the Cumulative Distribution Function (CDF) of the Normal distribution.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param mu The mean ($\mu$). Default is 0.0.
/// @param sigma The standard deviation ($\sigma$). Default is 1.0.
/// @return An expression representing the Normal CDF.
template <typename T>
numc::func<T> norm_cdf(const numc::func<T>& x, T mu = T(0.0), T sigma = T(1.0)) {
  return numc::func<T>([x, mu, sigma](T t) {
    T val = x(t);
    return T(0.5) * (T(1.0) + std::erf((val - mu) / (sigma * T(numc::SQRT2))));
  });
}

/// @brief Computes the PDF of the Continuous Uniform distribution.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param a The lower bound of the support interval. Default is 0.0.
/// @param b The upper bound of the support interval. Default is 1.0.
/// @return An expression representing the Uniform PDF.
template <typename T>
numc::func<T> uniform_pdf(const numc::func<T>& x, T a = T(0.0), T b = T(1.0)) {
  return numc::func<T>([x, a, b](T t) {
    T val = x(t);
    if (val < a || val > b) return T(0.0);
    return T(1.0) / (b - a);
  });
}

/// @brief Computes the CDF of the Continuous Uniform distribution.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param a The lower bound. Default is 0.0.
/// @param b The upper bound. Default is 1.0.
/// @return An expression representing the Uniform CDF.
template <typename T>
numc::func<T> uniform_cdf(const numc::func<T>& x, T a = T(0.0), T b = T(1.0)) {
  return numc::func<T>([x, a, b](T t) {
    T val = x(t);
    if (val < a) return T(0.0);
    if (val > b) return T(1.0);
    return (val - a) / (b - a);
  });
}

/// @brief Computes the PDF of the Exponential distribution.
/// Commonly used to model the time between independent events that happen at a constant average rate.
/// @tparam T Floating-point type.
/// @param x The input random variable expression ($x \ge 0$).
/// @param lambda The rate parameter ($\lambda$). Default is 1.0.
/// @return An expression representing the Exponential PDF.
template <typename T>
numc::func<T> exp_pdf(const numc::func<T>& x, T lambda = T(1.0)) {
  return numc::func<T>([x, lambda](T t) {
    T val = x(t);
    return (val >= T(0.0)) ? lambda * std::exp(-lambda * val) : T(0.0);
  });
}

/// @brief Computes the CDF of the Exponential distribution.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param lambda The rate parameter ($\lambda$). Default is 1.0.
/// @return An expression representing the Exponential CDF.
template <typename T>
numc::func<T> exp_cdf(const numc::func<T>& x, T lambda = T(1.0)) {
  return numc::func<T>([x, lambda](T t) {
    T val = x(t);
    return (val >= T(0.0)) ? T(1.0) - std::exp(-lambda * val) : T(0.0);
  });
}

/// @brief Computes the PDF of the Weibull distribution.
/// Widely used in reliability engineering and failure analysis.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param shape The shape parameter ($k > 0$).
/// @param scale The scale parameter ($\lambda > 0$). Default is 1.0.
/// @return An expression representing the Weibull PDF.
template <typename T>
numc::func<T> weibull_pdf(const numc::func<T>& x, T shape, T scale = T(1.0)) {
  return numc::func<T>([x, shape, scale](T t) {
    T val = x(t);
    if (val < T(0.0)) return T(0.0);
    return (shape / scale) * std::pow(val / scale, shape - T(1.0)) * std::exp(-std::pow(val / scale, shape));
  });
}

/// @brief Computes the CDF of the Weibull distribution.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param shape The shape parameter ($k > 0$).
/// @param scale The scale parameter ($\lambda > 0$). Default is 1.0.
/// @return An expression representing the Weibull CDF.
template <typename T>
numc::func<T> weibull_cdf(const numc::func<T>& x, T shape, T scale = T(1.0)) {
  return numc::func<T>([x, shape, scale](T t) {
    T val = x(t);
    if (val < T(0.0)) return T(0.0);
    return T(1.0) - std::exp(-std::pow(val / scale, shape));
  });
}

/// @brief Computes the PDF of the Rayleigh distribution.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param sigma The scale parameter ($\sigma > 0$). Default is 1.0.
/// @return An expression representing the Rayleigh PDF.
template <typename T>
numc::func<T> rayleigh_pdf(const numc::func<T>& x, T sigma = T(1.0)) {
  return numc::func<T>([x, sigma](T t) {
    T val = x(t);
    if (val < T(0.0)) return T(0.0);
    return (val / (sigma * sigma)) * std::exp(-(val * val) / (T(2.0) * sigma * sigma));
  });
}

/// @brief Computes the CDF of the Rayleigh distribution.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param sigma The scale parameter ($\sigma > 0$). Default is 1.0.
/// @return An expression representing the Rayleigh CDF.
template <typename T>
numc::func<T> rayleigh_cdf(const numc::func<T>& x, T sigma = T(1.0)) {
  return numc::func<T>([x, sigma](T t) {
    T val = x(t);
    if (val < T(0.0)) return T(0.0);
    return T(1.0) - std::exp(-(val * val) / (T(2.0) * sigma * sigma));
  });
}

/// @brief Computes the PDF of the Lognormal distribution.
/// @tparam T Floating-point type.
/// @param x The input random variable expression ($x > 0$).
/// @param mu The mean of the logarithm of the variable ($\mu$). Default is 0.0.
/// @param sigma The standard deviation of the logarithm of the variable ($\sigma$). Default is 1.0.
/// @return An expression representing the Lognormal PDF.
template <typename T>
numc::func<T> lognorm_pdf(const numc::func<T>& x, T mu = T(0.0), T sigma = T(1.0)) {
  return numc::func<T>([x, mu, sigma](T t) {
    T val = x(t);
    if (val <= T(0.0)) return T(0.0);
    return std::exp(-T(0.5) * std::pow((std::log(val) - mu) / sigma, T(2.0))) / (val * sigma * std::sqrt(T(2.0) * T(numc::PI)));
  });
}

/// @brief Computes the CDF of the Lognormal distribution.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param mu The parameter $\mu$. Default is 0.0.
/// @param sigma The parameter $\sigma$. Default is 1.0.
/// @return An expression representing the Lognormal CDF.
template <typename T>
numc::func<T> lognorm_cdf(const numc::func<T>& x, T mu = T(0.0), T sigma = T(1.0)) {
  return numc::func<T>([x, mu, sigma](T t) {
    T val = x(t);
    if (val <= T(0.0)) return T(0.0);
    return T(0.5) * (T(1.0) + std::erf((std::log(val) - mu) / (sigma * T(numc::SQRT2))));
  });
}

/// @brief Computes the PDF of the Laplace (double exponential) distribution.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param mu The location parameter ($\mu$). Default is 0.0.
/// @param b The scale parameter ($b > 0$). Default is 1.0.
/// @return An expression representing the Laplace PDF.
template <typename T>
numc::func<T> laplace_pdf(const numc::func<T>& x, T mu = T(0.0), T b = T(1.0)) {
  return numc::func<T>([x, mu, b](T t) {
    T val = x(t);
    return std::exp(-std::abs(val - mu) / b) / (T(2.0) * b);
  });
}

/// @brief Computes the CDF of the Laplace distribution.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param mu The location parameter. Default is 0.0.
/// @param b The scale parameter. Default is 1.0.
/// @return An expression representing the Laplace CDF.
template <typename T>
numc::func<T> laplace_cdf(const numc::func<T>& x, T mu = T(0.0), T b = T(1.0)) {
  return numc::func<T>([x, mu, b](T t) {
    T val = x(t);
    if (val < mu) return T(0.5) * std::exp((val - mu) / b);
    return T(1.0) - T(0.5) * std::exp(-(val - mu) / b);
  });
}

/// @brief Computes the PDF of the Cauchy distribution.
/// Note that its mean and variance are undefined.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param x0 The location parameter ($x_0$). Default is 0.0.
/// @param gamma The scale parameter ($\gamma > 0$). Default is 1.0.
/// @return An expression representing the Cauchy PDF.
template <typename T>
numc::func<T> cauchy_pdf(const numc::func<T>& x, T x0 = T(0.0), T gamma = T(1.0)) {
  return numc::func<T>([x, x0, gamma](T t) {
    T val = x(t);
    return T(1.0) / (T(numc::PI) * gamma * (T(1.0) + std::pow((val - x0) / gamma, T(2.0))));
  });
}

/// @brief Computes the CDF of the Cauchy distribution.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param x0 The location parameter. Default is 0.0.
/// @param gamma The scale parameter. Default is 1.0.
/// @return An expression representing the Cauchy CDF.
template <typename T>
numc::func<T> cauchy_cdf(const numc::func<T>& x, T x0 = T(0.0), T gamma = T(1.0)) {
  return numc::func<T>([x, x0, gamma](T t) {
    T val = x(t);
    return T(1.0) / T(numc::PI) * std::atan((val - x0) / gamma) + T(0.5);
  });
}

/// @brief Computes the PDF of the Pareto (power-law) distribution.
/// @tparam T Floating-point type.
/// @param x The input random variable expression ($x \ge x_m$).
/// @param xm The scale parameter (minimum possible value).
/// @param alpha The shape parameter (Pareto index).
/// @return An expression representing the Pareto PDF.
template <typename T>
numc::func<T> pareto_pdf(const numc::func<T>& x, T xm, T alpha) {
  return numc::func<T>([x, xm, alpha](T t) {
    T val = x(t);
    if (val < xm) return T(0.0);
    return alpha * std::pow(xm, alpha) / std::pow(val, alpha + T(1.0));
  });
}

/// @brief Computes the CDF of the Pareto distribution.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param xm The scale parameter.
/// @param alpha The shape parameter.
/// @return An expression representing the Pareto CDF.
template <typename T>
numc::func<T> pareto_cdf(const numc::func<T>& x, T xm, T alpha) {
  return numc::func<T>([x, xm, alpha](T t) {
    T val = x(t);
    if (val < xm) return T(0.0);
    return T(1.0) - std::pow(xm / val, alpha);
  });
}

/// @brief Computes the PDF of Student's t-distribution.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param nu The degrees of freedom parameter ($\nu > 0$).
/// @return An expression representing the Student's t PDF.
template <typename T>
numc::func<T> student_t_pdf(const numc::func<T>& x, T nu) {
  return numc::func<T>([x, nu](T t) {
    T val = x(t);
    T c = std::tgamma((nu + T(1.0)) / T(2.0)) / (std::sqrt(nu * T(numc::PI)) * std::tgamma(nu / T(2.0)));
    return c * std::pow(T(1.0) + (val * val) / nu, -(nu + T(1.0)) / T(2.0));
  });
}

/// @brief Computes the CDF of Student's t-distribution.
/// @note Currently returns 0.0 as it requires the regularized incomplete Beta function.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param nu The degrees of freedom.
/// @return An expression evaluating to 0.0 (placeholder).
template <typename T>
numc::func<T> student_t_cdf(const numc::func<T>& x, T nu) {
  return numc::func<T>([x](T t) {
    Log::Error("student_t_cdf requires the incomplete Beta function. Returning 0.");
    return T(0.0);
  });
}

/// @brief Computes the PDF of the Chi-squared ($\chi^2$) distribution.
/// Evaluated safely using logarithms to prevent floating-point overflow for large k.
/// @tparam T Floating-point type.
/// @param x The input random variable expression ($x > 0$).
/// @param k The degrees of freedom parameter ($k > 0$).
/// @return An expression representing the Chi-squared PDF.
template <typename T>
numc::func<T> chisq_pdf(const numc::func<T>& x, T k) {
  return numc::func<T>([x, k](T t) {
    T val = x(t);
    if (val <= T(0.0)) return T(0.0);
    return std::exp((k / T(2.0) - T(1.0)) * std::log(val) - val / T(2.0) - (k / T(2.0)) * std::log(T(2.0)) - std::lgamma(k / T(2.0)));
  });
}

/// @brief Computes the CDF of the Chi-squared distribution.
/// @note Currently returns 0.0 as it requires the lower incomplete Gamma function.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param k The degrees of freedom.
/// @return An expression evaluating to 0.0 (placeholder).
template <typename T>
numc::func<T> chisq_cdf(const numc::func<T>& x, T k) {
  return numc::func<T>([x](T t) {
    Log::Error("chisq_cdf requires the incomplete Gamma function. Returning 0.");
    return T(0.0);
  });
}

/// @brief Computes the PDF of the F-distribution (Snedecor's F).
/// @tparam T Floating-point type.
/// @param x The input random variable expression ($x > 0$).
/// @param d1 The numerator degrees of freedom ($d_1 > 0$).
/// @param d2 The denominator degrees of freedom ($d_2 > 0$).
/// @return An expression representing the F-distribution PDF.
template <typename T>
numc::func<T> f_pdf(const numc::func<T>& x, T d1, T d2) {
  return numc::func<T>([x, d1, d2](T t) {
    T val = x(t);
    if (val <= T(0.0)) return T(0.0);
    T log_num = (d1 / T(2.0)) * std::log(d1 * val) + (d2 / T(2.0)) * std::log(d2) - ((d1 + d2) / T(2.0)) * std::log(d1 * val + d2);
    return std::exp(log_num) / (val * std::beta(d1 / T(2.0), d2 / T(2.0)));
  });
}

/// @brief Computes the CDF of the F-distribution.
/// @note Currently returns 0.0 as it requires the regularized incomplete Beta function.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param d1 The numerator degrees of freedom.
/// @param d2 The denominator degrees of freedom.
/// @return An expression evaluating to 0.0 (placeholder).
template <typename T>
numc::func<T> f_cdf(const numc::func<T>& x, T d1, T d2) {
  return numc::func<T>([x](T t) {
    Log::Error("f_cdf requires the incomplete Beta function. Returning 0.");
    return T(0.0);
  });
}

/// @brief Computes the PDF of the Gamma distribution.
/// @tparam T Floating-point type.
/// @param x The input random variable expression ($x > 0$).
/// @param shape The shape parameter ($k$ or $\alpha > 0$).
/// @param scale The scale parameter ($\theta$ or $\beta > 0$). Default is 1.0.
/// @return An expression representing the Gamma PDF.
template <typename T>
numc::func<T> gamma_pdf(const numc::func<T>& x, T shape, T scale = T(1.0)) {
  return numc::func<T>([x, shape, scale](T t) {
    T val = x(t);
    if (val <= T(0.0)) return T(0.0);
    return std::exp((shape - T(1.0)) * std::log(val) - val / scale - shape * std::log(scale) - std::lgamma(shape));
  });
}

/// @brief Computes the CDF of the Gamma distribution.
/// @note Currently returns 0.0 as it requires the lower incomplete Gamma function.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param shape The shape parameter.
/// @param scale The scale parameter. Default is 1.0.
/// @return An expression evaluating to 0.0 (placeholder).
template <typename T>
numc::func<T> gamma_cdf(const numc::func<T>& x, T shape, T scale = T(1.0)) {
  return numc::func<T>([x](T t) {
    Log::Error("gamma_cdf requires the incomplete Gamma function. Returning 0.");
    return T(0.0);
  });
}

/// @brief Computes the PDF of the Beta distribution.
/// Defined on the interval [0, 1], highly useful for modeling proportions and probabilities.
/// @tparam T Floating-point type.
/// @param x The input random variable expression ($x \in [0, 1]$).
/// @param alpha The first shape parameter ($\alpha > 0$).
/// @param beta The second shape parameter ($\beta > 0$).
/// @return An expression representing the Beta PDF.
template <typename T>
numc::func<T> beta_pdf(const numc::func<T>& x, T alpha, T beta) {
  return numc::func<T>([x, alpha, beta](T t) {
    T val = x(t);
    if (val < T(0.0) || val > T(1.0)) return T(0.0);
    return std::pow(val, alpha - T(1.0)) * std::pow(T(1.0) - val, beta - T(1.0)) / std::beta(alpha, beta);
  });
}

/// @brief Computes the CDF of the Beta distribution.
/// @note Currently returns 0.0 as it requires the regularized incomplete Beta function.
/// @tparam T Floating-point type.
/// @param x The input random variable expression.
/// @param alpha The first shape parameter.
/// @param beta The second shape parameter.
/// @return An expression evaluating to 0.0 (placeholder).
template <typename T>
numc::func<T> beta_cdf(const numc::func<T>& x, T alpha, T beta) {
  return numc::func<T>([x](T t) {
    Log::Error("beta_cdf requires the incomplete Beta function. Returning 0.");
    return T(0.0);
  });
}

/// @} // End of distributions group

}  // namespace numc