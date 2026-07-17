#pragma once

#include "../common/vector.hpp"
#include "../common/tensor.hpp"
#include "../linear_algebra/solvers.hpp"
#include "../optimalization/solvers/roots.hpp"
#include "../inc.hpp"
#include "../utility/log.hpp"

namespace numc {
namespace linear_algebra {

/// @addtogroup linear_algebra
/// @{

/// @brief Computes eigenvalues and eigenvectors of a symmetric matrix using Jacobi rotation.
template <typename T = double>
std::pair<vector<T>, tensor<T>> jacobi(tensor<T> a, T tol = T(1e-9)) {
  size_t n = a.shape()[0];
  tensor<T> p = tensor<T>::eye(n);
  int max_rot = 5 * static_cast<int>(n * n);

  auto maxElem = [&](const tensor<T>& mat, size_t& k, size_t& l) -> T {
    T aMax = T(0.0);
    for (size_t i = 0; i < n - 1; ++i) {
      for (size_t j = i + 1; j < n; ++j) {
        T val = std::abs(mat(i, j));
        if (val >= aMax) {
          aMax = val;
          k = i;
          l = j;
        }
      }
    }
    return aMax;
  };

  auto rotate = [&](tensor<T>& mat, tensor<T>& trans, size_t k, size_t l) {
    T aDiff = mat(l, l) - mat(k, k);
    T t;
    if (std::abs(mat(k, l)) < std::abs(aDiff) * T(1e-36)) {
      t = mat(k, l) / aDiff;
    } else {
      T phi = aDiff / (T(2.0) * mat(k, l));
      t = T(1.0) / (std::abs(phi) + std::sqrt(phi * phi + T(1.0)));
      if (phi < T(0.0)) t = -t;
    }
    T c = T(1.0) / std::sqrt(t * t + T(1.0));
    T s = t * c;
    T tau = s / (T(1.0) + c);
    T temp = mat(k, l);
    mat(k, l) = T(0.0);
    mat(k, k) -= t * temp;
    mat(l, l) += t * temp;

    for (size_t i = 0; i < k; ++i) {
      T tmp = mat(i, k);
      mat(i, k) = tmp - s * (mat(i, l) + tau * tmp);
      mat(i, l) = mat(i, l) + s * (tmp - tau * mat(i, l));
    }
    for (size_t i = k + 1; i < l; ++i) {
      T tmp = mat(k, i);
      mat(k, i) = tmp - s * (mat(i, l) + tau * tmp);
      mat(i, l) = mat(i, l) + s * (tmp - tau * mat(i, l));
    }
    for (size_t i = l + 1; i < n; ++i) {
      T tmp = mat(k, i);
      mat(k, i) = tmp - s * (mat(l, i) + tau * tmp);
      mat(l, i) = mat(l, i) + s * (tmp - tau * mat(l, i));
    }
    for (size_t i = 0; i < n; ++i) {
      T tmp = trans(i, k);
      trans(i, k) = tmp - s * (trans(i, l) + tau * tmp);
      trans(i, l) = trans(i, l) + s * (tmp - tau * trans(i, l));
    }
  };

  for (int r = 0; r < max_rot; ++r) {
    size_t k = 0, l = 0;
    T aMax = maxElem(a, k, l);
    if (aMax < tol) {
      vector<T> eigenvalues(n);
      for (size_t i = 0; i < n; ++i) {
        eigenvalues[i] = a(i, i);
      }
      return {eigenvalues, p};
    }
    rotate(a, p, k, l);
  }
  Log::Warn("Jacobi method did not converge.");
  vector<T> eigenvalues(n);
  for (size_t i = 0; i < n; ++i) {
    eigenvalues[i] = a(i, i);
  }
  return {eigenvalues, p};
}

/// @brief Sorts eigenvalues and eigenvectors returned by Jacobi in ascending order.
template <typename T = double>
void sort_jacobi(vector<T>& lam, tensor<T>& x) {
  size_t n = lam.size();
  for (size_t i = 0; i < n - 1; ++i) {
    size_t index = i;
    T val = lam[i];
    for (size_t j = i + 1; j < n; ++j) {
      if (lam[j] < val) {
        index = j;
        val = lam[j];
      }
    }
    if (index != i) {
      std::swap(lam[i], lam[index]);
      for (size_t k = 0; k < n; ++k) {
        std::swap(x(k, i), x(k, index));
      }
    }
  }
}

/// @brief Transforms a generalized eigenvalue problem Ax = lam Bx to standard form Hz = lam z.
template <typename T = double>
std::pair<tensor<T>, tensor<T>> std_form(tensor<T> a, tensor<T> b) {
  size_t n = a.shape()[0];
  tensor<T> L = choleski(b);
  for (size_t j = 0; j < n - 1; ++j) {
    L(j, j) = T(1.0) / L(j, j);
    for (size_t i = j + 1; i < n; ++i) {
      T sum = T(0.0);
      for (size_t k = j; k < i; ++k) {
        sum += L(i, k) * L(k, j);
      }
      L(i, j) = -sum / L(i, i);
    }
  }
  L(n - 1, n - 1) = T(1.0) / L(n - 1, n - 1);

  tensor<T> LT({n, n});
  for (size_t i = 0; i < n; ++i) {
    for (size_t j = 0; j < n; ++j) {
      LT(i, j) = L(j, i);
    }
  }
  tensor<T> H = L * a * LT;
  return {H, LT};
}

/// @brief Inverse power method to find the eigenvalue closest to s and its eigenvector.
template <typename T = double>
std::pair<T, vector<T>> inverse_power(tensor<T> a, T s, T tol = T(1e-6)) {
  size_t n = a.shape()[0];
  tensor<T> a_star = a;
  for (size_t i = 0; i < n; ++i) {
    a_star(i, i) -= s;
  }
  a_star = lu_decomp(a_star);

  vector<T> x(n);
  std::srand(12345);
  for (size_t i = 0; i < n; ++i) {
    x[i] = T(std::rand()) / T(RAND_MAX);
  }
  x = x.normalize();
  T sign = T(1.0);
  T x_mag = T(1.0);

  for (int i = 0; i < 50; ++i) {
    vector<T> x_old = x;
    x = lu_solve(a_star, x);
    x_mag = x.norm();
    x = x / x_mag;
    if (x_old * x < T(0.0)) {
      sign = T(-1.0);
      x = -x;
    } else {
      sign = T(1.0);
    }
    if ((x_old - x).norm() < tol) {
      return {s + sign / x_mag, x};
    }
  }
  Log::Warn("Inverse power method did not converge.");
  return {s + sign / x_mag, x};
}

/// @brief Inverse power method for a symmetric tridiagonal matrix.
template <typename T = double>
std::pair<T, vector<T>> inverse_power3(vector<T> d, vector<T> c, T s, T tol = T(1e-6)) {
  size_t n = d.size();
  vector<T> d_star = d;
  for (size_t i = 0; i < n; ++i) {
    d_star[i] -= s;
  }
  vector<T> e(n - 1, T(0.0));
  for (size_t i = 0; i < n - 1; ++i) {
    e[i] = c[i];
  }
  vector<T> c_star = c;
  auto [c_dec, d_dec, e_dec] = lu_decomp3(c_star, d_star, e);

  vector<T> x(n);
  std::srand(12345);
  for (size_t i = 0; i < n; ++i) {
    x[i] = T(std::rand()) / T(RAND_MAX);
  }
  x = x.normalize();
  T sign = T(1.0);
  T x_mag = T(1.0);

  for (int i = 0; i < 50; ++i) {
    vector<T> x_old = x;
    x = lu_solve3(c_dec, d_dec, e_dec, x);
    x_mag = x.norm();
    x = x / x_mag;
    if (x_old * x < T(0.0)) {
      sign = T(-1.0);
      x = -x;
    } else {
      sign = T(1.0);
    }
    if ((x_old - x).norm() < tol) {
      return {s + sign / x_mag, x};
    }
  }
  Log::Warn("Inverse power method 3 did not converge.");
  return {s + sign / x_mag, x};
}

/// @brief Inverse power method for a symmetric pentadiagonal matrix Ax = lam Bx.
template <typename T = double>
std::pair<T, vector<T>> inverse_power5(
    std::function<vector<T>(const vector<T>&)> Bv,
    vector<T> d,
    vector<T> e,
    vector<T> f,
    T tol = T(1e-6)) {
  size_t n = d.size();
  auto [d_dec, e_dec, f_dec] = lu_decomp5(d, e, f);
  vector<T> x(n);
  std::srand(12345);
  for (size_t i = 0; i < n; ++i) {
    x[i] = T(std::rand()) / T(RAND_MAX);
  }
  x = x.normalize();
  T sign = T(1.0);
  T x_mag = T(1.0);

  for (int i = 0; i < 30; ++i) {
    vector<T> x_old = x;
    vector<T> bx = Bv(x_old);
    x = lu_solve5(d_dec, e_dec, f_dec, bx);
    x_mag = x.norm();
    x = x / x_mag;
    if (x_old * x < T(0.0)) {
      sign = T(-1.0);
      x = -x;
    } else {
      sign = T(1.0);
    }
    if ((x_old - x).norm() < tol) {
      return {sign / x_mag, x};
    }
  }
  Log::Warn("Inverse power method 5 did not converge.");
  return {sign / x_mag, x};
}

/// @brief Householder reduction of a symmetric matrix to symmetric tridiagonal form.
template <typename T = double>
std::pair<vector<T>, vector<T>> householder(tensor<T>& a) {
  size_t n = a.shape()[0];
  for (size_t k = 0; k < n - 2; ++k) {
    vector<T> u(n - k - 1);
    for (size_t i = k + 1; i < n; ++i) {
      u[i - k - 1] = a(i, k);
    }
    T u_norm = u.norm();
    if (u[0] < T(0.0)) u_norm = -u_norm;
    u[0] += u_norm;

    T h = (u * u) / T(2.0);
    if (std::abs(h) < T(1e-20)) continue;

    vector<T> v(n - k - 1);
    for (size_t i = k + 1; i < n; ++i) {
      T sum = T(0.0);
      for (size_t j = k + 1; j < n; ++j) {
        sum += a(i, j) * u[j - k - 1];
      }
      v[i - k - 1] = sum / h;
    }
    T g = (u * v) / (T(2.0) * h);
    v = v - g * u;

    for (size_t i = k + 1; i < n; ++i) {
      for (size_t j = k + 1; j < n; ++j) {
        a(i, j) -= v[i - k - 1] * u[j - k - 1] + u[i - k - 1] * v[j - k - 1];
      }
    }
    for (size_t i = k + 1; i < n; ++i) {
      a(i, k) = u[i - k - 1];
    }
    a(k, k + 1) = -u_norm;
    a(k + 1, k) = -u_norm;
  }
  vector<T> d(n);
  for (size_t i = 0; i < n; ++i) d[i] = a(i, i);
  vector<T> c(n - 1);
  for (size_t i = 0; i < n - 1; ++i) c[i] = a(i, i + 1);
  return {d, c};
}

/// @brief Computes accumulated transformation matrix P after Householder reduction.
template <typename T = double>
tensor<T> householder_P(const tensor<T>& a) {
  size_t n = a.shape()[0];
  tensor<T> p = tensor<T>::eye(n);
  for (size_t k = 0; k < n - 2; ++k) {
    vector<T> u(n - k - 1);
    for (size_t i = k + 1; i < n; ++i) {
      u[i - k - 1] = a(i, k);
    }
    T h = (u * u) / T(2.0);
    if (std::abs(h) < T(1e-20)) continue;

    vector<T> v(n);
    for (size_t i = 1; i < n; ++i) {
      T sum = T(0.0);
      for (size_t j = k + 1; j < n; ++j) {
        sum += p(i, j) * u[j - k - 1];
      }
      v[i] = sum / h;
    }
    for (size_t i = 1; i < n; ++i) {
      for (size_t j = k + 1; j < n; ++j) {
        p(i, j) -= v[i] * u[j - k - 1];
      }
    }
  }
  return p;
}

/// @brief Generates Sturm sequence for a symmetric tridiagonal matrix.
template <typename T = double>
vector<T> sturm_seq(const vector<T>& d, const vector<T>& c, T lam) {
  size_t n = d.size() + 1;
  vector<T> p(n);
  for (size_t i = 0; i < n; ++i) p[i] = T(1.0);
  p[1] = d[0] - lam;
  for (size_t i = 2; i < n; ++i) {
    p[i] = (d[i - 1] - lam) * p[i - 1] - (c[i - 2] * c[i - 2]) * p[i - 2];
  }
  return p;
}

/// @brief Counts the number of sign changes in the Sturm sequence.
template <typename T = double>
int num_lambdas(const vector<T>& p) {
  size_t n = p.size();
  int sign_old = 1;
  int num_lam = 0;
  for (size_t i = 1; i < n; ++i) {
    int sign_val = 0;
    if (p[i] > T(0.0))
      sign_val = 1;
    else if (p[i] < T(0.0))
      sign_val = -1;
    else
      sign_val = -sign_old;

    if (sign_val * sign_old < 0) {
      num_lam += 1;
    }
    sign_old = sign_val;
  }
  return num_lam;
}

/// @brief Determines lower and upper bounds of eigenvalues using Gerschgorin's theorem.
template <typename T = double>
std::pair<T, T> gerschgorin(const vector<T>& d, const vector<T>& c) {
  size_t n = d.size();
  T lamMin = d[0] - std::abs(c[0]);
  T lamMax = d[0] + std::abs(c[0]);
  for (size_t i = 1; i < n - 1; ++i) {
    T lam = d[i] - std::abs(c[i]) - std::abs(c[i - 1]);
    if (lam < lamMin) lamMin = lam;
    lam = d[i] + std::abs(c[i]) + std::abs(c[i - 1]);
    if (lam > lamMax) lamMax = lam;
  }
  T lam = d[n - 1] - std::abs(c[n - 2]);
  if (lam < lamMin) lamMin = lam;
  lam = d[n - 1] + std::abs(c[n - 2]);
  if (lam > lamMax) lamMax = lam;
  return {lamMin, lamMax};
}

/// @brief Brackets the N smallest eigenvalues of a symmetric tridiagonal matrix.
template <typename T = double>
vector<T> lam_range(const vector<T>& d, const vector<T>& c, int N) {
  auto [lamMin, lamMax] = gerschgorin(d, c);
  vector<T> r(N + 1);
  for (int i = 0; i <= N; ++i) r[i] = T(1.0);
  r[0] = lamMin;

  T current_lamMax = lamMax;
  for (int k = N; k > 0; --k) {
    T lam = (current_lamMax + lamMin) / T(2.0);
    T h = (current_lamMax - lamMin) / T(2.0);
    for (int i = 0; i < 1000; ++i) {
      vector<T> p = sturm_seq(d, c, lam);
      int numLam = num_lambdas(p);
      h = h / T(2.0);
      if (numLam < k) {
        lam = lam + h;
      } else if (numLam > k) {
        lam = lam - h;
      } else {
        break;
      }
    }
    current_lamMax = lam;
    r[k] = lam;
  }
  return r;
}

/// @brief Computes N smallest eigenvalues of a symmetric tridiagonal matrix using Brent's method and Sturm sequence.
template <typename T = double>
vector<T> eigenvals3(const vector<T>& d, const vector<T>& c, int N) {
  auto f_sturm = [&](T x) -> T {
    vector<T> p = sturm_seq(d, c, x);
    return p[p.size() - 1];
  };
  vector<T> lam(N);
  vector<T> r = lam_range(d, c, N);
  for (int i = 0; i < N; ++i) {
    lam[i] = optimalization::brent<T>(f_sturm, r[i], r[i + 1]);
  }
  return lam;
}

/// @}

}  // namespace linear_algebra
}  // namespace numc
