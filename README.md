<div align="center">

<img src="public/logo.png">

### The Ultimate Mathematics & Scientific Computing Framework for Modern C++

*Write math the way you think it — not the way your compiler demands it.*

[![C++20](https://img.shields.io/badge/C%2B%2B-20-blue?style=flat-square&logo=cplusplus)](https://en.cppreference.com/w/cpp/20)
[![Header-Only](https://img.shields.io/badge/Header--Only-99%25-brightgreen?style=flat-square)](.)
[![License](https://img.shields.io/badge/license-MIT-orange?style=flat-square)](LICENSE)
[![Platform](https://img.shields.io/badge/platform-Windows%20%7C%20Linux%20%7C%20macOS-lightgrey?style=flat-square)](.)

</div>

---

## 🌟 What is NumC?

**NumC** is a modern, high-performance, purely declarative mathematical framework written entirely in **C++20** with **zero external dependencies**. It brings the expressive elegance of Python's NumPy/SciPy/Matplotlib ecosystem into the world of compiled, statically-typed C++ — without sacrificing a single cycle of performance.

Everything — from neural networks to FFT to publication-quality SVG rendering — is implemented **from scratch in pure standard C++**.

```cpp
#include <numc>
using namespace numc;

int main() {
    auto x = funcv::x;

    // Declare math exactly like you'd write it on paper
    func<double> f = sin(x) * exp(-0.1 * x) + (x^2) / 10.0;

    // Evaluate, differentiate, integrate — one-liners
    double y     = f(3.14);                         // Evaluate
    double slope = derivative(f, 2.0);              // f'(2.0)
    double area  = integral(f, 0.0, PI);            // ∫₀^π f(x) dx

    // Visualize in a native window with one call
    plt::show(f, -10, 10, "My Function");
}
```

---

## ✨ Key Features

| Feature | Description |
|---------|-------------|
| 🧠 **Lazy Evaluation (AST)** | Mathematical expressions build expression trees internally. Functions are evaluated only when needed, enabling composition, differentiation, and visualization without manual loops. |
| 🎨 **Publication-Quality Plots** | Built-in SVG renderer (Light/Dark themes for papers/screens) **and** native Win32 GDI+ windows — no Python, no Matplotlib, no dependencies. |
| ⚡ **Zero Dependencies** | Everything is implemented in pure C++20. No Eigen, no Boost, no FFTW. Just `#include <numc>` and go. |
| 📐 **Comprehensive Math** | 14+ modules covering Linear Algebra, Calculus, Statistics, ML, Differential Equations, Graphs, Geometry, Number Theory, and more. |
| 🔬 **Research-Grade** | Levenberg-Marquardt fitting, Jacobi eigenvalue decomposition, adaptive Runge-Kutta ODE solvers, Miller-Rabin primality, PCA — all battle-tested. |

---

## 📦 Architecture Overview

```
numc/
├── common/          → vector, tensor, func, complex, polynomial, sparse_matrix, DataFrame, constants
├── analysis/        → differentiation, integration, interpolation, FFT, calculus, special functions
│   ├── functions/   → elementary (sin, cos, exp, ...), special (Γ, β, erf, ζ), orthogonal (Bessel, Legendre, Hermite), signal (Heaviside, sinc, rect)
│   └── multivariable/ → nabla (∇), gradient, divergence, curl, Jacobian, Hessian
├── linear_algebra/  → Gauss, LU, Cholesky, Jacobi eigenvalues, matrix operations
├── optimalization/  → root finding (bisection, Newton, Halley, secant), minimization (golden search), least squares (Levenberg-Marquardt)
├── differential_eq/ → ODE (RK4, RK5, Taylor), BVP (shooting), PDE (heat, wave, Laplace)
├── statistics/      → mean, variance, regression, histograms, distributions, fitting, error propagation, param
├── machine_learning/→ Neural Networks (MLP + Adam), K-Means, PCA, Decision Trees, KNN
├── graphs/          → Dijkstra, Bellman-Ford, A*, BFS, Floyd-Warshall, MST (Kruskal, Prim)
├── geometry/        → Convex Hull, Bézier curves, B-Splines, point-in-polygon, distances
├── number_theory/   → Sieve of Eratosthenes, Miller-Rabin, modular arithmetic, GCD/LCM, CRT
└── visualization/   → SVG renderer, native Win32 GDI+ windows, Dark/Light themes
```

---

## 🚀 Quick Start

### The Magic of `func<T>` — Declarative Mathematics

The heart of NumC is the `func<T>` class. It lets you compose arbitrarily complex mathematical expressions using natural syntax. Under the hood, it builds a lazy-evaluated expression tree.

```cpp
#include <numc>
using namespace numc;

int main() {
    // Global symbolic variables — just like in a math textbook
    auto x = funcv::x;
    auto y = funcv::y;
    auto z = funcv::z;

    // Compose functions naturally (lazy — nothing is computed yet)
    func<double> gaussian = (1.0 / sqrt(2.0 * PI)) * exp(-(x^2) / 2.0);

    // Evaluate at a point
    double peak = gaussian(0.0);   // ≈ 0.3989

    // Multivariable functions work the same way
    func<double> field = sin(x) * cos(y) + z;
    double val = field(1.0, 2.0, 3.0);   // sin(1)*cos(2) + 3

    // Map over a vector (element-wise evaluation)
    vector<double> xs = {0.0, 1.0, 2.0, 3.0};
    vector<double> ys = gaussian.map(xs);   // Evaluate on every point
}
```

### Operator Overloading — Math as Code

NumC overloads all arithmetic operators for `func<T>`, so you write code that reads like a formula:

```cpp
func<double> a = sin(x) + cos(x);          // f(x) = sin(x) + cos(x)
func<double> b = a * 2.0 - x;              // g(x) = 2·(sin(x) + cos(x)) - x
func<double> c = (x^3) - 3.0 * (x^2) + 1; // h(x) = x³ - 3x² + 1
func<double> d = exp(-x) * sin(10.0 * x);  // Damped oscillation
```

---

## 🧩 Modules Deep Dive

### 1. 📈 Visualization (`numc::plot`, `numc::plt`)

NumC ships with a complete, built-from-scratch visualization engine. No external tools required.

**Three rendering backends:**
- `PlotBackend::WINDOW` — Native Win32 GDI+ window (default). Resizable, anti-aliased, with ESC to close.
- `PlotBackend::SVG` — Generates publication-quality SVG files. Perfect for LaTeX, papers, and reports.
- `PlotBackend::WEB` — Interactive HTML charts using Plotly.js CDN.

**Two themes:**
- `PlotTheme::DARK` — Vibrant GitHub-Dark inspired palette. Default for on-screen analysis.
- `PlotTheme::LIGHT` — Clean white background with precise grays. Built for scientific publication.

```cpp
auto x = funcv::x;
func<double> f1 = sin(x) * (x / 2.0);
func<double> f2 = exp(-0.1 * x) * cos(x) * 5.0;

vector<double> data_x = {-8.0, -5.5, -2.1, 0.0, 3.14, 6.28, 9.0};
vector<double> data_y = {3.5, -2.1, -1.0, 5.2, -4.8, 5.1, -4.0};

// Full dashboard with functions + scatter data
plot dashboard("Vibration Analysis");
dashboard.set_xlabel("Time [ms]").set_ylabel("Amplitude [V]");
dashboard.add(f1, -10.0, 10.0, "Theoretical model");
dashboard.add(f2, -10.0, 10.0, "Damped model");
dashboard.add_scatter(data_x, data_y, "Sensor measurements");

dashboard.show();                                      // Native GDI+ window
dashboard.show(PlotBackend::SVG, PlotTheme::LIGHT);    // Publication SVG
dashboard.show(PlotTheme::LIGHT);                      // Window with light theme

// One-liner quick plot
plt::show<double>({f1, f2}, -10.0, 10.0, "Quick Plot", "x", "f(x)");
```

---

### 2. 🧮 Linear Algebra (`numc::linear_algebra`)

Full suite of solvers, decompositions, and matrix operations.

**Solvers:** Gauss elimination, LU decomposition & solve, Cholesky (LLᵀ), Gauss-Seidel iteration, conjugate gradient.

**Eigenvalues:** Jacobi rotation for symmetric matrices, inverse power iteration.

**Matrix ops:** Transpose, determinant, inverse, matrix multiplication, identity.

```cpp
// Solve a 3×3 system: Ax = b
tensor<double> A = {{4, -1, 1}, {-1, 4, -1}, {1, -1, 4}};
vector<double> b = {12, -1, 5};

vector<double> x = linear_algebra::gauss_elimin(A, b);  // x ≈ {3, 1, 1}

// LU decomposition
tensor<double> LU = linear_algebra::lu_decomp(A);
vector<double> x2 = linear_algebra::lu_solve(LU, b);

// Cholesky for positive-definite systems
tensor<double> L = linear_algebra::cholesky(A);
vector<double> x3 = linear_algebra::cholesky_solve(L, b);

// Eigenvalues of a symmetric matrix
auto [eigenvalues, eigenvectors] = linear_algebra::jacobi(A);
```

---

### 3. 📐 Calculus & Analysis (`numc::analysis`)

A full numerical analysis toolkit — from derivatives to Taylor series.

```cpp
auto x = funcv::x;
func<double> f = sin(x) * (x^2);

// ═══════ DERIVATIVES ═══════
double d1 = analysis::derivative_central(f, 2.0);     // f'(2.0), O(h²)
double d2 = analysis::derivative_central_O4(f, 2.0);  // f'(2.0), O(h⁴)
double d3 = analysis::derivative_second_central(f, 2.0); // f''(2.0)

// Generate f'(x) as a new func — composable!
func<double> df = derivative(f);
double slope_at_pi = df(PI);

// ═══════ INTEGRALS ═══════
double area1 = analysis::simpson_1_3(f, 0.0, PI);     // Simpson's 1/3 rule
double area2 = analysis::romberg(f, 0.0, PI).first;   // Romberg (high precision)
double area3 = analysis::gauss_legendre<double, 5>(f, 0.0, PI); // 5-point Gauss-Legendre

// Built-in on func itself (Adaptive Simpson)
double area4 = f.integral(0.0, PI);

// 2D integrals
func<double> g = sin(x) * cos(funcv::y);
double volume = g.integral_2d(0.0, PI, 0.0, PI);

// Monte Carlo N-dimensional integration
auto [mc_val, mc_err] = f.integral_mc({0.0}, {PI}, 100000);

// ═══════ LIMITS ═══════
func<double> sinc_f = sin(x) / x;
double lim_val = limit(sinc_f, 0.0);   // → 1.0

// ═══════ SERIES ═══════
func<double> term = 1.0 / (x^2);   // 1/n²
double basel = analysis::series_sum_convergent(term, 1);  // → π²/6 ≈ 1.6449

// ═══════ TAYLOR EXPANSION ═══════
vector<double> coeffs = taylor_coeffs(f, 0.0, 5);  // 5th-order Taylor around x=0
double approx = taylor_eval(coeffs, 0.0, 0.5);     // Evaluate approximation

// ═══════ GEOMETRY OF CURVES ═══════
double length = arc_length(f, 0.0, PI);       // Arc length
double kappa = analysis::curvature(f, 1.0);   // Curvature at x=1
```

---

### 4. 🔢 Special Functions

NumC provides a complete library of special mathematical functions, all composable with `func<T>`.

```cpp
auto x = funcv::x;

// Error functions
func<double> err = erf(x);         // Gauss error function
func<double> errc = erfc(x);      // Complementary error function

// Gamma & Beta
func<double> gamma_f = tgamma(x); // Γ(x) — generalized factorial
func<double> beta_f = beta(x, x); // B(x, x)

// Bessel functions (cylindrical & spherical)
func<double> J0 = bessel_j(0.0, x);   // J₀(x)
func<double> J1 = bessel_j(1.0, x);   // J₁(x)
func<double> Y0 = bessel_y(0.0, x);   // Y₀(x) — Neumann function

// Orthogonal polynomials
func<double> P3 = legendre(3, x);     // Legendre P₃(x)
func<double> H4 = hermite(4, x);      // Hermite H₄(x)
func<double> L2 = laguerre(2, x);     // Laguerre L₂(x)
func<double> T5 = chebyshev_t(5, x);  // Chebyshev T₅(x)

// Elliptic integrals
func<double> K = comp_ellint_1(x);    // Complete elliptic integral K(k)
func<double> E = comp_ellint_2(x);    // Complete elliptic integral E(k)

// Number theory meets analysis
func<double> zeta = riemann_zeta(x);  // ζ(x) — Riemann zeta function
func<double> Ei = expint(x);          // Ei(x) — Exponential integral
```

---

### 5. 📡 Signal Processing Functions

```cpp
auto x = funcv::x;

func<double> step = heaviside(x);              // Heaviside step H(x)
func<double> impulse = dirac(x);               // Dirac delta δ(x)
func<double> r = ramp(x);                      // Ramp R(x) = max(0, x)
func<double> sc = sinc(x);                     // sinc(x) = sin(x)/x
func<double> rect_pulse = rect(x, 0.0, 1.0);  // Rectangular pulse [0, 1]
func<double> tri = triangle(x, 1.0);           // Triangle wave (period 1)
func<double> sq = square_wave(x, 1.0);         // Square wave (period 1)
func<double> sw = sawtooth(x, 1.0);            // Sawtooth wave (period 1)
```

---

### 6. 🔍 Root Finding & Optimization (`numc::optimalization`)

Seven root-finding algorithms, plus minimization and least-squares fitting.

```cpp
auto x = funcv::x;
func<double> f = (x^3) - 2.0 * x - 5.0;  // Find where f(x) = 0

// ═══════ ROOT FINDING ═══════
double r1 = bisection(f, 2.0, 3.0);           // Bisection method
double r2 = regula_falsi(f, 2.0, 3.0);        // False position
double r3 = newton(f, 2.5);                    // Newton-Raphson
double r4 = halley(f, 2.5);                    // Halley's method (cubic convergence!)
double r5 = secant(f, 2.0, 3.0);              // Secant method
double r6 = inverse_quadratic(f, 2.0, 2.5, 3.0); // Inverse quadratic interpolation
double r7 = optimalization::fixed_point(cos(x), 0.5); // Fixed-point iteration

// ═══════ MINIMIZATION ═══════
func<double> g = (x - 3.0)^2 + 1.0;

auto bracket = optimalization::gold_bracket(g, 0.0, 0.5);  // Bracket the minimum
auto [x_min, f_min] = optimalization::gold_search(g, bracket.first, bracket.second);
// x_min ≈ 3.0, f_min ≈ 1.0
```

---

### 7. 🌊 Differential Equations (`numc::differential_equations`)

Solve ODEs, systems of ODEs, boundary value problems, and PDEs.

```cpp
// ═══════ ODE: RUNGE-KUTTA 4th ORDER ═══════
// Solve y' = -0.1*y, y(0) = 100 (exponential decay)
auto F = [](double t, const vector<double>& y) -> vector<double> {
    return {-0.1 * y[0]};
};
auto [t_vals, y_vals] = differential_equations::run_kut4(F, 0.0, {100.0}, 50.0, 0.5);

// ═══════ ADAPTIVE RK5 (Cash-Karp) ═══════
auto [t2, y2] = differential_equations::run_kut5(F, 0.0, {100.0}, 50.0, 0.5);

// ═══════ PDE: 1D HEAT EQUATION ═══════
// u_t = α · u_xx with boundary conditions
vector<double> u0(50);   // Initial temperature distribution
for (size_t i = 0; i < 50; ++i) u0[i] = std::sin(PI * i / 49.0);

auto solution = differential_equations::heat_explicit_1d(u0, 0.01, 0.02, 0.001, 1000);

// ═══════ PDE: 1D WAVE EQUATION ═══════
auto wave_sol = differential_equations::wave_explicit_1d(u0, u0, 1.0, 0.02, 0.01, 500);

// ═══════ PDE: 2D LAPLACE EQUATION ═══════
auto laplace = differential_equations::laplace_2d(20, 20);
```

---

### 8. 📊 Statistics & Data Analysis (`numc::statistics`)

Comprehensive descriptive statistics, histograms, curve fitting, and error propagation.

```cpp
vector<double> data = {2.3, 4.5, 1.2, 6.7, 3.8, 5.1, 2.9};

// ═══════ DESCRIPTIVE STATISTICS ═══════
double m   = mean(data);            // Arithmetic mean
double med = median(data);          // Median
double s   = std_dev(data);         // Standard deviation
double v   = variance(data);        // Variance
double gm  = statistics::geometric_mean(data);   // Geometric mean
double hm  = statistics::harmonic_mean(data);     // Harmonic mean
double sk  = statistics::skewness(data);          // Skewness
double ku  = statistics::kurtosis(data);          // Kurtosis

// ═══════ REGRESSION ═══════
vector<double> xs = {1, 2, 3, 4, 5};
vector<double> ys = {2.1, 3.9, 6.2, 8.1, 9.8};
auto [slope, intercept] = linear_regression(xs, ys);
// slope ≈ 2.0, intercept ≈ 0.1

// ═══════ HISTOGRAMS (ROOT TH1D-style) ═══════
statistics::Histogram1D<double> hist(50, -5.0, 5.0);
for (double val : data) hist.fill(val);
hist.print();  // ASCII art histogram in the console!

// ═══════ CURVE FITTING (Levenberg-Marquardt) ═══════
using statistics::param;
param<double> A("A", 1.0);   // Initial guess for amplitude
param<double> B("B", 0.5);   // Initial guess for decay constant

auto model = A * exp(-B * funcv::x);
statistics::fit(model, xs, ys, {&A, &B});

std::cout << A << std::endl;  // Prints: A = <fitted_value> +/- <error>
std::cout << B << std::endl;  // Prints: B = <fitted_value> +/- <error>

// ═══════ ERROR PROPAGATION ═══════
tensor<double> cov_X = {{0.04, 0.01}, {0.01, 0.09}};
tensor<double> T_mat = {{2.0, 1.0}, {0.0, 3.0}};
tensor<double> cov_Y = statistics::propagate_error(cov_X, T_mat);
```

---

### 9. 🤖 Machine Learning (`numc::machine_learning`)

Train models in pure C++. No Python overhead, no TensorFlow, no PyTorch.

```cpp
// ═══════ NEURAL NETWORK (MLP with Adam optimizer) ═══════
using namespace machine_learning;

MLP<double> net({2, 16, 16, 1}, Activation::Sigmoid);  // 2 inputs → 16 → 16 → 1 output

// XOR training data
std::vector<vector<double>> X = {{0,0}, {0,1}, {1,0}, {1,1}};
std::vector<vector<double>> Y = {{0},   {1},   {1},   {0}};

net.train(X, Y, /*epochs=*/5000, /*lr=*/0.01);  // Adam optimizer

vector<double> prediction = net.forward({1.0, 0.0});  // → ~1.0
std::cout << "XOR(1,0) = " << prediction[0] << std::endl;

// ═══════ K-MEANS CLUSTERING ═══════
std::vector<vector<double>> points = {{1,2}, {1.5,1.8}, {5,8}, {8,8}, {1,0.6}, {9,11}};
auto result = kmeans(points, /*k=*/2);
// result.labels → cluster assignments
// result.centroids → cluster centers

// ═══════ PCA (Principal Component Analysis) ═══════
auto pca_result = pca(points, /*n_components=*/1);
// pca_result.explained_variance → variance per component
// pca_result.transformed_data → projected data

// ═══════ KNN (K-Nearest Neighbors) ═══════
auto label = knn_classify(points, labels, query_point, /*k=*/3);
```

---

### 10. 🕸️ Graph Algorithms (`numc::graphs`)

Weighted directed/undirected graphs with classic pathfinding and minimum spanning tree algorithms.

```cpp
using namespace graphs;

Graph<double> g(6);
g.add_edge(0, 1, 7.0);
g.add_edge(0, 2, 9.0);
g.add_edge(1, 2, 10.0);
g.add_edge(1, 3, 15.0);
g.add_edge(2, 3, 11.0);
g.add_edge(3, 4, 6.0);
g.add_edge(4, 5, 9.0);

// ═══════ SHORTEST PATHS ═══════
auto [dist_d, prev_d] = dijkstra(g, 0);        // Dijkstra's algorithm
auto [dist_b, prev_b] = bellman_ford(g, 0);     // Bellman-Ford (handles negatives)
auto dist_fw = floyd_warshall(g);                // All-pairs shortest paths

// ═══════ A* SEARCH ═══════
auto heuristic = [](size_t node) { return 0.0; };
auto [dist_a, path] = a_star(g, 0, 5, heuristic);

// ═══════ MINIMUM SPANNING TREE ═══════
auto mst_edges = kruskal_mst(g);   // Kruskal's MST
auto mst_edges2 = prim_mst(g);     // Prim's MST
```

---

### 11. 📐 Geometry (`numc::geometry`)

Computational geometry primitives: hulls, curves, and spatial queries.

```cpp
using namespace geometry;

// ═══════ CONVEX HULL ═══════
std::vector<Point<double>> pts = {{0,0}, {1,0}, {0.5,1}, {0,1}, {1,1}, {0.5,0.5}};
auto hull = convex_hull(pts);   // Returns hull vertices in CCW order

// ═══════ POINT IN POLYGON ═══════
bool inside = point_in_polygon({0.5, 0.5}, hull);  // true

// ═══════ BÉZIER CURVES ═══════
std::vector<Point<double>> control = {{0,0}, {1,3}, {3,3}, {4,0}};
Point<double> mid = bezier_curve(control, 0.5);  // Point at t = 0.5

// ═══════ B-SPLINE ═══════
std::vector<double> knots = {0, 0, 0, 0, 1, 1, 1, 1};
Point<double> spline_pt = b_spline(control, knots, 3, 0.5);

// ═══════ DISTANCES ═══════
double d = distance_point_to_segment({1,1}, {0,0}, {2,0});  // Distance to segment
```

---

### 12. 🔢 Number Theory (`numc::number_theory`)

Efficient algorithms for primes, modular arithmetic, and combinatorics.

```cpp
using namespace number_theory;

// ═══════ PRIMES ═══════
auto primes = sieve_of_eratosthenes(1000000);     // All primes up to 10⁶
bool is_p = is_prime_miller_rabin(999999937ULL);   // Deterministic Miller-Rabin

// ═══════ MODULAR ARITHMETIC ═══════
uint64_t result = mod_pow(2, 1000, 1000000007);   // 2^1000 mod 10⁹+7
uint64_t inv = mod_inverse(3, 1000000007);         // Modular inverse
uint64_t g = gcd(48, 18);                          // GCD = 6
uint64_t l = lcm(12, 18);                          // LCM = 36

// ═══════ EXTENDED EUCLIDEAN ═══════
auto [gcd_val, x_coeff, y_coeff] = extended_gcd(35, 15);
// 35*x + 15*y = gcd(35,15) = 5

// ═══════ CHINESE REMAINDER THEOREM ═══════
auto solution = chinese_remainder_theorem({2, 3, 2}, {3, 5, 7});
// Find x such that x ≡ 2(mod3), x ≡ 3(mod5), x ≡ 2(mod7) → x = 23

// ═══════ POLLARD'S RHO FACTORIZATION ═══════
auto factors = factorize(1234567890ULL);
```

---

### 13. 🌐 Vector Calculus (`numc::analysis`, Nabla ∇)

The `nabla` operator enables gradient, divergence, and curl computations on arbitrary functions.

```cpp
auto x = funcv::x;
auto y = funcv::y;
auto z = funcv::z;

// Scalar field
func<double> phi = sin(x) * cos(y) + z * z;

// ═══════ GRADIENT ═══════
auto grad = analysis::nabla(phi);   // Returns a vector_field<double>
vector<double> g = grad({1.0, 2.0, 3.0});  // ∇φ at (1, 2, 3)

// ═══════ DIVERGENCE ═══════
std::vector<func<double>> F = {x * y, y * z, z * x};  // Vector field
func<double> div_F = analysis::nabla * F;              // ∇·F (scalar)
double d = div_F(1.0, 2.0, 3.0);

// ═══════ CURL ═══════
auto curl_F = analysis::nabla ^ F;   // ∇×F (vector field)
vector<double> c = curl_F({1.0, 2.0, 3.0});

// ═══════ LAPLACIAN ═══════
func<double> laplacian = analysis::laplacian(phi);  // ∇²φ (scalar)
```

---

### 14. 🧪 Complex Numbers (`numc::complex`)

Full-featured complex number class with arithmetic, polar form, and complex analysis functions.

```cpp
complex<double> z1(3.0, 4.0);     // 3 + 4i
complex<double> z2 = complex<double>::i();  // Pure imaginary unit

// Arithmetic
auto z3 = z1 + z2;
auto z4 = z1 * z2;
auto z5 = z1 / z2;
auto z6 = z1 ^ 3.0;   // z1³

// Properties
double mag = z1.abs();     // |z| = 5.0
double arg = z1.arg();     // arg(z) ≈ 0.9273 rad
auto conj = z1.conjugate();
auto inv = z1.inverse();

// Complex analysis functions
auto z_exp = numc::cexp(z1);    // e^z
auto z_log = numc::clog(z1);    // ln(z)
auto z_sin = numc::csin(z1);    // sin(z)
auto z_sqrt = numc::csqrt(z1);  // √z
```

---

### 15. 📋 Data Handling

#### Polynomials

```cpp
polynomial<double> p = {1.0, -3.0, 2.0};  // p(x) = 1 - 3x + 2x²
double val = p(2.0);                        // p(2) = 1 - 6 + 8 = 3

// Convert to func for use in calculus
func<double> f = p.to_func();
double deriv = derivative(f, 1.0);

// Polynomial arithmetic
auto sum = p + polynomial<double>({0, 1});  // Add: 1 - 2x + 2x²
auto prod = p * polynomial<double>({1, 1}); // Multiply
```

#### Sparse Matrices

```cpp
sparse_matrix<double> S(1000, 1000, {
    {0, 0, 5.0}, {0, 1, -1.0},
    {1, 0, -1.0}, {1, 1, 5.0},
    // ... only non-zero entries
});

vector<double> b = {4.0, 4.0, /* ... */};
vector<double> x = S.solve_cg(b);   // Conjugate gradient solver
```

#### DataFrames (Pandas-style)

```cpp
auto df = DataFrame<double>::from_csv("measurements.csv");

vector<double>& temps = df["temperature"];
vector<double>& pressure = df["pressure"];

// Filter rows
auto filtered = df.filter([](const auto& row) {
    return row.at("temperature") > 25.0;
});
```

#### Interpolation

```cpp
vector<double> xd = {0, 1, 2, 3, 4};
vector<double> yd = {0, 0.84, 0.91, 0.14, -0.76};

// Lagrange interpolation
double y_interp = analysis::lagrange_interp(xd, yd, 2.5);

// Newton's divided differences
auto coeffs = analysis::newton_coeffts(xd, yd);
double y2 = analysis::newton_eval(coeffs, xd, 2.5);

// Natural cubic spline
auto spline = analysis::cubic_spline(xd, yd);
double y3 = analysis::spline_eval(spline, xd, 2.5);
```

#### Fourier Transform (FFT)

```cpp
vector<double> signal = {1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0};

// DFT
auto spectrum = analysis::dft(signal);

// Inverse DFT
vector<double> recovered = analysis::idft(spectrum);

// FFT (Cooley-Tukey, power-of-2 sizes)
auto fast_spectrum = analysis::fft(signal);
auto fast_recovered = analysis::ifft(fast_spectrum);
```

---

## 🎓 Real-World Example: Complete Analysis Pipeline

Here's a full example combining multiple modules — the kind of analysis you'd do in a physics lab, but entirely in compiled C++:

```cpp
#include <numc>
#include <iostream>
using namespace numc;

int main() {
    auto x = funcv::x;

    // 1. Define theoretical model: damped harmonic oscillation
    func<double> theory = 5.0 * exp(-0.15 * x) * cos(2.0 * x);

    // 2. "Experimental" data with noise
    vector<double> t_exp = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0};
    vector<double> y_exp = {5.0, 3.2, 0.8, -1.9, -1.5, 0.7, 1.1, -0.3, -0.5};

    // 3. Fit an exponentially damped cosine to the data
    statistics::param A("A", 5.0);
    statistics::param gamma("γ", 0.1);
    statistics::param omega("ω", 2.0);

    auto model = A * exp(-gamma * x) * cos(omega * x);
    statistics::fit(model, t_exp, y_exp, {&A, &gamma, &omega});

    std::cout << "Fitted parameters:" << std::endl;
    std::cout << "  " << A << std::endl;
    std::cout << "  " << gamma << std::endl;
    std::cout << "  " << omega << std::endl;

    // 4. Compute some analysis on the fitted model
    double energy_integral = model.integral(0.0, 10.0);
    double max_slope = derivative(model, 0.0);

    std::cout << "Total energy ∝ ∫f(t)dt = " << energy_integral << std::endl;
    std::cout << "Initial slope f'(0) = " << max_slope << std::endl;

    // 5. Visualize everything
    plot dashboard("Damped Oscillation — Experiment vs. Fit");
    dashboard.set_xlabel("Time [s]").set_ylabel("Displacement [cm]");
    dashboard.add(theory, 0.0, 8.0, "Theory (exact)");
    dashboard.add(model, 0.0, 8.0, "Fitted model");
    dashboard.add_scatter(t_exp, y_exp, "Experimental data");

    dashboard.show();                                    // Interactive window
    dashboard.show(PlotBackend::SVG, PlotTheme::LIGHT);  // Publication SVG
}
```

---

## 🛠️ Building

### Prerequisites
- **C++20** compatible compiler (GCC 10+, Clang 12+, MSVC 19.29+)
- **CMake 3.16** or higher

### Build Steps

```bash
git clone https://github.com/mslotwinski-dev/NumC.git
cd NumC
mkdir build && cd build
cmake ..
cmake --build .
```

### Running Examples

```bash
# On Windows
.\build\examples\example_plots.exe

# On Linux/macOS
./build/examples/example_plots
```

### Integration into Your Project

NumC is ~99% header-only. Just point your compiler's include path to the `src/` directory:

```cmake
add_subdirectory(NumC)
target_link_libraries(your_target PRIVATE numc)
```

Or simply copy the `src/` folder into your project and `#include "numc.h"`.

---

## 📚 Mathematical Constants

NumC provides high-precision (34+ significant digits) constants:

| Constant | Symbol | Value |
|----------|--------|-------|
| `PI` | π | 3.14159265358979… |
| `TWO_PI` | 2π / τ | 6.28318530717958… |
| `E` | e | 2.71828182845904… |
| `PHI` | φ | 1.61803398874989… |
| `SQRT2` | √2 | 1.41421356237309… |
| `SQRT3` | √3 | 1.73205080756887… |
| `EGAMMA` | γ | 0.57721566490153… |
| `LN2` | ln 2 | 0.69314718055994… |
| `LN10` | ln 10 | 2.30258509299404… |
| `CATALAN` | G | 0.91596559417721… |

---

## 📖 Module Quick Reference

| Module | Namespace | Key Features |
|--------|-----------|-------------|
| **Core Types** | `numc` | `vector`, `tensor`, `func`, `complex`, `polynomial`, `sparse_matrix`, `DataFrame` |
| **Elementary** | `numc` | `sin`, `cos`, `tan`, `exp`, `ln`, `log2`, `log10`, `sqrt`, `abs`, `floor`, `ceil` |
| **Special** | `numc` | `erf`, `tgamma`, `beta`, `riemann_zeta`, `expint`, elliptic integrals |
| **Orthogonal** | `numc` | `bessel_j/y/i/k`, `legendre`, `hermite`, `laguerre`, `chebyshev_t/u` |
| **Signal** | `numc` | `heaviside`, `dirac`, `ramp`, `sinc`, `rect`, `triangle`, `square_wave`, `sawtooth` |
| **Calculus** | `numc::analysis` | `derivative`, `integral`, `limit`, `taylor_coeffs`, `arc_length`, `curvature`, `tabulate` |
| **Integration** | `numc::analysis` | `simpson_1_3/3_8`, `romberg`, `gauss_legendre`, `adaptive_simpson`, `boole` |
| **Differentiation** | `numc::analysis` | `derivative_forward/backward/central/O4`, `derivative_second_central` |
| **Interpolation** | `numc::analysis` | `lagrange_interp`, `newton_coeffts`, `neville`, `cubic_spline` |
| **FFT** | `numc::analysis` | `dft`, `idft`, `fft`, `ifft` |
| **∇ (Nabla)** | `numc::analysis` | gradient, divergence, curl, Laplacian, Jacobian, Hessian |
| **LinAlg** | `numc::linear_algebra` | `gauss_elimin`, `lu_decomp`, `cholesky`, `jacobi`, `inverse`, `transpose`, `det` |
| **Root Finding** | `numc::optimalization` | `bisection`, `newton`, `halley`, `secant`, `regula_falsi`, `inverse_quadratic` |
| **Minimization** | `numc::optimalization` | `gold_bracket`, `gold_search` |
| **Least Squares** | `numc::optimalization` | `levenberg_marquardt` |
| **ODE** | `numc::differential_equations` | `taylor`, `run_kut4`, `run_kut5` (adaptive) |
| **BVP** | `numc::differential_equations` | `shooting_method` |
| **PDE** | `numc::differential_equations` | `heat_explicit_1d`, `heat_implicit_1d`, `wave_explicit_1d`, `laplace_2d` |
| **Statistics** | `numc::statistics` | `mean`, `median`, `variance`, `std_dev`, `skewness`, `kurtosis`, `correlation` |
| **Regression** | `numc::statistics` | `linear_regression`, `polynomial_regression` |
| **Histograms** | `numc::statistics` | `Histogram1D` (ROOT TH1D-style, with `fill`, `print`, `normalize`) |
| **Fitting** | `numc::statistics` | `fit` (Levenberg-Marquardt with `param` objects) |
| **Distributions** | `numc::statistics` | `dist` (base class), custom distributions with normalization & sampling |
| **Errors** | `numc::statistics` | `propagate_error`, `extract_uncertainties` |
| **Neural Nets** | `numc::machine_learning` | `MLP`, `DenseLayer`, activations (Sigmoid/ReLU/Tanh), Adam optimizer |
| **Clustering** | `numc::machine_learning` | `kmeans` |
| **PCA** | `numc::machine_learning` | `pca` |
| **KNN** | `numc::machine_learning` | `knn_classify`, `knn_regress` |
| **Decision Tree** | `numc::machine_learning` | `decision_tree_classify` |
| **Graphs** | `numc::graphs` | `dijkstra`, `bellman_ford`, `a_star`, `bfs_shortest`, `floyd_warshall` |
| **MST** | `numc::graphs` | `kruskal_mst`, `prim_mst` |
| **Geometry** | `numc::geometry` | `convex_hull`, `point_in_polygon`, `bezier_curve`, `b_spline`, `distance_point_to_segment` |
| **Number Theory** | `numc::number_theory` | `sieve_of_eratosthenes`, `is_prime_miller_rabin`, `mod_pow`, `gcd`, `lcm`, `extended_gcd`, `mod_inverse`, `chinese_remainder_theorem`, `factorize` |
| **Visualization** | `numc` | `plot`, `plt::show`, `PlotBackend`, `PlotTheme` |
| **I/O** | `numc::utility` | `read_data_file`, `unpack_columns` |

---

<div align="center">

### *"If it can be written as math, it can be written in NumC."*

Built with ❤️ for elegant C++ and beautiful mathematics.

</div>
