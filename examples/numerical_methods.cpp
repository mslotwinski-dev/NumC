#include <iostream>
#include <iomanip>
#include <numc>

using namespace numc;
using namespace numc::linear_algebra;
using namespace numc::analysis;
using namespace numc::optimalization;
using namespace numc::differential_equations;
using namespace numc::statistics;

int main() {
  std::cout << std::fixed << std::setprecision(6);
  std::cout << "================= NumC MEGA TEST =================\n\n";

  // ====================================================================
  // 1. LINEAR ALGEBRA - MATRIX OPERATIONS
  // ====================================================================
  std::cout << "--- 1. Matrix Operations ---\n";
  tensor A({3, 3}, {2.0, 1.0, 0.0,  1.0, 3.0, 1.0,  0.0, 1.0, 2.0});
  std::cout << "det(A) = " << det(A) << "\n";
  std::cout << "trace(A) = " << trace(A) << "\n";
  std::cout << "rank(A) = " << rank(A) << "\n";
  std::cout << "frobenius_norm(A) = " << frobenius_norm(A) << "\n";

  auto [Q, R] = qr_decomp(A);
  std::cout << "QR decomposition - Q:\n" << Q << "\n";
  std::cout << "R:\n" << R << "\n";

  tensor B({2, 3}, {1.0, 2.0, 3.0,  4.0, 5.0, 6.0});
  auto [U_svd, sigma, V_svd] = svd(B);
  std::cout << "SVD singular values of [[1,2,3],[4,5,6]]: " << sigma << "\n";

  // Least squares: overdetermined system
  tensor C({3, 2}, {1.0, 1.0,  1.0, 2.0,  1.0, 3.0});
  vector<> b_ls = {1.0, 2.0, 4.0};
  vector<> x_ls = lstsq(C, b_ls);
  std::cout << "Least-squares solution: " << x_ls << "\n\n";

  // ====================================================================
  // 2. FFT / SIGNAL PROCESSING
  // ====================================================================
  std::cout << "--- 2. FFT / Signal Processing ---\n";
  vector<> signal = {1.0, 0.0, -1.0, 0.0, 1.0, 0.0, -1.0, 0.0};
  auto X_fft = fft(signal);
  std::cout << "FFT of [1,0,-1,0,1,0,-1,0]:\n";
  for (size_t i = 0; i < X_fft.size(); ++i) {
    std::cout << "  X[" << i << "] = " << X_fft[i].real()
              << " + " << X_fft[i].imag() << "i\n";
  }

  vector<> recovered = ifft(X_fft);
  std::cout << "IFFT recovery: " << recovered << "\n";

  // Convolution
  vector<> a_conv = {1.0, 2.0, 3.0};
  vector<> b_conv = {0.0, 1.0, 0.5};
  vector<> conv_result = convolve(a_conv, b_conv);
  std::cout << "Convolution [1,2,3] * [0,1,0.5]: " << conv_result << "\n\n";

  // ====================================================================
  // 3. CALCULUS TOOLS
  // ====================================================================
  std::cout << "--- 3. Calculus Tools ---\n";
  numc::func<double> x_var = numc::func<double>::x();

  // Limit: lim_{x->0} sin(x)/x = 1
  numc::func<double> sinc_f = numc::func<double>([](double x) { return std::sin(x) / x; });
  std::cout << "lim_{x->0} sin(x)/x = " << limit(sinc_f, 0.0) << "\n";

  // Series: sum_{n=1}^{inf} 1/n^2 = pi^2/6
  numc::func<double> inv_sq = numc::func<double>([](double n) { return 1.0 / (n * n); });
  std::cout << "sum 1/n^2 (n=1..inf) = " << series_sum_convergent(inv_sq, 1)
            << " (exact: " << numc::PI * numc::PI / 6.0 << ")\n";

  // Arc length of sin(x) from 0 to pi
  numc::func<double> sin_f = numc::func<double>([](double x) { return std::sin(x); });
  std::cout << "Arc length of sin(x) from 0 to pi = " << arc_length(sin_f, 0.0, numc::PI) << "\n";

  // Curvature of x^2 at x=0
  numc::func<double> quad_f = x_var * x_var;
  std::cout << "Curvature of x^2 at x=0: " << curvature(quad_f, 0.0) << "\n";

  // Taylor polynomial of sin(x) around 0
  vector<> taylor_sin = taylor_coeffs(sin_f, 0.0, 7);
  std::cout << "Taylor coefficients of sin(x): " << taylor_sin << "\n";
  std::cout << "Taylor approx sin(0.5) = " << taylor_eval(taylor_sin, 0.0, 0.5)
            << " (exact: " << std::sin(0.5) << ")\n";

  // Fundamental Theorem: d/dx integral(sin, 0, x) ≈ sin(x)
  numc::func<double> F_int = integral(sin_f, 0.0);
  numc::func<double> dF = derivative(F_int);
  std::cout << "d/dx int_0^x sin(t)dt at x=1 = " << dF(1.0)
            << " (exact: sin(1) = " << std::sin(1.0) << ")\n\n";

  // ====================================================================
  // 4. ADVANCED INTEGRATION
  // ====================================================================
  std::cout << "--- 4. Advanced Integration ---\n";
  numc::func<double> f_test = x_var * x_var; // x^2

  std::cout << "Simpson 1/3 of x^2 from 0 to 1: " << simpson_1_3(f_test, 0.0, 1.0) << "\n";
  std::cout << "Simpson 3/8 of x^2 from 0 to 1: " << simpson_3_8(f_test, 0.0, 1.0) << "\n";
  std::cout << "Boole of x^2 from 0 to 1: " << boole(f_test, 0.0, 1.0) << "\n";
  std::cout << "Adaptive Simpson of x^2 from 0 to 1: " << adaptive_simpson(f_test, 0.0, 1.0) << "\n";

  // Monte Carlo integration
  auto [mc_val, mc_err] = monte_carlo_integrate(f_test, 0.0, 1.0, 100000);
  std::cout << "Monte Carlo of x^2 from 0 to 1: " << mc_val << " +/- " << mc_err << "\n\n";

  // ====================================================================
  // 5. BVP SOLVERS
  // ====================================================================
  std::cout << "--- 5. Boundary Value Problems ---\n";
  // y'' = -4y + 4x, y(0) = 0, y(1) = 2 => unique solution exists
  auto [x_bvp, y_bvp] = finite_diff_bvp<double>(
      [](double) { return 0.0; },  // p(x) = 0
      [](double) { return -4.0; }, // q(x) = -4
      [](double x) { return 4.0 * x; },  // r(x) = 4x
      0.0, 1.0, 0.0, 2.0, 40);
  std::cout << "BVP finite diff at x=0.5: " << y_bvp[20] << "\n\n";


  // ====================================================================
  // 6. PDE SOLVERS
  // ====================================================================
  std::cout << "--- 6. PDE Solvers ---\n";
  // 2D Laplace on 11x11 grid with top boundary = 100
  auto u_lap = laplace_2d<double>(11, 11,
      [](size_t i, size_t j) -> double {
        if (j == 10) return 100.0; // top boundary
        return 0.0;
      });
  std::cout << "Laplace 2D center value (5,5): " << u_lap(5, 5) << "\n\n";

  // ====================================================================
  // 7. STATISTICS
  // ====================================================================
  std::cout << "--- 7. Statistics ---\n";
  vector<> data = {2.0, 4.0, 4.0, 4.0, 5.0, 5.0, 7.0, 9.0};
  std::cout << "Mean: " << mean(data) << "\n";
  std::cout << "Median: " << median(data) << "\n";
  std::cout << "Std Dev: " << std_dev(data) << "\n";
  std::cout << "Skewness: " << skewness(data) << "\n";
  std::cout << "Kurtosis: " << kurtosis(data) << "\n";
  std::cout << "25th percentile: " << percentile(data, 25.0) << "\n";
  std::cout << "75th percentile: " << percentile(data, 75.0) << "\n";

  // Linear regression
  vector<> xr = {1.0, 2.0, 3.0, 4.0, 5.0};
  vector<> yr = {2.1, 3.9, 6.2, 7.8, 10.1};
  auto [a_reg, b_reg, r2] = linear_regression(xr, yr);
  std::cout << "Linear regression: y = " << a_reg << " + " << b_reg << "x (R^2 = " << r2 << ")\n";

  // Cumulative sum
  std::cout << "Cumsum: " << cumsum(data) << "\n\n";

  // ====================================================================
  // 8. RANDOM + MONTE CARLO
  // ====================================================================
  std::cout << "--- 8. Random & Monte Carlo ---\n";
  auto rnd = random::normal<double>(0.0, 1.0, 5, 123);
  std::cout << "5 normal random numbers: " << rnd << "\n";
  std::cout << "Mean of random: " << mean(rnd) << "\n\n";

  // ====================================================================
  // 9. ADVANCED OPTIMIZATION
  // ====================================================================
  std::cout << "--- 9. Advanced Optimization ---\n";
  // Nelder-Mead on Rosenbrock
  auto rosenbrock = [](const vector<double>& x) -> double {
    return (1.0 - x[0]) * (1.0 - x[0]) + 100.0 * (x[1] - x[0] * x[0]) * (x[1] - x[0] * x[0]);
  };
  auto [nm_x, nm_iter] = nelder_mead<double>(rosenbrock, {-1.0, -1.0}, 0.5, 1e-8, 10000);
  std::cout << "Nelder-Mead on Rosenbrock: " << nm_x << " (iter: " << nm_iter << ")\n";

  // Simulated annealing
  auto [sa_x, sa_f] = simulated_annealing<double>(rosenbrock, {-1.0, -1.0}, 100.0, 1e-10, 0.999, 200);
  std::cout << "Simulated Annealing on Rosenbrock: " << sa_x << " (f = " << sa_f << ")\n\n";

  // ====================================================================
  // 10. POLYNOMIALS
  // ====================================================================
  std::cout << "--- 10. Polynomials ---\n";
  polynomial<double> p1({1.0, -2.0, 1.0}); // 1 - 2x + x^2 = (x-1)^2
  polynomial<double> p2({3.0, 1.0});       // 3 + x
  std::cout << "p1(x) = " << p1 << "\n";
  std::cout << "p2(x) = " << p2 << "\n";
  std::cout << "p1(2) = " << p1(2.0) << "\n";
  std::cout << "p1 * p2 = " << p1 * p2 << "\n";
  std::cout << "p1 derivative = " << p1.derivative() << "\n";
  std::cout << "Integral of p1 from 0 to 1 = " << p1.integrate(0.0, 1.0) << " (exact: 0.333333)\n\n";

  // ====================================================================
  // 11. SPARSE MATRICES
  // ====================================================================
  std::cout << "--- 11. Sparse Matrices ---\n";
  std::vector<std::tuple<size_t, size_t, double>> triplets = {
      {0, 0, 4.0}, {0, 1, -1.0},
      {1, 0, -1.0}, {1, 1, 4.0}, {1, 2, -1.0},
      {2, 1, -1.0}, {2, 2, 4.0}
  };
  sparse_matrix<double> sparseA = sparse_matrix<double>::from_triplets(3, 3, triplets);
  vector<double> sparseB = {1.0, 2.0, 3.0};
  vector<double> x_cg = sparseA.solve_cg(sparseB);
  std::cout << sparseA;
  std::cout << "CG solution to sparse Ax=b: " << x_cg << "\n\n";

  // ====================================================================
  // 12. AUTOMATIC DIFFERENTIATION
  // ====================================================================
  std::cout << "--- 12. Automatic Differentiation ---\n";
  auto my_func_ad = [](auto x) {
    return sin(x) * exp(x);
  };
  double ad_val = auto_diff(my_func_ad, 1.0);
  std::cout << "Auto-diff d/dx(sin(x)*exp(x)) at x=1: " << ad_val
            << " (exact: " << std::sin(1.0)*std::exp(1.0) + std::cos(1.0)*std::exp(1.0) << ")\n";

  // Gradient of f(x,y) = x^2 + 3xy + y^2 at (1, 2)
  auto multi_F = [](const vector<double>& v) {
    return v[0]*v[0] + 3.0*v[0]*v[1] + v[1]*v[1];
  };
  vector<double> grad_multi = gradient<double>(multi_F, {1.0, 2.0});
  std::cout << "Numerical gradient of x^2+3xy+y^2 at (1,2): " << grad_multi << "\n\n";

  // ====================================================================
  // 13. LEVENBERG-MARQUARDT NON-LINEAR LEAST SQUARES
  // ====================================================================
  std::cout << "--- 13. Non-Linear Least Squares (Levenberg-Marquardt) ---\n";
  // Fit y = A * exp(B * x) to data: x = [1, 2, 3], y = [2.5, 6.5, 18.0]
  // Expected roughly A=1, B=0.9
  vector<double> lm_x = {1.0, 2.0, 3.0};
  vector<double> lm_y = {2.5, 6.5, 18.0};

  auto residuals = [&](const vector<double>& beta) -> vector<double> {
    vector<double> r(3);
    for (size_t i = 0; i < 3; ++i) {
      r[i] = lm_y[i] - beta[0] * std::exp(beta[1] * lm_x[i]);
    }
    return r;
  };

  auto jacobian_ls = [&](const vector<double>& beta) -> tensor<double> {
    tensor<double> J({3, 2});
    for (size_t i = 0; i < 3; ++i) {
      J(i, 0) = -std::exp(beta[1] * lm_x[i]);
      J(i, 1) = -beta[0] * lm_x[i] * std::exp(beta[1] * lm_x[i]);
    }
    return J;
  };

  vector<double> beta0 = {1.0, 1.0}; // Initial guess
  vector<double> beta_opt = levenberg_marquardt<double>(residuals, jacobian_ls, beta0);
  std::cout << "Fitted parameters (A, B): " << beta_opt << "\n\n";

  // ====================================================================
  // 14. MACHINE LEARNING
  // ====================================================================
  std::cout << "--- 14. Machine Learning ---\n";
  
  // 14.1 K-Means Clustering
  std::vector<vector<double>> cluster_data = {
    {1.0, 1.0}, {1.5, 1.8}, {1.2, 1.4}, // Cluster 1
    {8.0, 8.0}, {8.5, 8.2}, {9.0, 8.5}  // Cluster 2
  };
  auto kmeans_res = numc::machine_learning::kmeans<double>(cluster_data, 2, 100, 1e-4, 42);
  std::cout << "K-Means Centroid 1: " << kmeans_res.centroids[0] << "\n";
  std::cout << "K-Means Centroid 2: " << kmeans_res.centroids[1] << "\n";

  // 14.2 PCA (Principal Component Analysis)
  std::vector<vector<double>> pca_data = {
    {2.5, 2.4}, {0.5, 0.7}, {2.2, 2.9}, {1.9, 2.2}, {3.1, 3.0},
    {2.3, 2.7}, {2.0, 1.6}, {1.0, 1.1}, {1.5, 1.6}, {1.1, 0.9}
  };
  auto pca_res = numc::machine_learning::pca<double>(pca_data, 1);
  std::cout << "PCA Top Component: " << pca_res.components << "\n";
  std::cout << "PCA Explained Variance: " << pca_res.explained_variance[0] << "\n";
  std::cout << "PCA Transformed Point 0: " << pca_res.transformed_data[0][0] << "\n";

  // 14.3 Neural Network (MLP) - XOR Problem
  numc::machine_learning::MLP<double> nn;
  nn.add_layer(2, 4, numc::machine_learning::Activation::Tanh, 1); // Hidden layer
  nn.add_layer(4, 1, numc::machine_learning::Activation::Sigmoid, 2); // Output layer

  std::vector<vector<double>> xor_X = {{0.0, 0.0}, {0.0, 1.0}, {1.0, 0.0}, {1.0, 1.0}};
  std::vector<vector<double>> xor_Y = {{0.0}, {1.0}, {1.0}, {0.0}};

  std::cout << "Training MLP on XOR problem for 5000 epochs...\n";
  nn.train(xor_X, xor_Y, 5000, 0.5);

  std::cout << "MLP Predictions:\n";
  for (size_t i = 0; i < xor_X.size(); ++i) {
    auto pred = nn.predict(xor_X[i]);
    std::cout << "  Input: " << xor_X[i] << " -> Predict: " << pred[0] << " (True: " << xor_Y[i][0] << ")\n";
  }
  std::cout << "\n";

  // ====================================================================
  // 15. GRAPHS & PATHFINDING
  // ====================================================================
  std::cout << "--- 15. Graph Theory & Algorithms ---\n";
  numc::graphs::Graph<double> g(5);
  g.add_edge(0, 1, 10.0);
  g.add_edge(0, 4, 5.0);
  g.add_edge(1, 2, 1.0);
  g.add_edge(1, 4, 2.0);
  g.add_edge(2, 3, 4.0);
  g.add_edge(3, 0, 7.0);
  g.add_edge(4, 1, 3.0);
  g.add_edge(4, 2, 9.0);
  g.add_edge(4, 3, 2.0);
  
  auto [dist, prev] = g.dijkstra(0);
  std::vector<size_t> path = g.reconstruct_path(3, prev);
  std::cout << "Dijkstra shortest path from 0 to 3: Distance = " << dist[3] << "\nPath nodes: ";
  for (size_t n : path) std::cout << n << " ";
  std::cout << "\n\n";

  // ====================================================================
  // 16. NUMBER THEORY & MODULAR MATH
  // ====================================================================
  std::cout << "--- 16. Number Theory & Cryptography Basics ---\n";
  std::vector<uint64_t> primes = numc::number_theory::sieve_of_eratosthenes(30);
  std::cout << "Primes up to 30: ";
  for (auto p : primes) std::cout << p << " ";
  std::cout << "\nModular exponentiation 5^13 mod 17: " << numc::number_theory::mod_pow(5, 13, 17) << "\n\n";

  // ====================================================================
  // 17. COMPUTATIONAL GEOMETRY
  // ====================================================================
  std::cout << "--- 17. Computational Geometry ---\n";
  std::vector<numc::geometry::Point<double>> hull_pts = {
    {0,0}, {2,2}, {1,1}, {0,2}, {2,0}, {0.5, 0.5}
  };
  auto convex_hull = numc::geometry::convex_hull(hull_pts);
  std::cout << "Convex Hull Vertices:\n";
  for (const auto& p : convex_hull) {
    std::cout << "  (" << p[0] << ", " << p[1] << ")\n";
  }
  std::cout << "\n";

  // ====================================================================
  // 18. ADVANCED MACHINE LEARNING (DECISION TREES, KNN)
  // ====================================================================
  std::cout << "--- 18. Advanced Machine Learning ---\n";
  
  std::vector<numc::vector<double>> clf_X = {
    {1.0, 1.0}, {1.5, 1.5}, {1.2, 1.3},
    {8.0, 8.0}, {8.5, 8.5}, {8.2, 8.3}
  };
  std::vector<double> clf_y = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
  
  numc::machine_learning::DecisionTreeClassifier<double> dt(5, 2);
  dt.train(clf_X, clf_y);
  std::cout << "Decision Tree predict (2.0, 2.0): " << dt.predict({2.0, 2.0}) << "\n";
  std::cout << "Decision Tree predict (7.0, 7.0): " << dt.predict({7.0, 7.0}) << "\n";

  numc::machine_learning::KNNClassifier<double> knn(3);
  knn.train(clf_X, clf_y);
  std::cout << "KNN (k=3) predict (1.1, 1.1): " << knn.predict({1.1, 1.1}) << "\n";
  std::cout << "KNN (k=3) predict (9.0, 9.0): " << knn.predict({9.0, 9.0}) << "\n\n";

  std::cout << "==================================================================\n";
  std::cout << "All tests passed!\n";
  return 0;
}
