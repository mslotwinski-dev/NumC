#include "../src/numc.h"
#include <iostream>
#include <vector>
#include <tuple>

using namespace numc;

int main() {
  std::cout << "=== Sparse Matrix Example ===" << std::endl;

  // A 3x3 sparse matrix
  // [ 4.0  -1.0   0.0 ]
  // [-1.0   4.0  -1.0 ]
  // [ 0.0  -1.0   3.0 ]
  // A 3x3 sparse matrix using clean initializer syntax
  sparse_matrix<double> A(3, 3, {
    {0, 0, 4.0}, {0, 1, -1.0},
    {1, 0, -1.0}, {1, 1, 4.0}, {1, 2, -1.0},
    {2, 1, -1.0}, {2, 2, 3.0}
  });
  std::cout << A << std::endl;

  vector<double> b = {1.0, 2.0, 3.0};
  std::cout << "Vector b: " << b[0] << ", " << b[1] << ", " << b[2] << std::endl;

  // Solve Ax = b using Conjugate Gradient
  vector<double> x_cg = A.solve_cg(b);
  std::cout << "Solution x (CG): " << x_cg[0] << ", " << x_cg[1] << ", " << x_cg[2] << std::endl;

  // Solve using Gauss-Seidel
  vector<double> x_gs = A.solve_gauss_seidel(b);
  std::cout << "Solution x (Gauss-Seidel): " << x_gs[0] << ", " << x_gs[1] << ", " << x_gs[2] << std::endl;

  return 0;
}
