#include "../src/numc.h"
#include <iostream>
#include <iomanip>

using namespace numc;
using namespace numc::analysis;

int main() {
  std::cout << "=== Autodiff Example ===" << std::endl;

  // 1D Derivative using Calculus module (since func<double> is evaluated numerically)
  // f(x) = sin(x) * x^2
  auto f_1d = sin(x) * (x ^ 2.0);

  double x_val = 2.0;
  double deriv = derivative_central_O4(f_1d, x_val);
  
  std::cout << "f(x) = sin(x) * x^2 at x = " << x_val << std::endl;
  std::cout << "f'(x) = " << deriv << std::endl;

  // ND Mathematical Gradient
  std::cout << "\n--- Vector Calculus with Nabla ---" << std::endl;
  
  // Declarative ND function: f(x, y) = sin(x) * cos(y)
  func<double> f_nd = sin(x) * cos(y);
  
  // Calculate the gradient vector field using the nabla operator
  auto grad_field = nabla(f_nd);
  
  vector<double> point = {3.14159 / 4.0, 3.14159 / 4.0};
  
  // Evaluate the gradient field at the point
  vector<double> grad = grad_field(point);
  
  std::cout << "f(x, y) = sin(x) * cos(y)" << std::endl;
  std::cout << "Gradient at (pi/4, pi/4) = [" << grad[0] << ", " << grad[1] << "]" << std::endl;

  // Numerical Hessian
  std::cout << "\n--- Numerical Hessian ---" << std::endl;
  tensor<double> H = hessian(f_nd, point);
  std::cout << "Hessian Matrix:" << std::endl;
  std::cout << H << std::endl;

  return 0;
}
