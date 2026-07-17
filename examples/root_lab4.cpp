#include "../src/numc.h"
#include <iostream>
#include <vector>
#include <cmath>

using namespace numc;

int main() {
  std::cout << "========== ROOT Lab 4 Replication (Beautiful API) ==========\n";
  
  std::string filename = "dane_lab4.txt";
  std::vector<std::vector<double>> data;
  try {
    data = utility::read_data_file<double>(filename);
  } catch(const std::exception& e) {
    std::cerr << "Pomiń: " << e.what() << "\n(Należy wygenerować plik dane_lab4.txt)\n";
    return 1;
  }

  if (data.empty()) return 1;

  vector<double> x1(data.size()), x2(data.size()), x3(data.size());
  for (size_t i = 0; i < data.size(); ++i) {
    if (data[i].size() >= 3) {
      x1[i] = data[i][0];
      x2[i] = data[i][1];
      x3[i] = data[i][2];
    }
  }

  // Macierz kowariancji
  std::vector<vector<double>> dataset = {x1, x2, x3};
  auto cov_X = statistics::cov_matrix(dataset);
  std::cout << "Macierz kowariancji X (X1, X2, X3):\n" << cov_X << "\n";

  double e_x1 = statistics::mean(x1);
  double e_x2 = statistics::mean(x2);
  double e_x3 = statistics::mean(x3);
  std::cout << "E(X1) = " << e_x1 << ", E(X2) = " << e_x2 << ", E(X3) = " << e_x3 << "\n";

  double expected_y1 = 21.9967;
  double expected_y2 = -3.77383;

  double a_approx = (expected_y2 - 0.5 * e_x1 + 3.0 * e_x3) / e_x2;
  int a = std::round(a_approx);
  double b_approx = expected_y1 - a * e_x3;
  int b = std::round(b_approx);

  std::cout << "Znalezione stałe: a = " << a << ", b = " << b << "\n";

  // Transformation Matrix T
  tensor<double> T({2, 3});
  T(0, 0) = 0.0;  T(0, 1) = 0.0;          T(0, 2) = static_cast<double>(a);
  T(1, 0) = 0.5;  T(1, 1) = static_cast<double>(a); T(1, 2) = -3.0;

  auto cov_Y = statistics::propagate_error(cov_X, T);
  std::cout << "\nMacierz kowariancji Y (Y1, Y2):\n" << cov_Y << "\n";

  auto u_X = statistics::extract_uncertainties(cov_X);
  auto u_Y = statistics::extract_uncertainties(cov_Y);

  std::cout << "\nNiepewności u(X):\n";
  std::cout << "u(X1) = " << u_X[0] << "\n";
  std::cout << "u(X2) = " << u_X[1] << "\n";
  std::cout << "u(X3) = " << u_X[2] << "\n";

  std::cout << "\nNiepewności u(Y):\n";
  std::cout << "u(Y1) = " << u_Y[0] << "\n";
  std::cout << "u(Y2) = " << u_Y[1] << "\n";

  // WIZUALIZACJA: relacje pomiędzy zmiennymi z datasetu
  plot dashboard("ROOT Lab 4: Relacje X1-X2 oraz X2-X3");
  dashboard.set_xlabel("X1 / X2").set_ylabel("X2 / X3");
  
  dashboard.add_scatter(x1, x2, "Zależność X1 vs X2");
  dashboard.add_scatter(x2, x3, "Zależność X2 vs X3");
  
  dashboard.show();

  return 0;
}
