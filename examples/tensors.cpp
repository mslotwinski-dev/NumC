#include <iostream>
#include <numc>

using namespace numc;

int main() {
  tensor A({2, 2}, {4.0, 7.0, 2.0, 6.0});
  tensor B({2, 2}, {1.0, 2.0, 3.0, 4.0});
  tensor Z({2, 3});
  tensor I = tensor<>::eye(2);

  Z(0, 1) = 5.0;
  Z(1, 2) = 9.0;

  std::cout << "A:\n" << A << "\n\n";
  std::cout << "B:\n" << B << "\n\n";
  std::cout << "Z:\n" << Z << "\n\n";
  std::cout << "I:\n" << I << "\n\n";

  std::cout << "A + B =\n" << A + B << "\n\n";
  std::cout << "A - B =\n" << A - B << "\n\n";
  std::cout << "2A =\n" << 2 * A << "\n\n";

  std::cout << "A * B =\n" << A * B << "\n\n";
  std::cout << "Mnozenie Hadamarda =\n" << A.hadamard(B) << "\n\n";
  std::cout << "Macierz odwrotna (A^-1) =\n" << 1 / A << "\n\n";
  std::cout << "Dowod: A * A^-1 (Powinno dac macierz jednostkowa):\n" << A * (1 / A) << "\n\n";

  std::cout << "Wymiary Z: " << Z.dimensions() << "D\n";
  std::cout << "Ksztalt Z: (" << Z.shape()[0] << ", " << Z.shape()[1] << ")\n";
  std::cout << "Laczna liczba elementow Z: " << Z.total_size() << "\n";

  return 0;
}