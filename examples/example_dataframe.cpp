#include "../src/numc.h"
#include <iostream>

using namespace numc;

int main() {
  std::cout << "=== DataFrame Example ===" << std::endl;

  DataFrame<double> df;
  df.add_column("x", {1.0, 2.0, 3.0, 4.0, 5.0});
  df.add_column("y", {2.0, 4.0, 6.0, 8.0, 10.0});
  df.add_column("z", {1.5, 2.5, 3.5, 4.5, 5.5});

  std::cout << "DataFrame has " << df.get_row_count() << " rows." << std::endl;
  
  std::cout << "Column 'y' values: ";
  for (double val : df["y"]) {
    std::cout << val << " ";
  }
  std::cout << std::endl;

  // Declarative Filtering using boolean masks
  auto filtered_df = df.filter(df["x"] > 2.0);

  std::cout << "Filtered DataFrame (x > 2.0) has " << filtered_df.get_row_count() << " rows." << std::endl;

  std::cout << "\n--- utility::unpack_columns ---" << std::endl;
  std::vector<std::vector<double>> raw_data = {
    {1.0, 0.1, 10.0},
    {2.0, 0.2, 20.0},
    {3.0, 0.3, 30.0}
  };

  // Automatically unpack matrix rows into columnar vectors
  auto [col1, col2, col3] = utility::unpack_columns<3>(raw_data);
  
  std::cout << "Column 1: " << col1[0] << ", " << col1[1] << ", " << col1[2] << std::endl;
  std::cout << "Column 2: " << col2[0] << ", " << col2[1] << ", " << col2[2] << std::endl;
  std::cout << "Column 3: " << col3[0] << ", " << col3[1] << ", " << col3[2] << std::endl;

  return 0;
}
