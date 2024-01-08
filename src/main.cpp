#include <iostream>
#include "../include/matrixAVX.h"
#include <vector>

int main(void){
  /*std::cout << "Insira o número de linhas e colunas: ";
  size_t lin, col;
  std::cin >> lin >> col;

  std::vector<std::vector<float>> dados(lin, std::vector<float>(col, 0.0));

  std::cout << "Insira a matriz linha a linha:\n";
    for (size_t i = 0; i < lin; ++i) {
      for (size_t j = 0; j < col; ++j) {
        std::cout << "Matriz[" << i << "][" << j << "]: ";
        std::cin >> dados[i][j];
      }
    }*/
  const size_t rows = 2;
  const size_t cols = 2;

  std::vector<std::vector<float>> data(rows, std::vector<float>(cols, 0.0f));

    // Fill the matrix with random values
  for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < cols; ++j) {
          data[i][j] = static_cast<float>(rand()) / RAND_MAX;
        }
    }
  MatrizAVX matriz1(data);
  matriz1.print();
  MatrizAVX matriz2(data);
  MatrizAVX result1 = matriz1 + matriz2;
  result1.print();
  MatrizAVX result2 = matriz1 - matriz2;
  result2.print();
  MatrizAVX result3 = matriz1 * 2;
  result3.print();
  MatrizAVX result4 = matriz1 / matriz2;
  result4.print();

  return 0;
}
