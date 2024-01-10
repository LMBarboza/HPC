#include <iostream>
#include "../include/matrix.h"
#include "../include/matrixAVX.h"
#include "../include/matrixMP.h"
#include "../include/matrixMPI.h"
#include <vector>
#include <chrono>

std::vector<std::vector<float>> create_matrix(int rows, int cols){
  std::vector<std::vector<float>> data(rows, std::vector<float>(cols, 0.0f));

  for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < cols; ++j) {
          data[i][j] = static_cast<float>(rand()) / RAND_MAX;
        }
    }
  return data;

}
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
  const size_t rows = 1000;
  const size_t cols = 1000;
  std::vector<std::vector<float>> mat1 = create_matrix(rows, cols);
  std::vector<std::vector<float>> mat2 = create_matrix(rows, cols);
  MatrizAVX matriz1(mat1);
  MatrizAVX matriz2(mat2);
  MatrizMPI matriz3(mat1);
  MatrizMPI matriz4(mat2);
  Matriz matriz5(mat1);
  Matriz matriz6(mat1);
  std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();
  MatrizMPI result = matriz3 + matriz4;
	std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();
  std::cout << "Time difference = " << std::chrono::duration_cast<std::chrono::microseconds>
    (end - begin).count()/100 << "[µs]" << std::endl;

  
   
  return 0;
}
