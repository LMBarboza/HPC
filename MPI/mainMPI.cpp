#include <vector>
#include <iostream>
#include <mpi.h>

#define MASTER_RANK 0

#define MASTER_TAG 1
#define WORKER_TAG 2

#define M_SIZE 100

MPI_Status status;

std::vector<std::vector<float>> aMatrix(rows, std::vector<float>(cols, 0.0f));
std::vector<std::vector<float>> bMatrix(rows, std::vector<float>(cols, 0.0f));
std::vector<std::vector<float>> cMatrix(rows, std::vector<float>(cols, 0.0f));


std::vector<std::vector<float>> create_matrix(int rows, int cols){
  std::vector<std::vector<float>> data(rows, std::vector<float>(cols, 0.0f));
  for (size_t i = 0; i < rows; ++i) {
      for (size_t j = 0; j < cols; ++j) {
          data[i][j] = static_cast<float>(rand()) / RAND_MAX;
        }
    }
  return data;
}

void matrix_multiplication(void){

}

void main(void){
  

}

