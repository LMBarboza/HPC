#include <vector>
#include <iostream>
#include <mpi.h>
#include <iterator>

void create_matrix(float data[], int rows, int cols){
  for (int i = 0; i < rows * cols; ++i) {
          data[i] = static_cast<float>(rand()) / RAND_MAX;
    }
}
void create_zero_matrix(float data[], int rows, int cols){
  for (int i = 0; i < rows * cols; ++i) {
          data[i] = 0;
    }
}


void print_matrix(float data[], int rows, int cols){
  for (int i = 0; i < rows * cols; ++i) {
    std::cout << data[i] << " ";
    }


}
 

int main(void){
  MPI_Status status;
  int rank, numProcesses;
  double start, end;
  MPI_Init(NULL, NULL);

  MPI_Barrier(MPI_COMM_WORLD); 
  start = MPI_Wtime();
  MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  int globalRows = 500;
  int globalCols = 500;
  int rowsProcess = globalRows / numProcesses;
  int remainingRows = globalRows % numProcesses;
  float matrizArray[globalRows * globalCols];
  float localVectorA[globalCols * rowsProcess];
  float localVectorB[globalCols * rowsProcess];
  float localVectorC[globalCols * rowsProcess];
  float globalVectorC[globalRows * globalCols];
  create_zero_matrix(globalVectorC, globalRows, globalCols);
  if (rank == 0){
    create_matrix(matrizArray, globalRows, globalCols);
    }
  MPI_Scatter(matrizArray, globalCols*rowsProcess, MPI_FLOAT, localVectorB,
        globalCols*rowsProcess, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Scatter(globalVectorC, globalCols*rowsProcess, MPI_FLOAT, localVectorC,
        globalCols*rowsProcess, MPI_FLOAT, 0, MPI_COMM_WORLD);
  MPI_Scatter(matrizArray, globalCols*rowsProcess, MPI_FLOAT, localVectorA,
        globalCols*rowsProcess, MPI_FLOAT, 0, MPI_COMM_WORLD);
  for (int i = 0; i < globalCols * rowsProcess; i++){
    localVectorC[i] += localVectorA[i] * localVectorB[i];

    MPI_Gather(localVectorC, globalCols * rowsProcess, MPI_FLOAT, globalVectorC, 
        globalCols * rowsProcess, MPI_FLOAT, 0, MPI_COMM_WORLD);
    
  }
  if (rank == 0){
    for (int i = globalCols * rowsProcess * numProcesses; 
        i < (globalCols * rowsProcess * numProcesses) + (remainingRows * globalCols); i++){
      globalVectorC[i] += matrizArray[i] * matrizArray[i];
    } 
  }
  MPI_Barrier(MPI_COMM_WORLD); 
  end = MPI_Wtime();
  MPI_Finalize();
  if (rank == 0) { /* use time on master node */
    printf("Runtime = %f\n", end-start);
}

  return 0;

}
