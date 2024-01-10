#include "../include/matrixMPI.h"
#include <iostream>
#include <mpi.h>
#include <algorithm>
#define MASTER_RANK 0

#define MASTER_TAG 1
#define WORKER_TAG 2


MatrizMPI::MatrizMPI(std::vector<std::vector<float>>& dados) : elementos(dados) {}

void MatrizMPI::print() {
    for (std::vector<float>& linha : elementos) {
        for (float& elemento : linha) {
            std::cout << elemento << " ";
        }
        std::cout << "\n";
    }
}
MatrizMPI MatrizMPI::operator+(MatrizMPI& obj) {
    if (obj.elementos.size() != elementos.size() || obj.elementos[0].size() != elementos[0].size()) {
        throw std::invalid_argument("1");
    }
    MPI_Status status;
    int numRow, offset, i, j, processId, rowsPerProcess, remainingRows, numProcesses, rank;

    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(elementos[0].size(), 0.0));
    MPI_Init(NULL, NULL); 
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0){
      rowsPerProcess = elementos.size() / (numProcesses - 1);
      remainingRows = elementos.size() % (numProcesses - 1);
      offset = 0;
      for (processId = 1; processId < numProcesses; processId++){
        numRow = processId <= remainingRows ? rowsPerProcess + 1 : rowsPerProcess;
        MPI_Send(&offset, 1, MPI_INT, processId, 1, MPI_COMM_WORLD);
        MPI_Send(&numRow, 1, MPI_INT, processId, 1, MPI_COMM_WORLD);
        MPI_Send(&elementos[offset][0], numRow * elementos.size(), MPI_FLOAT, processId, 1, MPI_COMM_WORLD);
        MPI_Send(&obj.elementos, obj.elementos.size() * obj.elementos[0].size(), 
            MPI_FLOAT, processId, 1, MPI_COMM_WORLD);
        offset += numRow;
      }
      for (int processId = 1; processId < numProcesses; processId++){
        MPI_Recv(&offset, 1, MPI_INT, processId, 2, MPI_COMM_WORLD, &status);
        MPI_Recv(&numRow, 1, MPI_INT, processId, 2, MPI_COMM_WORLD, &status);
        MPI_Recv(&res[offset][0], numRow * elementos[0].size(), MPI_FLOAT, processId, 2, MPI_COMM_WORLD, &status);
      }
    }

    if (rank != 0){
      MPI_Recv(&offset, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		  MPI_Recv(&numRow, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		  MPI_Recv(&elementos, numRow * elementos[0].size(), MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);
		  MPI_Recv(&obj.elementos, obj.elementos.size() * obj.elementos[0].size(), MPI_FLOAT, 0, 1,
          MPI_COMM_WORLD, &status);

      for (i = 0; i < elementos.size(); ++i) {
        for (j = 0; j < numRow; ++j) {
            res[i][j] = elementos[i][j] + obj.elementos[i][j];
        }
      }
      MPI_Send(&offset, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
		  MPI_Send(&numRow, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
	  	MPI_Send(&res, numRow * elementos[0].size(), MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return MatrizMPI(res);
}
/*
MatrizMPI MatrizMPI::operator+(MatrizMPI& obj) {
    if (obj.elementos.size() != elementos.size() || obj.elementos[0].size() != elementos[0].size()) {
        throw std::invalid_argument("1");
    }
    MPI_Status status;
    int numRow, offset, i, j;
    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(elementos[0].size(), 0.0));
    MPI_Init(NULL, NULL);
    int numProcesses, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0){
      int rowsPerProcess = elementos.size() / (numProcesses - 1);
      int remainingRows = elementos.size() % (numProcesses - 1);
      offset = 0;
      for (int processId = 1; processId < numProcesses; processId++){
        numRow = processId <= remainingRows ? rowsPerProcess + 1 : rowsPerProcess;
        MPI_Send(&offset, 1, MPI_INT, processId, 1, MPI_COMM_WORLD);
        MPI_Send(&numRow, 1, MPI_INT, processId, 1, MPI_COMM_WORLD);
        MPI_Send(&elementos[offset][0], numRow * elementos[0].size(), MPI_FLOAT, processId, 1, MPI_COMM_WORLD);
        MPI_Send(&obj.elementos, obj.elementos.size() * obj.elementos[0].size(), 
            MPI_FLOAT, processId, 1, MPI_COMM_WORLD);
        offset += numRow;
      }
      for (int processId = 1; processId < numProcesses; processId++){
        MPI_Recv(&offset, 1, MPI_INT, processId, 2, MPI_COMM_WORLD, &status);
        MPI_Recv(&numRow, 1, MPI_INT, processId, 2, MPI_COMM_WORLD, &status);
        MPI_Recv(&res[offset][0], numRow * elementos[0].size(), MPI_FLOAT, processId, 2, MPI_COMM_WORLD, &status);
      }
    }

    if (rank != 0){
      MPI_Recv(&offset, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		  MPI_Recv(&numRow, 1, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		  MPI_Recv(&elementos, numRow * elementos[0].size(), MPI_FLOAT, 0, 1, MPI_COMM_WORLD, &status);
		  MPI_Recv(&obj.elementos, obj.elementos.size() * obj.elementos[0].size(), MPI_FLOAT, 0, 1,
          MPI_COMM_WORLD, &status);

      for (i = 0; i < elementos.size(); ++i) {
        for (j = 0; j < numRow; ++j) {
            res[i][j] = elementos[i][j] + obj.elementos[i][j];
        }
      }
      MPI_Send(&offset, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
		  MPI_Send(&numRow, 1, MPI_INT, 0, 2, MPI_COMM_WORLD);
	  	MPI_Send(&res, numRow * elementos[0].size(), MPI_FLOAT, 0, 2, MPI_COMM_WORLD);
    }
    MPI_Finalize();
    return MatrizMPI(res);
}*/
/*MatrizMPI MatrizMPI::operator+(const MatrizMPI& obj) {
    if (obj.elementos.size() != elementos.size() || obj.elementos[0].size() != elementos[0].size()) {
        throw std::invalid_argument("1");
    }

    MPI_Init(NULL, NULL);
    int numProcesses, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(elementos[0].size(), 0.0f));
    std::vector<std::vector<float>> resGather(elementos.size(), std::vector<float>(elementos[0].size(), 0.0f));

    int rowsPerProcess = elementos.size() / numProcesses;
    int remainingRows = elementos.size() % numProcesses;

    int startIndex = rank * rowsPerProcess + (rank < remainingRows ? rank : remainingRows);

    int endIndex = startIndex + rowsPerProcess + (rank < remainingRows ? 1 : 0);

    std::vector<std::vector<float>> localMatrixA(rowsPerProcess, std::vector<float>(elementos[0].size(), 0.0f));
    std::vector<std::vector<float>> localMatrixB(rowsPerProcess, std::vector<float>(elementos[0].size(), 0.0f));

    MPI_Scatter(elementos.data(), rowsPerProcess * elementos[0].size(), MPI_FLOAT,
        localMatrixA.data(), rowsPerProcess * elementos[0].size(), MPI_FLOAT, 0, MPI_COMM_WORLD);
    MPI_Scatter(obj.elementos.data(), rowsPerProcess * elementos[0].size(), MPI_FLOAT,
        localMatrixB.data(), rowsPerProcess * elementos[0].size(), MPI_FLOAT, 0, MPI_COMM_WORLD);
    
    for (size_t i = 0; i < localMatrixA.size(); ++i) {
        for (size_t j = 0; j < localMatrixA[0].size(); ++j) {
            res[i][j] = localMatrixA[i][j] + localMatrixB[i][j];
        }
    }

    MPI_Gather(res.data(), rowsPerProcess * elementos[0].size(), MPI_FLOAT,
        resGather.data(), rowsPerProcess * elementos[0].size(), MPI_FLOAT, 0, MPI_COMM_WORLD);
    


    MPI_Finalize();
    return MatrizMPI(res);
}
*/
MatrizMPI MatrizMPI::operator-(MatrizMPI& obj) {
    if (obj.elementos.size() != elementos.size() || obj.elementos[0].size() != elementos[0].size()) {
        throw std::invalid_argument("1");
    }

    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(elementos[0].size(), 0.0));

    for (size_t i = 0; i < elementos.size(); ++i) {
        for (size_t j = 0; j < elementos[0].size(); ++j) {
            res[i][j] = elementos[i][j] - obj.elementos[i][j];
        }
    }

    return MatrizMPI(res);
}

MatrizMPI MatrizMPI::operator*(float a) {
    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(elementos[0].size(), 0.0));

    for (size_t i = 0; i < elementos.size(); ++i) {
        for (size_t j = 0; j < elementos[0].size(); ++j) {
            res[i][j] = elementos[i][j] * a;
        }
    }

    return MatrizMPI(res);
}

MatrizMPI MatrizMPI::operator/(MatrizMPI& obj){
  std::vector<std::vector<float>> a = elementos;
  std::vector<std::vector<float>> b = obj.elementos;
  int M_SIZE = elementos.size();
  std::vector<std::vector<float>> c(M_SIZE, std::vector<float>(M_SIZE)); 
  int communicator_size;
	int process_rank;
	int process_id;
	int offset;
	int rows_num;
	int workers_num;
	int remainder;
	int whole_part;
	int message_tag;
	int i;
	int j;
	int k;
  MPI_Status status;
  
	MPI_Init(NULL, NULL);
  
	MPI_Comm_size(MPI_COMM_WORLD, &communicator_size);
	MPI_Comm_rank(MPI_COMM_WORLD, &process_rank);
  if(process_rank == MASTER_RANK){
		workers_num = communicator_size - 1;
		whole_part = elementos.size() / workers_num;
		remainder = elementos.size() % workers_num;
		offset = 0;

		message_tag = MASTER_TAG;
		for(process_id = 1; process_id <= workers_num; process_id ++ ){
			rows_num = process_id <= remainder ? whole_part + 1 : whole_part;
			MPI_Send(&offset, 1, MPI_INT, process_id, message_tag, MPI_COMM_WORLD);
			MPI_Send(&rows_num, 1, MPI_INT, process_id, message_tag, MPI_COMM_WORLD);
		  MPI_Send(&a[offset][0], rows_num * M_SIZE, MPI_FLOAT, 
          process_id, message_tag, MPI_COMM_WORLD);
			MPI_Send(&b,M_SIZE*M_SIZE, MPI_FLOAT,
          process_id, message_tag, MPI_COMM_WORLD);

			offset += rows_num;
		}

		message_tag = WORKER_TAG;
		for(process_id = 1; process_id <= workers_num; process_id ++){
			MPI_Recv(&offset, 1, MPI_INT, process_id, message_tag, MPI_COMM_WORLD, &status);
			MPI_Recv(&rows_num, 1, MPI_INT, process_id, message_tag, MPI_COMM_WORLD, &status);
		  MPI_Recv(&c[offset][0], rows_num * M_SIZE, MPI_FLOAT,
          process_id, message_tag, MPI_COMM_WORLD, &status);
		}
			} 

if(process_rank != MASTER_RANK){
		message_tag = MASTER_TAG;
		MPI_Recv(&offset, 1, MPI_INT, MASTER_RANK, message_tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&rows_num, 1, MPI_INT, MASTER_RANK, message_tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&a, rows_num * M_SIZE, MPI_FLOAT,
        MASTER_RANK, message_tag, MPI_COMM_WORLD, &status);
		MPI_Recv(&b, M_SIZE * M_SIZE, MPI_FLOAT, MASTER_RANK,
        message_tag, MPI_COMM_WORLD, &status);

		for(k = 0; k < M_SIZE; k ++){
			for(i = 0; i < rows_num; i ++){
        c[i][k] = 0;
				for(j = 0; j < M_SIZE; j ++){
					c[i][k] += a[i][j] * b[j][k];
				}
			}
		}

		message_tag = WORKER_TAG;
		MPI_Send(&offset, 1, MPI_INT, MASTER_RANK, message_tag, MPI_COMM_WORLD);
		MPI_Send(&rows_num, 1, MPI_INT, MASTER_RANK, message_tag, MPI_COMM_WORLD);
		MPI_Send(&c, rows_num * M_SIZE, MPI_FLOAT, MASTER_RANK, message_tag, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return MatrizMPI(c);

}


/*MatrizMPI MatrizMPI::operator/(MatrizMPI& obj) {
    if (elementos[0].size() != obj.elementos.size()) {
        throw std::invalid_argument("1");
    }

    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(obj.elementos[0].size(), 0.0));

    for (size_t i = 0; i < elementos.size(); ++i) {
        for (size_t j = 0; j < obj.elementos[0].size(); ++j) {
            for (size_t k = 0; k < elementos[0].size(); ++k) {
                res[i][j] += elementos[i][k] * obj.elementos[k][j];
            }
        }
    }

    return MatrizMPI(res);
}
*/
MatrizMPI MatrizMPI::transpor() {
    std::vector<std::vector<float>> transposta(elementos[0].size(), std::vector<float>(elementos.size(), 0.0));

    for (size_t i = 0; i < elementos.size(); ++i) {
        for (size_t j = 0; j < elementos[0].size(); ++j) {
            transposta[j][i] = elementos[i][j];
        }
    }

    return MatrizMPI(transposta);
}
