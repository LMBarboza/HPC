#include "../include/matrixMPI.h"
#include <iostream>
#include <mpi.h>
MatrizMPI::MatrizMPI(const std::vector<std::vector<float>>& dados) : elementos(dados) {}

void MatrizMPI::print() {
    for (const std::vector<float>& linha : elementos) {
        for (const float& elemento : linha) {
            std::cout << elemento << " ";
        }
        std::cout << "\n";
    }
}

MatrizMPI MatrizMPI::operator+(const MatrizMPI& obj) {
    if (obj.elementos.size() != elementos.size() || obj.elementos[0].size() != elementos[0].size()) {
        throw std::invalid_argument("1");
    }

    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(elementos[0].size(), 0.0));
    MPI_Init(NULL, NULL);
    int numProcesses, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &numProcesses);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);


    for (size_t i = 0; i < elementos.size(); ++i) {
        for (size_t j = 0; j < elementos[0].size(); ++j) {
            res[i][j] = elementos[i][j] + obj.elementos[i][j];
        }
    }
    MPI_Finalize();
    return MatrizMPI(res);
}

MatrizMPI MatrizMPI::operator-(const MatrizMPI& obj) {
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

MatrizMPI MatrizMPI::operator/(const MatrizMPI& obj) {
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

MatrizMPI MatrizMPI::transpor() {
    std::vector<std::vector<float>> transposta(elementos[0].size(), std::vector<float>(elementos.size(), 0.0));

    for (size_t i = 0; i < elementos.size(); ++i) {
        for (size_t j = 0; j < elementos[0].size(); ++j) {
            transposta[j][i] = elementos[i][j];
        }
    }

    return MatrizMPI(transposta);
}
