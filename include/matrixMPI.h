#ifndef MATRIXMPI_H
#define MATRIXMPI_H

#include <vector>
#include <stdexcept>
#include <mpi.h>

class MatrizMPI {
public:
    std::vector<std::vector<float>> elementos;

    MatrizMPI(const std::vector<std::vector<float>>& dados);

    void print();

    MatrizMPI operator+(const MatrizMPI& obj);

    MatrizMPI operator-(const MatrizMPI& obj);

    MatrizMPI operator*(float a);

    MatrizMPI operator/(const MatrizMPI& obj);

    MatrizMPI transpor();

};
#endif
