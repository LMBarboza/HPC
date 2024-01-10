#ifndef MATRIXMP_H
#define MATRIXMP_H

#include <vector>
#include <stdexcept>

class MatrizMP {
public:
    std::vector<std::vector<float>> elementos;

    MatrizMP(const std::vector<std::vector<float>>& dados);

    void print();

    MatrizMP operator+(const MatrizMP& obj);

    MatrizMP operator-(const MatrizMP& obj);

    MatrizMP operator*(float a);

    MatrizMP operator/(const MatrizMP& obj);

    MatrizMP transpor();

};
#endif
