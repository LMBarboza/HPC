#ifndef MATRIXAVX_H
#define MATRIXAVX_H

#include <vector>
#include <stdexcept>
#include <immintrin.h>

class MatrizAVX {
public:
    std::vector<std::vector<float>> elementos;

    MatrizAVX(std::vector<std::vector<float>>& dados);

    void print();

    MatrizAVX operator+(MatrizAVX& obj);

    MatrizAVX operator-(MatrizAVX& obj);

    MatrizAVX operator*(float a);

    MatrizAVX operator/(MatrizAVX& obj);

    MatrizAVX transpor();

};
#endif
