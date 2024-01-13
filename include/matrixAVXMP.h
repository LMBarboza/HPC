#ifndef MATRIXAVXMP_H 
#define MATRIXAVXMP_H
#include <vector>
#include <stdexcept>
#include <immintrin.h>

class MatrizAVXMP {
public:
    std::vector<std::vector<float>> elementos;

    MatrizAVXMP(std::vector<std::vector<float>>& dados);

    void print();

    MatrizAVXMP operator+(MatrizAVXMP& obj);

    MatrizAVXMP operator-(MatrizAVXMP& obj);

    MatrizAVXMP operator*(float a);

    MatrizAVXMP operator/(MatrizAVXMP& obj);

    MatrizAVXMP transpor();

};
#endif
