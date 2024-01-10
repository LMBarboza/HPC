#ifndef MATRIXMP_H
#define MATRIXMP_H

#include <vector>
#include <stdexcept>

class MatrizMP {
public:
    std::vector<std::vector<float>> elementos;

    MatrizMP(std::vector<std::vector<float>>& dados);

    void print();

    MatrizMP operator+(MatrizMP& obj);

    MatrizMP operator-(MatrizMP& obj);

    MatrizMP operator*(float a);

    MatrizMP operator/(MatrizMP& obj);

    MatrizMP transpor();

};
#endif
