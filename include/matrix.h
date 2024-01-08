#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <stdexcept>

class Matriz {
public:
    std::vector<std::vector<float>> elementos;

    Matriz(const std::vector<std::vector<float>>& dados);

    void print();

    Matriz operator+(const Matriz& obj);

    Matriz operator-(const Matriz& obj);

    Matriz operator*(float a);

    Matriz operator/(const Matriz& obj);

    Matriz transpor();

};
#endif 
