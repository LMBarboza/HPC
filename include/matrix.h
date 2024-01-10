#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <stdexcept>

class Matriz {
public:
    std::vector<std::vector<float>> elementos;

    Matriz(std::vector<std::vector<float>>& dados);

    void print();

    Matriz operator+(Matriz& obj);

    Matriz operator-(Matriz& obj);

    Matriz operator*(float a);

    Matriz operator/(Matriz& obj);

    Matriz transpor();

};
#endif 
