#include "../include/matrixMP.h"
#include <iostream>
#include <omp.h>

MatrizMP::MatrizMP(const std::vector<std::vector<float>>& dados) : elementos(dados) {}

void MatrizMP::print() {
    for (const std::vector<float>& linha : elementos) {
        for (const float& elemento : linha) {
            std::cout << elemento << " ";
        }
        std::cout << "\n";
    }
}

MatrizMP MatrizMP::operator+(const MatrizMP& obj) {
    if (obj.elementos.size() != elementos.size() || obj.elementos[0].size() != elementos[0].size()) {
        throw std::invalid_argument("1");
    }

    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(elementos[0].size(), 0.0));
    #pragma omp parallel for 
    for (size_t i = 0; i < elementos.size(); ++i) {
        for (size_t j = 0; j < elementos[0].size(); ++j) {
            res[i][j] = elementos[i][j] + obj.elementos[i][j];
        }
    }

    return MatrizMP(res);
}

MatrizMP MatrizMP::operator-(const MatrizMP& obj) {
    if (obj.elementos.size() != elementos.size() || obj.elementos[0].size() != elementos[0].size()) {
        throw std::invalid_argument("1");
    }

    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(elementos[0].size(), 0.0));
    #pragma omp parallel for 
    for (size_t i = 0; i < elementos.size(); ++i) {
        for (size_t j = 0; j < elementos[0].size(); ++j) {
            res[i][j] = elementos[i][j] - obj.elementos[i][j];
        }
    }

    return MatrizMP(res);
}

MatrizMP MatrizMP::operator*(float a) {
    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(elementos[0].size(), 0.0));
    #pragma omp parallel for 
    for (size_t i = 0; i < elementos.size(); ++i) {
        for (size_t j = 0; j < elementos[0].size(); ++j) {
            res[i][j] = elementos[i][j] * a;
        }
    }

    return MatrizMP(res);
}

MatrizMP MatrizMP::operator/(const MatrizMP& obj) {
    if (elementos[0].size() != obj.elementos.size()) {
        throw std::invalid_argument("1");
    }

    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(obj.elementos[0].size(), 0.0));
    #pragma omp parallel for 
    for (size_t i = 0; i < elementos.size(); ++i) {
        for (size_t j = 0; j < obj.elementos[0].size(); ++j) {
            for (size_t k = 0; k < elementos[0].size(); ++k) {
                res[i][j] += elementos[i][k] * obj.elementos[k][j];
            }
        }
    }

    return MatrizMP(res);
}

MatrizMP MatrizMP::transpor() {
    std::vector<std::vector<float>> transposta(elementos[0].size(), std::vector<float>(elementos.size(), 0.0));
    #pragma omp parallel for
    for (size_t i = 0; i < elementos.size(); ++i) {
        for (size_t j = 0; j < elementos[0].size(); ++j) {
            transposta[j][i] = elementos[i][j];
        }
    }

    return MatrizMP(transposta);
}
