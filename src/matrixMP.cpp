#include "../include/matrixMP.h"
#include <iostream>
#include <omp.h>

MatrizMP::MatrizMP(std::vector<std::vector<float>>& dados) : elementos(dados) {}

void MatrizMP::print() {
    for (std::vector<float>& linha : elementos) {
        for (float& elemento : linha) {
            std::cout << elemento << " ";
        }
        std::cout << "\n";
    }
}

MatrizMP MatrizMP::operator+(MatrizMP& obj) {
    if (obj.elementos.size() != elementos.size() || obj.elementos[0].size() != elementos[0].size()) {
        throw std::invalid_argument("1");
    }
    size_t i, j;
    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(elementos[0].size(), 0.0));
    #pragma omp parallel for private(i, j) 
    for (i = 0; i < elementos.size(); ++i) {
        for (j = 0; j < elementos[0].size(); ++j) {
            res[i][j] = elementos[i][j] + obj.elementos[i][j];
        }
    }

    return MatrizMP(res);
}

MatrizMP MatrizMP::operator-(MatrizMP& obj) {
    if (obj.elementos.size() != elementos.size() || obj.elementos[0].size() != elementos[0].size()) {
        throw std::invalid_argument("1");
    }
    size_t i, j;
    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(elementos[0].size(), 0.0));
    #pragma omp parallel for private(i, j)
    for (i = 0; i < elementos.size(); ++i) {
        for (j = 0; j < elementos[0].size(); ++j) {
            res[i][j] = elementos[i][j] - obj.elementos[i][j];
        }
    }

    return MatrizMP(res);
}

MatrizMP MatrizMP::operator*(float a) {
    size_t i, j;
    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(elementos[0].size(), 0.0));
    #pragma omp parallel for private (i, j)
    for (i = 0; i < elementos.size(); ++i) {
        for (j = 0; j < elementos[0].size(); ++j) {
            res[i][j] = elementos[i][j] * a;
        }
    }

    return MatrizMP(res);
}

MatrizMP MatrizMP::operator/(MatrizMP& obj) {
    if (elementos[0].size() != obj.elementos.size()) {
        throw std::invalid_argument("1");
    }
    size_t i, j, k;
    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(obj.elementos[0].size(), 0.0));
    #pragma omp parallel for private(i, j, k)
    for (i = 0; i < elementos.size(); ++i) {
        for (j = 0; j < obj.elementos[0].size(); ++j) {
            for (k = 0; k < elementos[0].size(); ++k) {
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
