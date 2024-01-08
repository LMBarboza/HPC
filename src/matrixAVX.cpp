#include "../include/matrixAVX.h"
#include <iostream>
#include <immintrin.h>


MatrizAVX::MatrizAVX(const std::vector<std::vector<float>>& dados) : elementos(dados) {}

void MatrizAVX::print() {
    for (const std::vector<float>& linha : elementos) {
        for (const float& elemento : linha) {
            std::cout << elemento << " ";
        }
        std::cout << "\n";
    }
}

MatrizAVX MatrizAVX::operator+(const MatrizAVX& obj) {
    if (obj.elementos.size() != elementos.size() || obj.elementos[0].size() != elementos[0].size()) {
        throw std::invalid_argument("1");
    }

    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(elementos[0].size(), 0.0));

    for (size_t i = 0; i < elementos.size(); ++i) {
        const int vectorSamples = (elementos[0].size()/8)*8;
        int j = 0;
        for (; j < vectorSamples; j += 8) {
            __m256 aRegister = _mm256_loadu_ps(&elementos[i][j]);
            __m256 bRegister = _mm256_loadu_ps(&obj.elementos[i][j]);
            __m256 resRegister = _mm256_add_ps(aRegister, bRegister);
            _mm256_storeu_ps(&res[i][j], resRegister);

        }
        for (; j < elementos.size(); j++){
            res[i][j] = elementos[i][j] + obj.elementos[i][j];

        }
    }

    return MatrizAVX(res);
}

MatrizAVX MatrizAVX::operator-(const MatrizAVX& obj) {
    if (obj.elementos.size() != elementos.size() || obj.elementos[0].size() != elementos[0].size()) {
        throw std::invalid_argument("1");
    }

    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(elementos[0].size(), 0.0));
    for (size_t i = 0; i < elementos.size(); ++i) {
        const int vectorSamples = (elementos[0].size()/8)*8;
        int j = 0;
        for (; j < vectorSamples; j += 8) {
            __m256 aRegister = _mm256_loadu_ps(&elementos[i][j]);
            __m256 bRegister = _mm256_loadu_ps(&obj.elementos[i][j]);
            __m256 resRegister = _mm256_sub_ps(aRegister, bRegister);
            _mm256_storeu_ps(&res[i][j], resRegister);

        }
        for (; j < elementos.size(); j++){
            res[i][j] = elementos[i][j] - obj.elementos[i][j];

        }
    }

    return MatrizAVX(res);
    }

MatrizAVX MatrizAVX::operator*(float a) {
    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(elementos[0].size(), 0.0));
    __m256 scalarRegister = _mm256_set1_ps(a);
    for (size_t i = 0; i < elementos.size(); ++i) {
        const int vectorSamples = (elementos[0].size()/8)*8;
        int j = 0;
        for (; j < vectorSamples; j += 8) {
            __m256 aRegister = _mm256_loadu_ps(&elementos[i][j]);
            __m256 resRegister = _mm256_mul_ps(aRegister, scalarRegister);
            _mm256_storeu_ps(&res[i][j], resRegister);

        }
        for (; j < elementos.size(); j++){
            res[i][j] = elementos[i][j] * a;

        }
    }

    return MatrizAVX(res);
}


MatrizAVX MatrizAVX::operator/(const MatrizAVX& obj) {
    if (elementos.empty() || obj.elementos.empty() || elementos[0].size() != obj.elementos.size()) {
        throw std::invalid_argument("Invalid matrix dimensions");
    }

    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(obj.elementos[0].size(), 0.0));

    for (size_t i = 0; i < elementos.size(); ++i) {
        for (size_t j = 0; j < obj.elementos[0].size(); ++j) {
            __m256 resRegister = _mm256_setzero_ps();
            const size_t vectorSamples = (elementos[0].size() / 8) * 8;
            size_t k = 0;

            for (; k < vectorSamples; k += 8) {
                if (k < elementos[i].size() && k < obj.elementos.size()) {
                    __m256 aRegister = _mm256_loadu_ps(&elementos[i][k]);
                    __m256 bRegister = _mm256_loadu_ps(&obj.elementos[k][j]);

                    resRegister = _mm256_add_ps(resRegister, _mm256_mul_ps(aRegister, bRegister));
                }
            }

            // Corrected the syntax for storing the result
            if (i < res.size() && j < res[0].size()) {
                _mm256_storeu_ps(&res[i][j], resRegister);
            }

            for (; k < elementos[0].size(); ++k) {
                if (i < elementos.size() && k < elementos[i].size() && k < obj.elementos.size()) {
                    res[i][j] += elementos[i][k] * obj.elementos[k][j];
                }
            }
        }
    }

    return MatrizAVX(res);
}

/*MatrizAVX MatrizAVX::operator/(const MatrizAVX& obj) {
    if (elementos[0].size() != obj.elementos.size()) {
        throw std::invalid_argument("1");
    }

    std::vector<std::vector<float>> res(elementos.size(), std::vector<float>(obj.elementos[0].size(), 0.0));

    for (size_t i = 0; i < elementos.size(); ++i) {
        for (size_t j = 0; j < obj.elementos[0].size(); ++j) {
            __m256 resRegister = _mm256_setzero_ps();
            const size_t vectorSamples = (elementos[0].size() / 8) * 8;
            size_t k = 0;
            for (; k < vectorSamples; k += 8) {
                __m256 aRegister = _mm256_loadu_ps(&elementos[i][k]);
                __m256 bRegister = _mm256_loadu_ps(&obj.elementos[k][j]);

                resRegister = _mm256_add_ps(resRegister, _mm256_mul_ps(aRegister, bRegister));
              }
            _mm256_storeu_ps(&res[i][j], resRegister);

            for (; k < elementos[0].size(); k++){
                res[i][j] += elementos[i][k] * obj.elementos[k][j];

            }
        }
    }

    return MatrizAVX(res);
}
*/
MatrizAVX MatrizAVX::transpor() {
    std::vector<std::vector<float>> transposta(elementos[0].size(), std::vector<float>(elementos.size(), 0.0));

    for (size_t i = 0; i < elementos.size(); ++i) {
        for (size_t j = 0; j < elementos[0].size(); ++j) {
            transposta[j][i] = elementos[i][j];
        }
    }

    return MatrizAVX(transposta);
}
