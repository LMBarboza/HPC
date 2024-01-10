#include <iostream>
#include "../include/matrix.h"
#include "../include/matrixAVX.h"
#include "../include/matrixMP.h"
#include <vector>
#include <chrono>

std::vector<std::vector<float>> create_matrix(int rows, int cols) {
    std::vector<std::vector<float>> data(rows, std::vector<float>(cols, 0.0f));

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            data[i][j] = static_cast<float>(rand()) / RAND_MAX;
        }
    }
    return data;
}

template <class T>
void benchmark_operation(T& matA, T& matB, const std::string& operation, float a) {
    std::chrono::steady_clock::time_point begin = std::chrono::steady_clock::now();

    if (operation == "+") {
        T result = matA + matB;
        std::cout << "Benchmark Matrix Addition ";
    } else if (operation == "*") {
        T result = matA * a; 
        std::cout << "Benchmark Multiplication by a scalar ";
    } else if (operation == "/") {
        T result = matA / matB; 
        std::cout << "Benchmark Matrix Multiplication ";
    }

    std::chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

    std::cout << "Time = "
              << std::chrono::duration_cast<std::chrono::microseconds>(end - begin).count() / 100 << "[µs]"
              << std::endl;
}
int main(void) {
    const size_t rows = 1000;
    const size_t cols = 1000;
    std::cout << "Matrix Size:\n";
    std::cout << rows << "x" << cols << "\n";
    std::vector<std::vector<float>> mat1 = create_matrix(rows, cols);
    std::vector<std::vector<float>> mat2 = create_matrix(rows, cols);
    float a = 2.0;
    Matriz matriz1(mat1);
    Matriz matriz2(mat2);
    MatrizAVX matriz3(mat1);
    MatrizAVX matriz4(mat2);
    MatrizMP matriz5(mat1);
    MatrizMP matriz6(mat2);
    
    std::cout << "Basic: \n";
    benchmark_operation(matriz1, matriz2, "+", a);
    std::cout << "AVX: \n";
    benchmark_operation(matriz3, matriz4, "+", a);
    std::cout << "OpenMP: \n";
    benchmark_operation(matriz5, matriz6, "+", a);

    std::cout << "Basic: \n";
    benchmark_operation(matriz1, matriz2, "*", a);
    std::cout << "AVX: \n";
    benchmark_operation(matriz3, matriz4, "*", a);
    std::cout << "OpenMP: \n";
    benchmark_operation(matriz5, matriz6, "*", a);

    std::cout << "Basic: \n";
    benchmark_operation(matriz1, matriz2, "/", a);
    std::cout << "AVX: \n"; 
    benchmark_operation(matriz3, matriz4, "/", a);
    std::cout << "OpenMP: \n";
    benchmark_operation(matriz5, matriz6, "/", a);

    return 0;
}
