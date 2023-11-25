#include <iostream>
#include <omp.h>
#include <stdio.h>
#include "matrix.h"


int main(){
    printf("Hello world\n");
    typedef short value_t;
    int dim_A_x = 40000, dim_A_y = 50000, dim_B_x = 40000, dim_B_y = 50000;
    // int dim_A_x = 4, dim_A_y = 5, dim_B_x = 4, dim_B_y = 5;
    value_t initial_A = 1, initial_B = 2;
    int threads = 2;

     // Matrix A
    size_t size_A = dim_A_x * dim_A_y * sizeof(value_t);

    // Matrix B
    size_t size_B = dim_B_x * dim_B_y * sizeof(value_t);

    // Initial size of Matrix C (before resizing in matrix<value_t>)
    size_t size_C_initial = dim_A_x * dim_A_y * sizeof(value_t);

    std::cout << "Estimated RAM usage:\n";
    std::cout << "Matrix A: " << size_A / (1024 * 1024) << " MB\n";
    std::cout << "Matrix B: " << size_B / (1024 * 1024) << " MB\n";
    std::cout << "Matrix C (initial): " << size_C_initial / (1024 * 1024) << " MB\n";

    // What will matrix do? dunno
    std::vector<std::vector<value_t>> A = std::vector<std::vector<value_t>>(dim_A_x,std::vector<value_t>(dim_A_y,initial_A));
    std::vector<std::vector<value_t>> B = std::vector<std::vector<value_t>>(dim_B_x,std::vector<value_t>(dim_B_y,initial_B));
    std::vector<std::vector<value_t>> C;
    
    matrix<value_t> obj = matrix<value_t>(&A,&B,&C, threads);

    obj.subtraction();

    return 0;
}