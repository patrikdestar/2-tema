#include <iostream>
#include <omp.h>
#include <stdio.h>
#include "matrix.h"
#include "Complex.h"

int main(){
    printf("Hello world\n");
    typedef double value_t;
    // int dim_A_x = 4000, dim_A_y = 800, dim_B_x = 300, dim_B_y = 4000;
    int dim_A_x = 4, dim_A_y = 4, dim_B_x = 4, dim_B_y = 4;
    // value_t initial_A = value_t(1,1), initial_B = value_t(3, -1);
    value_t initial_A = 3.14, initial_B = 6.28;
    int threads = 8;

     // Matrix A
    size_t size_A = dim_A_x * dim_A_y * sizeof(value_t);

    // Matrix B
    size_t size_B = dim_B_x * dim_B_y * sizeof(value_t);

    // Initial size of Matrix C (before resizing in matrix<value_t>)
    size_t size_C_initial = dim_A_y * dim_B_x * sizeof(value_t);

    std::cout << "Estimated RAM usage:\n";
    std::cout << "Matrix A: " << size_A / (1024) << " kB\n";
    std::cout << "Matrix B: " << size_B / (1024) << " kB\n";
    std::cout << "Matrix C (initial): " << size_C_initial / (1024) << " kB\n";

    // What will matrix do? dunno
    std::vector<std::vector<value_t>> A = std::vector<std::vector<value_t>>(dim_A_x,std::vector<value_t>(dim_A_y,initial_A));
    std::vector<std::vector<value_t>> B = std::vector<std::vector<value_t>>(dim_B_x,std::vector<value_t>(dim_B_y,initial_B));
    std::vector<std::vector<value_t>> C = std::vector<std::vector<value_t>>(dim_B_x,std::vector<value_t>(dim_A_y,0));;
    
    /* Note: Matrix C has to be initialised before constructing object*/
    matrix<value_t> obj = matrix<value_t>(&A,&B,&C, threads);

    obj.multiplication(true);

    // // klasė veikia su compleksiniais skaiciais
    // for (int y = 0; y < dim_A_y; ++y) {
    //     for (int x = 0; x < dim_B_x; ++x) {
    //         //std::cout << C.size();
    //         C.at(x).at(y).print();
    //         std::cout << " ";
    //     }
    //     std::cout << std::endl;
    // }

    return 0;
}