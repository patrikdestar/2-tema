#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include <omp.h>
#include <stdio.h>
#include <ctime>
#include <string>
#include <chrono>



/**
 * This class is what we're supposed to work with
*/
template <typename T>
class matrix
{
private:
    typedef std::vector<std::vector<T>> vec;
    vec* matrix_A; // pirma matrica
    vec* matrix_B; // antra matrica
    vec* matrix_C; // rezultatas
    long dim_A_x, dim_A_y, dim_B_x, dim_B_y; // dimensions A[x][y]
    int threads;
    std::chrono::duration<double> duration; // laikas

    /* Gets stuff */
    void get_stuff(vec* A, vec* B, vec* C){
        matrix_A = A;
        matrix_B = B;
        matrix_C = C;
        dim_A_x = matrix_A->size();
        dim_A_y = matrix_A->at(0).size();
        dim_B_x = matrix_B->size();
        dim_B_y = matrix_B->at(0).size();
        omp_set_num_threads(threads);
    }
    // ---------------------------------------------------------------------------------------------
public:
    /** Contructor
     * \param A A matrica
     * \param B B matrica
     * \param C C matrica (not initialised)
    */
    matrix(vec* A, vec* B, vec* C)
    {
        threads = omp_get_max_threads();
        get_stuff(A,B,C);
    }
    // ---------------------------------------------------------------------------------------------
    /** Constructor
     * \param A A matrica
     * \param B B matrica
     * \param C C matrica (not initialised)
     * \param thr threads count (default max)
    */
    matrix(vec* A, vec* B, vec* C,
        int thr)
    {
        threads = thr;
        get_stuff(A,B,C);
    }
    // ---------------------------------------------------------------------------------------------
    /* Printer (cant print complex for sure)*/
    void print(std::string which_one, long xx, long yy)
    {
        if (which_one == "sum"){
            printf("\nMatrcų A ir B suma C:\n\n");
            for (int y = 0; y < yy; ++y) {
                for (int x = 0; x < xx; ++x) {
                    std::cout << matrix_C->at(x).at(y) << " ";
                }
                std::cout << std::endl;
            }
        } else if (which_one == "sub"){
            printf("\nMatrcų A ir B skirtumas C:\n\n");
            for (int y = 0; y < yy; ++y) {
                for (int x = 0; x < xx; ++x) {
                    std::cout << matrix_C->at(x).at(y) << " ";
                }
                std::cout << std::endl;
            }
        } else if (which_one == "mul"){
            printf("\nMatrcų A ir B sandauga C:\n\n");
            for (int y = 0; y < yy; ++y) {
                for (int x = 0; x < xx; ++x) {
                    std::cout << matrix_C->at(x).at(y) << " ";
                }
                std::cout << std::endl;
            }
        }
        std::cout << std::endl;
    }
    // ---------------------------------------------------------------------------------------------
    /* Allocates */
    void allocate(){
        long xx = dim_B_x, yy = dim_A_y;
        matrix_C->resize(xx);
        for (long i = 0; i < xx; i++){
            matrix_C->at(i).resize(yy);
        }

    }
    // ---------------------------------------------------------------------------------------------
    /* Sudeda matricas */
    void summation(bool do_print = false){
        if (dim_A_x == dim_B_x && dim_A_y == dim_B_y){
            long xx = dim_A_x, yy = dim_A_y; // dimensions
            allocate();
            printf("Sumavimas prasideda:\n");
            auto start = std::chrono::high_resolution_clock::now();
            #pragma omp parallel for collapse(2)
            for (long x = 0; x < xx; ++x) {
                for (long y = 0; y < yy; ++y) {
                    matrix_C->at(x).at(y) = matrix_A->at(x).at(y) + matrix_B->at(x).at(y);
                }
            }
            duration = std::chrono::high_resolution_clock::now() - start;
            std::cout << "Elapsed time: " << duration.count() << " seconds\n";
            if (do_print) print("sum", xx, yy);
        } else {
            printf("\nNegalimas veiksmas, matricų dimensijos nesutampa!!!\n\n");
        }
    }
    // ---------------------------------------------------------------------------------------------
    /* Atima matricas */
    void subtraction(bool do_print = false){
        if (dim_A_x == dim_B_x && dim_A_y == dim_B_y){
            long xx = dim_A_x, yy = dim_A_y; // dimensions
            allocate();
            printf("Atemimas prasideda:\n");
            auto start = std::chrono::high_resolution_clock::now();
            #pragma omp parallel for collapse(2)
            for (long x = 0; x < xx; ++x) {
                for (long y = 0; y < yy; ++y) {
                    matrix_C->at(x).at(y) = matrix_A->at(x).at(y) - matrix_B->at(x).at(y);
                }
            }
            duration = std::chrono::high_resolution_clock::now() - start;
            std::cout << "Elapsed time: " << duration.count() << " seconds\n";
            if (do_print) print("sub", xx, yy);
        } else {
            printf("\nNegalimas veiksmas, matricų dimensijos nesutampa!!!\n\n");
        }
    }
    // ---------------------------------------------------------------------------------------------
    /* Daugina matricas*/
    void multiplication(bool do_print = false){
        if (dim_A_x == dim_B_y){
            allocate();
            printf("Daugyba praideda:\n");
            long dim = dim_A_x;
            auto start = std::chrono::high_resolution_clock::now();
            #pragma omp parallel for collapse(2)
            for (long yA = 0; yA < dim_A_y; ++yA) {
                for (long xB = 0; xB < dim_B_x; ++xB) {
                    T element = 0;
                    for (long n = 0; n < dim; n++) {
                        element = element + matrix_A->at(n).at(yA)*matrix_B->at(xB).at(n);
                    }
                    matrix_C->at(xB).at(yA) = element;
                }
            }
            duration = std::chrono::high_resolution_clock::now() - start;
            std::cout << "Elapsed time: " << duration.count() << " seconds\n";
            if (do_print) print("mul", dim_B_x, dim_A_y);
        } else {
            printf("\nNegalimas veiksmas, matricų dimensijos nesutampa!!!\n\n");
        }
    }
};

#endif